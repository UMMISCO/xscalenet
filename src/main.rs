//! Command-line entry point for xscalenet.
//!
//! Usage:
//!   xscalenet <input.tsv> <output_dir> [options]
//!
//! Input may be a plain TSV or gzip (`.gz` extension or gzip magic bytes).
//!
//! Options (defaults match the README `scs` example):
//!   --methods aracne,bayes_hc     reconstruction methods to embed
//!   --var-perc 0.2                fraction of variables per subgraph
//!   --eigen-perc auto             fraction of eigenvectors (auto = elbow)
//!   --thresh 0.3,0.8              presence-frequency thresholds
//!   --epsilon 0.001               ARACNE DPI tolerance
//!   --restart 21                  hill-climbing random restarts
//!   --seed 6196                   PRNG seed
//!   --discretize                  PAM-discretize numeric columns first
//!   --transpose                   input has features in rows (feature id in col 1)
//!   --edges-only                  write only present edges (skip absent pairs)
//!   --max-clusters 5              max PAM clusters per discretized column
//!   --maxp N                      max parents per node in hill-climbing (default: unlimited)
//!   --max-iter 10000              hill-climbing iteration cap
//!   --iss 1.0                     BDeu imaginary sample size

use std::path::PathBuf;
use std::process::exit;

use xscalenet::atc;
use xscalenet::data::Dataset;
use xscalenet::discretize;
use xscalenet::metrics;
use xscalenet::output;
use xscalenet::pipeline::RunConfig;
use xscalenet::scalenet::{Method, Params};

struct Args {
    input: String,
    outdir: String,
    methods: Vec<Method>,
    params: Params,
    discretize: bool,
    transpose: bool,
    edges_only: bool,
    max_clusters: usize,
    metrics: bool,
    metrics_weighted: bool,
    atc_file: Option<String>,
    atc_col: Option<String>,
}

fn parse_args() -> Result<Args, String> {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 3 {
        return Err(format!(
            "usage: {} <input(.tsv|.gz)> <output_dir> [--methods aracne,bayes_hc] \
             [--var-perc 0.2] [--eigen-perc auto] [--thresh 0.3,0.8] \
             [--epsilon 0.001] [--restart 21] [--seed 6196] [--discretize] \
             [--transpose] [--edges-only] [--max-clusters 5] [--maxp N] \
             [--max-iter 10000] [--iss 1.0]",
            argv.first().map(|s| s.as_str()).unwrap_or("xscalenet")
        ));
    }
    let input = argv[1].clone();
    let outdir = argv[2].clone();

    let mut methods = vec![Method::Aracne, Method::BayesHc];
    let mut discretize = false;
    let mut transpose = false;
    let mut edges_only = false;
    let mut max_clusters = 5usize;
    let mut metrics = false;
    let mut metrics_weighted = false;
    let mut atc_file: Option<String> = None;
    let mut atc_col: Option<String> = None;
    let mut params = Params {
        var_perc: 0.2,
        eigen_perc: None,
        pres_freq_thresh: vec![0.3, 0.8],
        epsilon: 0.001,
        restart: 21,
        seed: 6196,
        max_iter: 10_000,
        iss: 1.0,
        maxp: usize::MAX,
    };

    let mut i = 3;
    while i < argv.len() {
        let a = argv[i].clone();
        // Fetch the value for options that take one.
        let mut value = || -> Result<String, String> {
            i += 1;
            argv.get(i)
                .cloned()
                .ok_or_else(|| format!("missing value for {a}"))
        };
        match a.as_str() {
            "--discretize" => discretize = true,
            "--transpose" => transpose = true,
            "--edges-only" => edges_only = true,
            "--max-clusters" => max_clusters = value()?.parse().map_err(|e| format!("{e}"))?,
            "--maxp" => params.maxp = value()?.parse().map_err(|e| format!("{e}"))?,
            "--max-iter" => params.max_iter = value()?.parse().map_err(|e| format!("{e}"))?,
            "--iss" => params.iss = value()?.parse().map_err(|e| format!("{e}"))?,
            "--metrics" => metrics = true,
            "--metrics-weighted" => {
                metrics = true;
                metrics_weighted = true;
            }
            "--atc" => atc_file = Some(value()?),
            "--atc-col" => atc_col = Some(value()?),
            "--methods" => {
                methods = value()?
                    .split(',')
                    .map(|m| match m.trim() {
                        "aracne" => Ok(Method::Aracne),
                        "bayes_hc" => Ok(Method::BayesHc),
                        other => Err(format!("unknown method: {other}")),
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "--var-perc" => params.var_perc = value()?.parse().map_err(|e| format!("{e}"))?,
            "--eigen-perc" => {
                let v = value()?;
                params.eigen_perc = if v == "auto" || v == "-1" {
                    None
                } else {
                    Some(v.parse().map_err(|e| format!("{e}"))?)
                };
            }
            "--thresh" => {
                params.pres_freq_thresh = value()?
                    .split(',')
                    .map(|t| t.trim().parse::<f64>().map_err(|e| format!("{e}")))
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "--epsilon" => params.epsilon = value()?.parse().map_err(|e| format!("{e}"))?,
            "--restart" => params.restart = value()?.parse().map_err(|e| format!("{e}"))?,
            "--seed" => params.seed = value()?.parse().map_err(|e| format!("{e}"))?,
            other => return Err(format!("unknown option: {other}")),
        }
        i += 1;
    }

    Ok(Args {
        input,
        outdir,
        methods,
        params,
        discretize,
        transpose,
        edges_only,
        max_clusters,
        metrics,
        metrics_weighted,
        atc_file,
        atc_col,
    })
}

/// Read a node-annotation TSV (header, node id in column 1) and map each network
/// node to a categorical theme index for a chosen column (default: column 2).
/// Returns per-node `Some(theme)`/`None` plus the distinct theme names.
fn read_annotation_labels(
    path: &str,
    node_names: &[String],
    col: Option<&str>,
) -> Result<(Vec<Option<usize>>, Vec<String>), String> {
    let text = read_input(path)?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or("empty annotation file")?
        .split('\t')
        .map(str::trim)
        .collect();
    let col_idx = match col {
        Some(name) => header
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| format!("annotation column '{name}' not found"))?,
        None => 1.min(header.len().saturating_sub(1)),
    };
    // node name -> raw annotation value
    let mut value_of: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        if f.len() > col_idx {
            value_of.insert(f[0].to_string(), f[col_idx].to_string());
        }
    }
    // distinct values -> theme indices (sorted for determinism)
    let mut themes: Vec<String> = value_of.values().cloned().collect();
    themes.sort();
    themes.dedup();
    let theme_idx: std::collections::HashMap<&String, usize> =
        themes.iter().enumerate().map(|(i, t)| (t, i)).collect();
    let labels = node_names
        .iter()
        .map(|n| value_of.get(n).and_then(|v| theme_idx.get(v).copied()))
        .collect();
    Ok((labels, themes))
}

/// Read a TSV, transparently decompressing gzip (`.gz` extension or gzip magic).
fn read_input(path: &str) -> Result<String, String> {
    use std::io::Read;
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    if path.ends_with(".gz") || bytes.starts_with(&[0x1f, 0x8b]) {
        let mut s = String::new();
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(|e| format!("cannot decompress {path}: {e}"))?;
        Ok(s)
    } else {
        String::from_utf8(bytes).map_err(|e| format!("{path} is not valid UTF-8: {e}"))
    }
}

/// Transpose a features-in-rows matrix (first column = feature id, header =
/// observation ids) into the features-in-columns layout xscalenet expects.
fn transpose_tsv(text: &str) -> Result<String, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines.next().ok_or("empty input")?;
    let nobs = header.split('\t').count().saturating_sub(1); // drop the corner label
    let mut feat_names: Vec<&str> = Vec::new();
    let mut rows: Vec<Vec<&str>> = Vec::new();
    for (r, line) in lines.enumerate() {
        let mut it = line.split('\t');
        let fname = it.next().unwrap_or("");
        let vals: Vec<&str> = it.collect();
        if vals.len() != nobs {
            return Err(format!(
                "feature row {} has {} values, expected {nobs}",
                r + 1,
                vals.len()
            ));
        }
        feat_names.push(fname);
        rows.push(vals);
    }
    let mut out = String::with_capacity(text.len() + feat_names.len());
    out.push_str(&feat_names.join("\t"));
    out.push('\n');
    for j in 0..nobs {
        for (i, row) in rows.iter().enumerate() {
            if i > 0 {
                out.push('\t');
            }
            out.push_str(row[j]);
        }
        out.push('\n');
    }
    Ok(out)
}

fn run() -> Result<(), String> {
    let args = parse_args()?;

    // Load and (optionally) discretize. The file text and the zero-copy string
    // slices are scoped so they are freed before the (heavier) reconstruction.
    let data = {
        let mut text = read_input(&args.input)?;
        if args.transpose {
            text = transpose_tsv(&text)?;
        }
        let (names, raw, nrow) = Dataset::parse_tsv_raw(&text)?;
        if args.discretize {
            eprintln!("# Discretizing numeric columns (pamk)...");
            discretize::discretize_to_dataset(names, &raw, nrow, args.max_clusters)
        } else {
            Dataset::from_raw_columns(names, &raw, nrow)
        }
    };
    eprintln!(
        "# Loaded {} variables x {} observations",
        data.ncol(),
        data.nrow
    );

    // Reconstruct + fuse via the library pipeline; formatting is the CLI's job.
    let config = RunConfig {
        methods: args.methods.clone(),
        params: args.params.clone(),
    };
    let result = xscalenet::pipeline::run(&data, &config, |msg| eprintln!("{msg}"))?;

    // Per-method global edge lists.
    for (mi, method) in result.methods.iter().enumerate() {
        for (ti, thr) in result.thresholds.iter().enumerate() {
            let path: PathBuf = [
                &args.outdir,
                "globalGraph",
                &format!("globalNet_presFreq{thr}"),
                &format!("edgesList.{}.txt", method.name()),
            ]
            .iter()
            .collect();
            output::write_edge_list(&path, &result.nodes, &result.per_method[mi][ti], args.edges_only)
                .map_err(|e| format!("write {}: {e}", path.display()))?;
        }
    }

    // Consensus per threshold.
    for (ti, thr) in result.thresholds.iter().enumerate() {
        let base: PathBuf = [
            &args.outdir,
            "consensusGraph",
            &format!("globalNet_presFreq{thr}"),
        ]
        .iter()
        .collect();
        let cons = &result.consensus[ti];
        output::write_edge_list(&base.join("edgesList.txt"), &result.nodes, &cons.network, args.edges_only)
            .map_err(|e| format!("write consensus: {e}"))?;
        output::write_raw_consensus(
            &base.join("edgesList.consensus.rawAvg.txt"),
            &result.nodes,
            &cons.ranked,
        )
        .map_err(|e| format!("write consensus raw: {e}"))?;
    }

    // Optional per-node analysis (topological metrics and/or ATC) on each
    // consensus network.
    if args.metrics || args.atc_file.is_some() {
        let atc_labels = match &args.atc_file {
            Some(path) => Some(read_annotation_labels(
                path,
                &result.nodes,
                args.atc_col.as_deref(),
            )?),
            None => None,
        };
        for (ti, thr) in result.thresholds.iter().enumerate() {
            let net = &result.consensus[ti].network;
            let g = metrics::Graph::from_network(net, result.nodes.len(), args.metrics_weighted);
            let gm = metrics::compute(&g);
            let atc = atc_labels.as_ref().map(|(labels, themes)| {
                let t0 = atc::membership_from_labels(labels, themes.len());
                let cfg = atc::AtcConfig {
                    weighted: args.metrics_weighted,
                    ..Default::default()
                };
                atc::compute(&g, &t0, &cfg).atc
            });
            eprintln!(
                "# Analysis @ presFreq {thr}: {} nodes, {} edges, {} communities, modularity {:.3}{}",
                gm.n_nodes,
                gm.n_edges,
                gm.n_communities,
                gm.modularity,
                atc.as_ref()
                    .map(|_| format!(", ATC over {} themes", atc_labels.as_ref().unwrap().1.len()))
                    .unwrap_or_default(),
            );
            let path: PathBuf = [
                &args.outdir,
                "consensusGraph",
                &format!("globalNet_presFreq{thr}"),
                "nodeAnalysis.txt",
            ]
            .iter()
            .collect();
            output::write_node_analysis(&path, &result.nodes, &gm, atc.as_deref())
                .map_err(|e| format!("write node analysis: {e}"))?;
        }
    }

    eprintln!("# Done. Outputs under {}", args.outdir);
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        exit(1);
    }
}
