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

use xscalenet::data::Dataset;
use xscalenet::discretize;
use xscalenet::output;
use xscalenet::scalenet::{self, Method, Params};
use xscalenet::scs::{self, MethodInfo};

struct Args {
    input: String,
    outdir: String,
    methods: Vec<Method>,
    params: Params,
    discretize: bool,
    transpose: bool,
    edges_only: bool,
    max_clusters: usize,
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
    })
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

fn info_for(method: Method) -> MethodInfo {
    match method {
        // aracne: no orientation, ranked by presence score
        Method::Aracne => MethodInfo {
            method,
            ort: false,
            weight_ecorr: false,
        },
        // bayes_hc: oriented, ranked by correlation magnitude
        Method::BayesHc => MethodInfo {
            method,
            ort: true,
            weight_ecorr: true,
        },
    }
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

    // Shared preprocessing (spectral embedding + subgraphs).
    let prep = scalenet::prepare(&data, &args.params)?;
    eprintln!(
        "# Kept k={} eigenvectors, {} subgraphs, m={} vars/subgraph",
        prep.spectral.k,
        prep.subgraphs.len(),
        prep.subset_m
    );

    // Precompute per-column ranks once (reused for every edge's ecorr).
    let ranks = data.rank_matrix();

    // Run each method, keeping per-threshold edge lists in memory.
    // methods_lists[method_i][thresh_i]
    let mut methods_lists = Vec::new();
    for &method in &args.methods {
        eprintln!("# Reconstructing with {}...", method.name());
        let mut lists = scalenet::run_method(&data, &prep, method, &args.params);
        for el in &mut lists {
            scalenet::fill_ecorr(&data, &ranks, el);
        }
        // Write per-method global edge lists.
        for (ti, thr) in args.params.pres_freq_thresh.iter().enumerate() {
            let path: PathBuf = [
                &args.outdir,
                "globalGraph",
                &format!("globalNet_presFreq{thr}"),
                &format!("edgesList.{}.txt", method.name()),
            ]
            .iter()
            .collect();
            output::write_edge_list(&path, &data, &lists[ti], args.edges_only)
                .map_err(|e| format!("write {}: {e}", path.display()))?;
        }
        methods_lists.push(lists);
    }

    // Consensus per threshold.
    for (ti, thr) in args.params.pres_freq_thresh.iter().enumerate() {
        let pairs: Vec<(MethodInfo, &scalenet::EdgeList)> = args
            .methods
            .iter()
            .enumerate()
            .map(|(mi, &m)| (info_for(m), &methods_lists[mi][ti]))
            .collect();
        let (cons, raw) = scs::consensus(&data, &pairs);

        let base: PathBuf = [
            &args.outdir,
            "consensusGraph",
            &format!("globalNet_presFreq{thr}"),
        ]
        .iter()
        .collect();
        output::write_edge_list(&base.join("edgesList.txt"), &data, &cons, args.edges_only)
            .map_err(|e| format!("write consensus: {e}"))?;
        output::write_raw_consensus(&base.join("edgesList.consensus.rawAvg.txt"), &data, &raw)
            .map_err(|e| format!("write consensus raw: {e}"))?;

        let n_edges = cons.len();
        eprintln!("# Consensus @ presFreq {thr}: {n_edges} edges");
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
