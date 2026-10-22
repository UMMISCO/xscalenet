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
//!   --metrics                     per-node topological metrics on the consensus network
//!   --atc annot.tsv               propagate a node annotation (ATC centrality)
//!   --atc-col NAME                annotation column to propagate (default: col 2)
//!   --theme-network               build the multi-layer theme-proximity network
//!   --json                        also export graph.json (nodes+metrics+edges+theme layer)
//!   --legacy-lrw                  reproduce the pre-fix, column-order-dependent spectral step
//!   --hc-orders 200               random variable orderings fused per subgraph (hill climbing)
//!   --hc-order-freq 0.0           drop arcs seen in fewer than this fraction of orderings
//!   (with --hc-orders > 1, --restart defaults to 0 unless given: see parse_args)

use std::path::PathBuf;
use std::process::exit;

/// Random variable orderings fused per subgraph by the Bayesian step.
///
/// This is a *robustness* knob, not the source of label-invariance: the base
/// order is already derived from the variables' own data (see
/// `scalenet::run_method`), so renaming or reordering the input gives an
/// identical network at any count. Fusing several orderings additionally stops
/// the result resting on one arbitrary scan order, at the price of ordinary
/// Monte-Carlo variation with the seed.
const HC_ORDERS_DEFAULT: usize = 200;

use xscalenet::analysis;
use xscalenet::data::Dataset;
use xscalenet::discretize;
use xscalenet::io;
use xscalenet::output;
use xscalenet::pipeline::RunConfig;
use xscalenet::scalefree::{self, Sign};
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
    theme_network: bool,
    json: bool,
    scale_free: bool,
}

fn parse_args() -> Result<Args, String> {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 3 {
        return Err(format!(
            "usage: {} <input(.tsv|.gz)> <output_dir> [--methods aracne,bayes_hc] \
             [--var-perc 0.2] [--eigen-perc auto] [--thresh 0.3,0.8] \
             [--epsilon 0.001] [--restart 21] [--seed 6196] [--discretize] \
             [--transpose] [--edges-only] [--max-clusters 5] [--maxp N] \
             [--max-iter 10000] [--iss 1.0] [--metrics] [--atc annot.tsv] \
             [--atc-col NAME] [--theme-network] [--json] [--scale-free] \
             [--iterative-eigen] [--legacy-lrw] [--hc-orders 200] \
             [--hc-order-freq 0.0]",
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
    let mut theme_network = false;
    let mut json = false;
    let mut scale_free = false;
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
        iterative_eigen: false,
        legacy_lrw: false,
        hc_orders: HC_ORDERS_DEFAULT,
        hc_order_freq: 0.0,
    };
    // Whether --restart was given explicitly; see the note below.
    let mut restart_set = false;

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
            "--theme-network" => theme_network = true,
            "--json" => json = true,
            "--scale-free" => scale_free = true,
            "--iterative-eigen" => params.iterative_eigen = true,
            "--legacy-lrw" => params.legacy_lrw = true,
            "--hc-orders" => params.hc_orders = value()?.parse().map_err(|e| format!("{e}"))?,
            "--hc-order-freq" => {
                params.hc_order_freq = value()?.parse().map_err(|e| format!("{e}"))?
            }
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
            "--restart" => {
                params.restart = value()?.parse().map_err(|e| format!("{e}"))?;
                restart_set = true;
            }
            "--seed" => params.seed = value()?.parse().map_err(|e| format!("{e}"))?,
            other => return Err(format!("unknown option: {other}")),
        }
        i += 1;
    }

    // Hill-climbing restarts and ordering consensus both exist to diversify the
    // greedy search, and under ordering consensus the restarts are redundant:
    // measured on the reference matrices, keeping 21 restarts alongside 400
    // orderings *lowers* rename reproducibility (0.974-0.978 vs 0.988-0.995)
    // while doubling the runtime, because each restart adds unaveraged random
    // variation. So when orderings are in use and --restart was not given
    // explicitly, default to none. An explicit --restart is always honoured, and
    // --hc-orders 1 (or --legacy-lrw) keeps the published default of 21.
    if params.hc_orders > 1 && !restart_set && !params.legacy_lrw {
        params.restart = 0;
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
        theme_network,
        json,
        scale_free,
    })
}


fn run() -> Result<(), String> {
    let args = parse_args()?;

    // Load and (optionally) discretize. The file text and the zero-copy string
    // slices are scoped so they are freed before the (heavier) reconstruction.
    let data = {
        let mut text = io::read_input(&args.input)?;
        if args.transpose {
            text = io::transpose_tsv(&text)?;
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
            Some(path) => Some(io::read_annotation_labels(
                &io::read_input(path)?,
                &result.nodes,
                args.atc_col.as_deref(),
            )?),
            None => None,
        };
        let acfg = analysis::AnalysisConfig {
            weighted: args.metrics_weighted,
            theme_network: args.theme_network,
        };
        for (ti, thr) in result.thresholds.iter().enumerate() {
            let net = &result.consensus[ti].network;
            let annot = atc_labels.as_ref().map(|(l, t)| (l.as_slice(), t.as_slice()));
            let an = analysis::analyze(net, result.nodes.len(), annot, &acfg);

            let base: PathBuf = [
                &args.outdir,
                "consensusGraph",
                &format!("globalNet_presFreq{thr}"),
            ]
            .iter()
            .collect();
            eprintln!(
                "# Analysis @ presFreq {thr}: {} nodes, {} edges, {} communities, modularity {:.3}{}",
                an.metrics.n_nodes,
                an.metrics.n_edges,
                an.metrics.n_communities,
                an.metrics.modularity,
                if an.atc.is_some() {
                    format!(", ATC over {} themes", an.themes.len())
                } else {
                    String::new()
                },
            );

            let atc_ref = match (&an.atc, &an.atc_last) {
                (Some(a), Some(al)) => Some((a.as_slice(), al.as_slice())),
                _ => None,
            };
            output::write_node_analysis(&base.join("nodeAnalysis.txt"), &result.nodes, &an.metrics, atc_ref)
                .map_err(|e| format!("write node analysis: {e}"))?;
            if !an.enrichment.is_empty() {
                output::write_enrichment(&base.join("enrichment.txt"), &an.enrichment)
                    .map_err(|e| format!("write enrichment: {e}"))?;
            }
            if let Some(tp) = &an.theme_network {
                output::write_theme_network(&base.join("themeNetwork.txt"), tp)
                    .map_err(|e| format!("write theme network: {e}"))?;
            }
            if args.json {
                output::write_json(&base.join("graph.json"), &result.nodes, net, &an)
                    .map_err(|e| format!("write json: {e}"))?;
            }
        }
    }

    // Optional scale-free topology analysis: sweep |ecorr| thresholds, fit a
    // power-law degree distribution at each (P(k) ~ k^-gamma), and propose the
    // threshold that best approximates scale-free topology (Barabási / WGCNA).
    if args.scale_free {
        for (ti, thr) in result.thresholds.iter().enumerate() {
            let net = &result.consensus[ti].network;
            let sweep = scalefree::sweep_ecorr(net, result.nodes.len(), 20, Sign::Both);
            eprintln!("# Scale-free topology @ presFreq {thr}   (P(k) ~ k^-gamma)");
            eprintln!("#   min|ecorr|  nodes  edges   gamma      R2   scale-free");
            for r in &sweep {
                match r.fit {
                    Some(f) => eprintln!(
                        "#       {:>6.3}  {:>5}  {:>5}   {:>5.2}   {:>5.3}   {}",
                        r.threshold, r.n_nodes, r.n_edges, f.gamma, f.r2,
                        if f.is_scale_free() { "yes" } else { "" }
                    ),
                    None => eprintln!(
                        "#       {:>6.3}  {:>5}  {:>5}       -       -",
                        r.threshold, r.n_nodes, r.n_edges
                    ),
                }
            }
            match scalefree::recommend(&sweep, 0.25) {
                Some(rec) => {
                    let f = rec.fit.unwrap();
                    eprintln!(
                        "#   -> recommended min|ecorr| ~ {:.3}: gamma={:.2}, R2={:.3}, keeps {:.0}% of nodes",
                        rec.threshold, f.gamma, f.r2, rec.frac_nodes * 100.0
                    );
                }
                None => eprintln!("#   -> no threshold kept enough of the network to fit"),
            }
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
