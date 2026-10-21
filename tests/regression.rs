//! Regression tests against bundled `pop2mat` fixtures. The golden ARACNE edge
//! lists were verified to match the original R `scalenet` package (present-edge
//! set and `epresenceScore` to ~13 significant digits; `ecorr` sign checked
//! against the original-value Spearman). These lock in that behaviour.

use std::collections::HashMap;
use std::path::PathBuf;

use xscalenet::data::Dataset;
use xscalenet::discretize;
use xscalenet::scalenet::{self, EdgeList, Method, Params};
use xscalenet::util::{fmt, pair_at};

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn params() -> Params {
    Params {
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
        // The R fixtures were produced by R's own order-dependent spectral step
        // (see src/spectral.rs); reproducing them requires the legacy path.
        legacy_lrw: true,
    }
}

fn load_pop2mat() -> Dataset {
    let text = std::fs::read_to_string(data_dir().join("pop2mat.txt")).unwrap();
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text).unwrap();
    discretize::discretize_to_dataset(names, &raw, nrow, 5)
}

/// Normalized "x\ty\tepresenceScore\tecorr" lines for present edges, sorted.
fn present_lines(data: &Dataset, el: &EdgeList) -> Vec<String> {
    let p = data.ncol();
    let mut out: Vec<String> = el
        .edges
        .iter()
        .map(|(&idx, rec)| {
            let (a, b) = pair_at(idx, p);
            let (xa, yb) = (&data.names[a], &data.names[b]);
            let (x, y) = if xa < yb { (xa, yb) } else { (yb, xa) };
            format!("{x}\t{y}\t{}\t{}", fmt(rec.epresence_score), fmt(rec.ecorr))
        })
        .collect();
    out.sort();
    out
}

#[test]
fn spectral_structure_matches_r() {
    let data = load_pop2mat();
    let prep = scalenet::prepare(&data, &params()).unwrap();
    assert_eq!(prep.spectral.k, 13, "k eigenvectors");
    assert_eq!(prep.subset_m, 9, "m vars per subgraph");
    assert_eq!(prep.subgraphs.len(), 24, "2*(k-1) subgraphs");
}

#[test]
fn aracne_matches_golden() {
    let data = load_pop2mat();
    let prep = scalenet::prepare(&data, &params()).unwrap();
    let ranks = data.rank_matrix();
    let mut lists = scalenet::run_method(&data, &prep, Method::Aracne, &params());
    for el in &mut lists {
        scalenet::fill_ecorr(&data, &ranks, el);
    }
    for (ti, thr) in ["0.3", "0.8"].iter().enumerate() {
        let golden: Vec<String> = std::fs::read_to_string(
            data_dir().join(format!("pop2mat.aracne.{thr}.golden")),
        )
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
        let got = present_lines(&data, &lists[ti]);
        assert_eq!(got, golden, "ARACNE present edges @ presFreq {thr}");
    }
}

#[test]
fn run_api_returns_typed_networks() {
    use xscalenet::network::{NetworkKind, Orientation};
    use xscalenet::pipeline::{self, RunConfig};

    let data = load_pop2mat();
    let config = RunConfig {
        methods: vec![Method::Aracne, Method::BayesHc],
        params: params(),
    };
    let result = pipeline::run(&data, &config, |_| {}).unwrap();

    assert_eq!(result.spectral.k, 13);
    assert_eq!(result.per_method.len(), 2); // two methods
    assert_eq!(result.consensus.len(), 2); // two thresholds
    assert_eq!(result.nodes.len(), 45);

    // ARACNE @0.3 present-edge count matches the golden, via the typed network.
    assert_eq!(result.per_method[0][0].edges.len(), 205);

    // ARACNE edges are undirected; consensus is a Consensus-kind network.
    assert!(result.per_method[0][0]
        .edges
        .iter()
        .all(|e| e.orientation == Orientation::Undirected));
    assert_eq!(result.consensus[0].network.kind, NetworkKind::Consensus);

    // Typed edges reference valid nodes with source < target.
    for e in &result.consensus[0].network.edges {
        assert!(e.source < e.target && e.target < result.nodes.len());
    }
}

#[test]
fn bayes_hc_is_deterministic() {
    let data = load_pop2mat();
    let prep = scalenet::prepare(&data, &params()).unwrap();
    let edge_set = |el: &EdgeList| -> HashMap<usize, f64> {
        el.edges.iter().map(|(&i, r)| (i, r.eorient)).collect()
    };
    let a = scalenet::run_method(&data, &prep, Method::BayesHc, &params());
    let b = scalenet::run_method(&data, &prep, Method::BayesHc, &params());
    for ti in 0..2 {
        assert_eq!(
            edge_set(&a[ti]),
            edge_set(&b[ti]),
            "bayes_hc reproducible across runs (threshold {ti})"
        );
    }
}
