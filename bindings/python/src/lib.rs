//! Python bindings for xscalenet: reconstruct a network and (optionally) run the
//! functional-layer analysis, returning a JSON string the caller parses with
//! `json.loads`. Mirrors the CLI on top of the same library entry points.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use ::xscalenet::analysis::{analyze, AnalysisConfig};
use ::xscalenet::data::Dataset;
use ::xscalenet::discretize;
use ::xscalenet::io;
use ::xscalenet::output::graph_json;
use ::xscalenet::pipeline::{run, RunConfig};
use ::xscalenet::scalenet::{Method, Params};

fn err(e: String) -> PyErr {
    PyValueError::new_err(e)
}

fn parse_methods(s: &str) -> Result<Vec<Method>, String> {
    s.split(',')
        .map(|m| match m.trim() {
            "aracne" => Ok(Method::Aracne),
            "bayes_hc" => Ok(Method::BayesHc),
            other => Err(format!("unknown method: {other}")),
        })
        .collect()
}

/// Reconstruct a network from `input` (TSV or gzip). With `metrics`/`annotation`
/// the functional layer (metrics, ATC, theme network) is added per threshold.
/// Returns a JSON string: `{ "nodes": [...], "thresholds": [...], "results": [
/// per-threshold graph ] }`.
#[pyfunction]
#[pyo3(signature = (
    input,
    methods = "aracne,bayes_hc",
    discretize = false,
    transpose = false,
    var_perc = 0.2,
    eigen_perc = None,
    thresh = "0.3,0.8",
    epsilon = 0.001,
    restart = 21,
    seed = 6196,
    max_iter = 10000,
    iss = 1.0,
    maxp = None,
    max_clusters = 5,
    metrics = false,
    annotation = None,
    atc_col = None,
    theme_network = false,
    weighted = false,
))]
#[allow(clippy::too_many_arguments)]
fn reconstruct(
    input: &str,
    methods: &str,
    discretize: bool,
    transpose: bool,
    var_perc: f64,
    eigen_perc: Option<f64>,
    thresh: &str,
    epsilon: f64,
    restart: usize,
    seed: u64,
    max_iter: usize,
    iss: f64,
    maxp: Option<usize>,
    max_clusters: usize,
    metrics: bool,
    annotation: Option<&str>,
    atc_col: Option<&str>,
    theme_network: bool,
    weighted: bool,
) -> PyResult<String> {
    // Build the dataset.
    let mut text = io::read_input(input).map_err(err)?;
    if transpose {
        text = io::transpose_tsv(&text).map_err(err)?;
    }
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text).map_err(err)?;
    let data = if discretize {
        discretize::discretize_to_dataset(names, &raw, nrow, max_clusters)
    } else {
        Dataset::from_raw_columns(names, &raw, nrow)
    };

    let thresholds: Vec<f64> = thresh
        .split(',')
        .map(|t| t.trim().parse::<f64>().map_err(|e| err(e.to_string())))
        .collect::<PyResult<Vec<_>>>()?;

    let params = Params {
        var_perc,
        eigen_perc,
        pres_freq_thresh: thresholds.clone(),
        epsilon,
        restart,
        seed,
        max_iter,
        iss,
        maxp: maxp.unwrap_or(usize::MAX),
    };
    let cfg = RunConfig {
        methods: parse_methods(methods).map_err(err)?,
        params,
    };

    let result = run(&data, &cfg, |_| {}).map_err(err)?;

    // Annotation labels (for ATC / enrichment / theme network).
    let annot = match annotation {
        Some(path) => {
            let atext = io::read_input(path).map_err(err)?;
            Some(io::read_annotation_labels(&atext, &result.nodes, atc_col).map_err(err)?)
        }
        None => None,
    };
    let want_analysis = metrics || annot.is_some();
    let acfg = AnalysisConfig {
        weighted,
        theme_network,
    };

    // Assemble the top-level JSON.
    let mut out = String::from("{\n\"nodes\": [");
    out.push_str(
        &result
            .nodes
            .iter()
            .map(|n| format!("{:?}", n))
            .collect::<Vec<_>>()
            .join(", "),
    );
    out.push_str("],\n\"thresholds\": [");
    out.push_str(
        &thresholds
            .iter()
            .map(|t| t.to_string())
            .collect::<Vec<_>>()
            .join(", "),
    );
    out.push_str("],\n\"results\": [\n");

    for (ti, _thr) in thresholds.iter().enumerate() {
        let net = &result.consensus[ti].network;
        if ti > 0 {
            out.push_str(",\n");
        }
        if want_analysis {
            let a = annot.as_ref().map(|(l, t)| (l.as_slice(), t.as_slice()));
            let an = analyze(net, result.nodes.len(), a, &acfg);
            out.push_str(&graph_json(&result.nodes, net, &an));
        } else {
            // Minimal: consensus edges only.
            out.push_str("{\"edges\": [");
            let edges = net
                .edges
                .iter()
                .map(|e| {
                    format!(
                        "{{\"source\": {:?}, \"target\": {:?}}}",
                        result.nodes[e.source], result.nodes[e.target]
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            out.push_str(&edges);
            out.push_str("]}");
        }
    }
    out.push_str("\n]}\n");
    Ok(out)
}

#[pymodule]
fn xscalenet(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reconstruct, m)?)?;
    Ok(())
}
