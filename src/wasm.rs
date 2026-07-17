//! In-browser reconstruction entry point (wasm-bindgen).
//!
//! Compiled only with `--features wasm` (which also disables `parallel`, so the
//! sequential shim in [`crate::par`] is used). Exposes [`reconstruct`], which
//! runs the validated pipeline on a raw matrix and returns a `graph.json` string
//! the web frontend can load directly — no server, all client-side.
#![cfg(feature = "wasm")]

use serde_json::json;
use wasm_bindgen::prelude::*;

use crate::data::Dataset;
use crate::pipeline::{run, RunConfig};
use crate::scalenet::{Method, Params};
use crate::{discretize, io};

/// Reconstruct a network from a raw matrix, in the browser.
///
/// * `matrix_text` — TSV/CSV, variables in columns & samples in rows (or tick
///   `transpose`).
/// * `methods_csv` — e.g. `"aracne"` or `"aracne,bayes_hc"` (default: aracne).
/// * `thresh_csv` — presence-frequency thresholds, e.g. `"0.3"` (default 0.3).
///
/// Returns a `graph.json` string, or a JS error.
#[wasm_bindgen]
pub fn reconstruct(
    matrix_text: &str,
    methods_csv: &str,
    thresh_csv: &str,
    transpose: bool,
    discretize_flag: bool,
) -> Result<String, JsValue> {
    reconstruct_inner(matrix_text, methods_csv, thresh_csv, transpose, discretize_flag)
        .map_err(|e| JsValue::from_str(&e))
}

fn reconstruct_inner(
    matrix_text: &str,
    methods_csv: &str,
    thresh_csv: &str,
    transpose: bool,
    discretize_flag: bool,
) -> Result<String, String> {
    // Load matrix (+ optional transpose / discretize) — the same path as the CLI.
    let text = if transpose {
        io::transpose_tsv(matrix_text)?
    } else {
        matrix_text.to_string()
    };
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text)?;
    let data = if discretize_flag {
        discretize::discretize_to_dataset(names, &raw, nrow, 5)
    } else {
        Dataset::from_raw_columns(names, &raw, nrow)
    };

    let methods: Vec<Method> = {
        let m: Vec<Method> = methods_csv
            .split(',')
            .filter_map(|s| match s.trim().to_ascii_lowercase().as_str() {
                "aracne" => Some(Method::Aracne),
                "bayes_hc" | "hc" | "bayeshc" => Some(Method::BayesHc),
                _ => None,
            })
            .collect();
        if m.is_empty() {
            vec![Method::Aracne]
        } else {
            m
        }
    };
    let thresh: Vec<f64> = {
        let t: Vec<f64> = thresh_csv
            .split(',')
            .filter_map(|s| s.trim().parse().ok())
            .collect();
        if t.is_empty() {
            vec![0.3]
        } else {
            t
        }
    };

    let params = Params {
        var_perc: 0.2,
        eigen_perc: None,
        pres_freq_thresh: thresh,
        epsilon: 0.001,
        restart: 21,
        seed: 6196,
        max_iter: 10_000,
        iss: 1.0,
        maxp: usize::MAX,
        iterative_eigen: false,
    };
    let config = RunConfig { methods, params };
    let result = run(&data, &config, |_| {})?;

    // Minimal graph.json: nodes + the consensus edges (last threshold). Enough to
    // visualize (colour by ecorr sign, filter by |ecorr|, layout, …); richer
    // annotation/metrics are added client-side or via the server path.
    let names = &result.nodes;
    let cons = result
        .consensus
        .last()
        .ok_or_else(|| "no consensus network produced".to_string())?;
    let edges: Vec<_> = cons
        .network
        .edges
        .iter()
        .map(|e| {
            json!({
                "source": names[e.source],
                "target": names[e.target],
                "ecorr": e.ecorr,
                "score": e.presence_score,
            })
        })
        .collect();
    let nodes: Vec<_> = names.iter().map(|n| json!({ "id": n })).collect();
    let g = json!({
        "meta": {
            "title": format!("Imported — {} features (wasm)", names.len()),
            "source": format!("in-browser reconstruction · {} edges · k={}", edges.len(), result.spectral.k),
            "directed": false,
            "k": result.spectral.k,
            "subset_m": result.spectral.subset_m,
        },
        "nodes": nodes,
        "edges": edges,
    });
    Ok(g.to_string())
}
