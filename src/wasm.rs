//! In-browser reconstruction entry point (wasm-bindgen).
//!
//! Compiled only with `--features wasm` (which also disables `parallel`, so the
//! sequential shim in [`crate::par`] is used). Exposes [`reconstruct`], which
//! runs the validated pipeline on a raw matrix and returns a `graph.json` string
//! the web frontend can load directly — no server, all client-side.
#![cfg(feature = "wasm")]

use wasm_bindgen::prelude::*;

use crate::data::Dataset;
use crate::pipeline::RunConfig;
use crate::scalenet::{Method, Params};
use crate::{discretize, io};

/// Reconstruct a network from a raw matrix, in the browser.
///
/// * `matrix_text` — TSV/CSV, variables in columns & samples in rows (or tick
///   `transpose`).
/// * `methods_csv` — e.g. `"aracne"` or `"aracne,bayes_hc"` (default: aracne).
/// * `thresh_csv` — presence-frequency thresholds, e.g. `"0.3"` (default 0.3).
/// * `annotation` — optional annotation-table text (id + attribute columns); empty = none.
/// * `color_by` — optional colour-by column name; empty = auto.
///
/// Returns the same **rich** `graph.json` the server produces (per-node metrics +
/// annotation + ATC, consensus edges, `atc_trace`, spectral subgraphs), or a JS error.
#[wasm_bindgen]
pub fn reconstruct(
    matrix_text: &str,
    methods_csv: &str,
    thresh_csv: &str,
    transpose: bool,
    discretize_flag: bool,
    annotation: &str,
    color_by: &str,
) -> Result<String, JsValue> {
    reconstruct_inner(
        matrix_text,
        methods_csv,
        thresh_csv,
        transpose,
        discretize_flag,
        annotation,
        color_by,
    )
    .map_err(|e| JsValue::from_str(&e))
}

#[allow(clippy::too_many_arguments)]
fn reconstruct_inner(
    matrix_text: &str,
    methods_csv: &str,
    thresh_csv: &str,
    transpose: bool,
    discretize_flag: bool,
    annotation: &str,
    color_by: &str,
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
        legacy_lrw: false,
        // Lower than the CLI default: the wasm build is single-threaded, so the
        // in-browser budget caps the orderings. Reproducibility is accordingly
        // weaker (see README); the server build uses the full default.
        hc_orders: 24,
        hc_order_freq: 0.0,
    };
    let config = RunConfig { methods, params };

    // Same rich graph.json as the server (metrics + annotation + ATC + subgraphs).
    let title = format!("Imported — {} features (wasm)", data.ncol());
    let annotation = if annotation.trim().is_empty() { None } else { Some(annotation) };
    let color_by = if color_by.trim().is_empty() { None } else { Some(color_by) };
    crate::graphjson::build(&title, &data, &config, false, 0, annotation, color_by)
}
