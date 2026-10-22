//! Rich `graph.json` builder — shared by the `xscalenet-web` server and the wasm
//! (`--features wasm`) in-browser reconstruction, so both emit an identical
//! interchange file: per-node topology metrics + annotation columns + ATC, the
//! consensus edges, the propagation `atc_trace`, and the spectral subgraphs.
#![cfg(feature = "graphjson")]

use std::collections::HashMap;

use serde_json::{json, Value};

use crate::atc::{self, AtcConfig};
use crate::data::Dataset;
use crate::metrics::{self, Graph};
use crate::network::Orientation;
use crate::pipeline::{run, RunConfig};

pub const PALETTE: [&str; 12] = [
    "#2f6df6", "#3fb950", "#db6d28", "#a371f7", "#e3b341", "#e5534b", "#39c5cf", "#db61a2",
    "#7ee787", "#ffa657", "#79c0ff", "#d2a8ff",
];

/// Values that mean "this node carries no annotation".
///
/// Kept in step with `ATC_MISSING` in `web/network.html`: the browser recomputes
/// ATC on the visible subgraph, and if the two disagree about what counts as
/// missing they build different theme axes and the recomputed scores are no
/// longer the same quantity as the engine's. `"unannotated"` was in neither list
/// despite being exactly that, so the bundled 10-node reference propagated over
/// 3 themes in the engine and 4 in the browser.
pub fn is_missing(s: &str) -> bool {
    matches!(
        s,
        "NA" | "Unclassified"
            | "unclassified"
            | ""
            | "unknown"
            | "Unknown"
            | "NaN"
            | "unannotated"
            | "Unannotated"
    )
}
fn round6(x: f64) -> f64 {
    if x.is_finite() {
        (x * 1e6).round() / 1e6
    } else {
        0.0
    }
}

/// Parse an annotation table (id in column 1) into (header, id -> row). Delimiter
/// auto-detected (tab / semicolon / comma).
pub fn parse_annotation(text: &str) -> (Vec<String>, HashMap<String, Vec<String>>) {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(head_line) = lines.next() else {
        return (Vec::new(), HashMap::new());
    };
    let delim = if head_line.contains('\t') {
        '\t'
    } else if head_line.contains(';') {
        ';'
    } else {
        ','
    };
    let clean = |s: &str| s.trim().trim_matches('"').to_string();
    let header: Vec<String> = head_line.split(delim).map(clean).collect();
    let mut map = HashMap::new();
    for line in lines {
        let row: Vec<String> = line.split(delim).map(clean).collect();
        if let Some(id) = row.first() {
            map.insert(id.clone(), row);
        }
    }
    (header, map)
}

/// Run the reconstruction pipeline and serialize the rich `graph.json`.
///
/// * `annotation` — optional annotation-table text (id + attribute columns).
/// * `color_by` — optional explicit colour-by column; otherwise a taxonomy rank
///   or the first attribute is chosen. Drives the themes and the ATC.
pub fn build(
    title: &str,
    data: &Dataset,
    config: &RunConfig,
    weighted: bool,
    threshold_index: usize,
    annotation: Option<&str>,
    color_by: Option<&str>,
) -> Result<String, String> {
    let result = run(data, config, |_| {})?;
    if result.consensus.is_empty() {
        return Err("no consensus network produced".into());
    }
    let ti = threshold_index.min(result.consensus.len() - 1);
    let net = &result.consensus[ti].network;
    let n = result.nodes.len();

    // Topology (degree / betweenness / community …).
    let g = Graph::from_network(net, n, weighted);
    let gm = metrics::compute(&g);

    // Full node annotation table (all columns), keyed by node id.
    let (an_header, an_map) = annotation.map(parse_annotation).unwrap_or_default();
    let attr_cols: Vec<&str> = an_header.iter().skip(1).map(String::as_str).collect();
    let col_index = |c: &str| an_header.iter().position(|h| h == c);
    let node_val = |i: usize, ci: usize| -> Option<&str> {
        an_map
            .get(&result.nodes[i])
            .and_then(|row| row.get(ci))
            .map(String::as_str)
    };

    // Colour-by column: explicit arg, else a taxonomy rank, else the first attribute.
    let color_by: Option<String> = color_by
        .map(str::to_string)
        .filter(|c| attr_cols.contains(&c.as_str()))
        .or_else(|| {
            ["phylum_group", "phylum", "genus", "family", "class", "order"]
                .iter()
                .find(|c| attr_cols.contains(c))
                .map(|s| s.to_string())
        })
        .or_else(|| attr_cols.first().map(|s| s.to_string()));

    // colour_groups / colour_map over ALL values of the colour-by column, then
    // themes = the non-missing groups (colours FROM color_map so legend, fill and
    // propagation pies agree); then ATC on those themes.
    let mut color_groups: Vec<String> = Vec::new();
    let mut color_map = serde_json::Map::new();
    let mut themes: Vec<String> = Vec::new();
    let mut theme_colors: Vec<String> = Vec::new();
    let mut atc = None;
    let mut atc_last = None;
    let mut trace = None;
    if let Some(cb) = color_by.as_ref().and_then(|c| col_index(c)) {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for i in 0..n {
            *freq.entry(node_val(i, cb).unwrap_or("NA").to_string()).or_default() += 1;
        }
        color_groups = freq.keys().cloned().collect();
        color_groups.sort_by_key(|k| std::cmp::Reverse(freq[k]));
        for (i, g0) in color_groups.iter().enumerate() {
            let col = if is_missing(g0) {
                "#c9ced6".to_string()
            } else {
                PALETTE[i % PALETTE.len()].to_string()
            };
            color_map.insert(g0.clone(), json!(col));
        }
        themes = color_groups.iter().filter(|g0| !is_missing(g0)).cloned().collect();
        theme_colors = themes.iter().map(|t| color_map[t].as_str().unwrap().to_string()).collect();

        let idx: HashMap<&String, usize> = themes.iter().enumerate().map(|(i, t)| (t, i)).collect();
        let labels: Vec<Option<usize>> = (0..n)
            .map(|i| node_val(i, cb).and_then(|v| idx.get(&v.to_string()).copied()))
            .collect();
        let t0 = atc::membership_from_labels(&labels, themes.len());
        let acfg = AtcConfig {
            weighted,
            self_loop: true,
            coinertia_threshold: 0.999,
            max_iter: 20,
            capture_trace: true,
        };
        let r = atc::compute(&g, &t0, &acfg);
        atc = Some(r.atc);
        atc_last = Some(r.atc_last_config);
        trace = r.trace;
    }

    // nodes: id + topology + every annotation column + ATC.
    let nodes_json: Vec<Value> = (0..n)
        .map(|i| {
            let m = &gm.nodes[i];
            let mut nd = serde_json::Map::new();
            nd.insert("id".into(), json!(result.nodes[i]));
            nd.insert("degree".into(), json!(m.degree));
            nd.insert("betweenness".into(), json!(round6(m.betweenness)));
            nd.insert("closeness".into(), json!(round6(m.closeness)));
            nd.insert("clustering".into(), json!(round6(m.clustering)));
            nd.insert("community".into(), json!(m.community));
            for (ci, col) in an_header.iter().enumerate().skip(1) {
                if let Some(v) = node_val(i, ci) {
                    nd.insert(col.clone(), json!(v));
                }
            }
            if let Some(a) = &atc {
                nd.insert("atc".into(), json!(round6(a[i])));
            }
            if let Some(a) = &atc_last {
                nd.insert("atc_last".into(), json!(round6(a[i])));
            }
            Value::Object(nd)
        })
        .collect();

    let orient = |o: Orientation| match o {
        Orientation::Forward => 1,
        Orientation::Backward => -1,
        Orientation::Undirected => 2,
    };
    let edges_json: Vec<Value> = net
        .edges
        .iter()
        .map(|e| {
            json!({ "source": result.nodes[e.source], "target": result.nodes[e.target],
                "ecorr": e.ecorr, "score": e.presence_score, "orient": orient(e.orientation),
                "orient_score": e.orient_score })
        })
        .collect();

    let atc_trace = trace.map(|tr| {
        let steps: Vec<Value> = (0..tr.configs.len())
            .map(|p| {
                let config: Vec<Vec<f64>> = tr.configs[p]
                    .iter()
                    .map(|row| row.iter().map(|&x| round6(x)).collect())
                    .collect();
                if p == 0 {
                    json!({ "step": p, "coinertia": Value::Null, "config": config, "delta": Value::Null })
                } else {
                    let delta: Vec<f64> = tr.delta[p - 1].iter().map(|&x| round6(x)).collect();
                    json!({ "step": p, "coinertia": round6(tr.coinertia[p - 1]), "config": config, "delta": delta })
                }
            })
            .collect();
        json!({ "annotation": color_by, "n_iterations": tr.coinertia.len(), "steps": steps })
    });

    // Spectral subgraphs fused into this network — "the subnetworks that made it".
    let subgraphs_json: Vec<Value> = result
        .subgraphs
        .iter()
        .enumerate()
        .map(|(j, sg)| {
            json!({ "eigenvector": j / 2 + 1, "sign": if j % 2 == 0 { "pos" } else { "neg" }, "nodes": sg })
        })
        .collect();

    let graph = json!({
        "meta": { "title": title, "source": "xscalenet reconstruction",
            "directed": false, "color_by": color_by, "n_nodes": n, "n_edges": net.edges.len(),
            "k": result.spectral.k, "subset_m": result.spectral.subset_m },
        "color_groups": color_groups,
        "color_map": Value::Object(color_map),
        "themes": themes,
        "theme_colors": theme_colors,
        "nodes": nodes_json,
        "edges": edges_json,
        "atc_trace": atc_trace,
        "subgraphs": subgraphs_json,
    });
    Ok(serde_json::to_string(&graph).unwrap())
}
