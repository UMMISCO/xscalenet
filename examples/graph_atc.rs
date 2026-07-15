//! Enrich an existing `graph.json` with ATC (annotation propagation) computed on
//! its own edges, annotated by one categorical node field.
//!
//! Adds, in place: `atc` / `atc_last` on every node, and top-level `themes`,
//! `theme_colors`, and an `atc_trace` block (the per-step propagation, same shape
//! the propagation player reads). Uses the validated [`xscalenet::atc`] core, so
//! the scores match the CLI / R reference.
//!
//! Usage:  cargo run --example graph_atc -- <graph.json> <node-field> [out.json]
//! e.g.    cargo run --example graph_atc -- web/scapis.graph.json phylum_group

use std::fs;

use serde_json::{json, Value};

use xscalenet::atc::{self, AtcConfig};
use xscalenet::metrics::Graph;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: graph_atc <graph.json> <field> [out.json]");
    let field = args.get(2).expect("need a node field to annotate by (e.g. phylum_group)");
    let out = args.get(3).cloned().unwrap_or_else(|| path.clone());

    let mut g: Value = serde_json::from_str(&fs::read_to_string(path).expect("read")).expect("parse json");

    // ---- nodes + index ------------------------------------------------------
    let nodes = g["nodes"].as_array().expect("nodes[]").clone();
    let n = nodes.len();
    let mut index = std::collections::HashMap::new();
    for (i, nd) in nodes.iter().enumerate() {
        index.insert(nd["id"].as_str().unwrap().to_string(), i);
    }

    // ---- themes: reuse color_groups/color_map when they describe this field --
    let (themes, theme_colors): (Vec<String>, Vec<String>) = if g.get("color_groups").is_some()
        && g["meta"].get("color_by").and_then(|v| v.as_str()) == Some(field.as_str())
    {
        let groups: Vec<String> = g["color_groups"].as_array().unwrap().iter()
            .map(|v| v.as_str().unwrap().to_string()).collect();
        let colors = groups.iter()
            .map(|gp| g["color_map"][gp].as_str().unwrap_or("#888888").to_string()).collect();
        (groups, colors)
    } else {
        // distinct values, most frequent first
        use std::collections::HashMap;
        let mut freq: HashMap<String, usize> = HashMap::new();
        for nd in &nodes {
            let v = nd[field].as_str().unwrap_or("NA").to_string();
            *freq.entry(v).or_default() += 1;
        }
        let mut gs: Vec<String> = freq.keys().cloned().collect();
        gs.sort_by_key(|k| std::cmp::Reverse(freq[k]));
        let pal = ["#2f6df6","#3fb950","#db6d28","#a371f7","#e3b341","#e5534b","#39c5cf","#db61a2"];
        let cs = (0..gs.len()).map(|i| pal[i % pal.len()].to_string()).collect();
        (gs, cs)
    };
    let theme_idx: std::collections::HashMap<&String, usize> =
        themes.iter().enumerate().map(|(i, t)| (t, i)).collect();
    let labels: Vec<Option<usize>> = nodes.iter()
        .map(|nd| theme_idx.get(&nd[field].as_str().unwrap_or("NA").to_string()).copied())
        .collect();

    // ---- graph from edges ---------------------------------------------------
    let mut adj = vec![Vec::new(); n];
    for e in g["edges"].as_array().expect("edges[]") {
        let (a, b) = (e["source"].as_str().unwrap(), e["target"].as_str().unwrap());
        if let (Some(&ia), Some(&ib)) = (index.get(a), index.get(b)) {
            adj[ia].push((ib, 1.0));
            adj[ib].push((ia, 1.0));
        }
    }
    let graph = Graph { n, adj };

    // ---- propagate ----------------------------------------------------------
    let t0 = atc::membership_from_labels(&labels, themes.len());
    let cfg = AtcConfig { weighted: false, self_loop: true, coinertia_threshold: 0.999,
        max_iter: 20, capture_trace: true };
    let res = atc::compute(&graph, &t0, &cfg);
    let trace = res.trace.expect("trace");

    // ---- write back: per-node scores + themes + trace -----------------------
    let mut nodes = nodes;
    for (i, nd) in nodes.iter_mut().enumerate() {
        nd["atc"] = json!(round6(res.atc[i]));
        nd["atc_last"] = json!(round6(res.atc_last_config[i]));
    }
    g["nodes"] = Value::Array(nodes);
    g["themes"] = json!(themes);
    g["theme_colors"] = json!(theme_colors);

    let steps: Vec<Value> = (0..trace.configs.len()).map(|p| {
        let config: Vec<Vec<f64>> = trace.configs[p].iter()
            .map(|row| row.iter().map(|&x| round6(x)).collect()).collect();
        if p == 0 {
            json!({ "step": p, "coinertia": Value::Null, "config": config, "delta": Value::Null })
        } else {
            let delta: Vec<f64> = trace.delta[p - 1].iter().map(|&x| round6(x)).collect();
            json!({ "step": p, "coinertia": round6(trace.coinertia[p - 1]), "config": config, "delta": delta })
        }
    }).collect();
    g["atc_trace"] = json!({ "annotation": field, "n_iterations": res.n_iterations, "steps": steps });
    if let Some(m) = g["meta"].as_object_mut() {
        m.insert("atc_annotation".into(), json!(field));
    }

    fs::write(&out, serde_json::to_string(&g).unwrap()).expect("write");
    println!("enriched {out}: ATC by '{field}' over {} themes, {} propagation steps, {} nodes",
        themes.len(), trace.configs.len() - 1, n);
}

fn round6(x: f64) -> f64 {
    if x.is_finite() { (x * 1e6).round() / 1e6 } else { 0.0 }
}
