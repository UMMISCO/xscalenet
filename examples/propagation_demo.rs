//! Generate the built-in propagation demo for the web app.
//!
//! Emits `web/propagation.graph.json` from the exact 10-node graph `J` and weighted
//! membership `m0` that xscalenet validates against the R `dynamicPropagation`
//! reference (see `src/atc.rs::matches_r_dynamic_propagation`). The JSON carries
//! the network, per-node metrics/ATC, and the full step-by-step ATC propagation
//! trace (`Tᵖ` configurations) so the frontend can animate each step.
//!
//! Run with: `cargo run --example propagation_demo`

use std::fmt::Write as _;
use std::fs;

use xscalenet::atc::{self, AtcConfig};
use xscalenet::metrics::{self, Graph};

fn num(x: f64) -> String {
    if x.is_finite() {
        // Trim to 6 significant decimals, drop trailing zeros.
        let s = format!("{x:.6}");
        let t = s.trim_end_matches('0').trim_end_matches('.');
        if t.is_empty() || t == "-0" { "0".into() } else { t.to_string() }
    } else {
        "null".into()
    }
}

fn arr(v: &[f64]) -> String {
    let inner: Vec<String> = v.iter().map(|&x| num(x)).collect();
    format!("[{}]", inner.join(", "))
}

fn main() {
    // Graph J (0-indexed). Node 3 (g4) is isolated — the floating pie in the paper figure.
    let edges = [
        (0, 1),
        (0, 2),
        (1, 2),
        (1, 4),
        (2, 4),
        (2, 5),
        (5, 6),
        (5, 7),
        (6, 7),
        (7, 8),
        (8, 9),
    ];
    let n = 10;
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in &edges {
        adj[a].push((b, 1.0));
        adj[b].push((a, 1.0));
    }
    let g = Graph { n, adj };

    // Weighted, multi-theme membership m0 over themes {brown, black, yellow}.
    let m0: Vec<Vec<f64>> = vec![
        vec![1., 0., 0.], // g1
        vec![0., 0., 0.], // g2
        vec![0., 1., 0.], // g3
        vec![0., 0., 4.], // g4 (isolated)
        vec![1., 1., 0.], // g5
        vec![0., 0., 0.], // g6
        vec![1., 0., 0.], // g7
        vec![0., 0., 1.], // g8
        vec![0., 0., 1.], // g9
        vec![0., 0., 0.], // g10
    ];
    let themes = ["brown", "black", "yellow"];
    let theme_colors = ["#9c6b3f", "#2b2b2b", "#e3b341"];

    // Run the propagation, capturing every configuration. Convergence (RV ≥ 0.9)
    // is reached by step 2, but we run further so the animation shows the full
    // approach to steady state.
    let cfg = AtcConfig {
        weighted: false,
        self_loop: true,
        coinertia_threshold: 0.9995,
        max_iter: 15,
        capture_trace: true,
    };
    let res = atc::compute(&g, &m0, &cfg);
    let trace = res.trace.expect("trace requested");

    // Topology (community / degree / betweenness) for layout + node info.
    let gm = metrics::compute(&g);
    let converged_at = res.coinertia.iter().position(|&rv| rv >= 0.9).map(|p| p + 1);

    let mut out = String::new();
    out.push_str("{\n");
    write!(
        out,
        "  \"meta\": {{\"title\": \"ATC propagation — 10-node test network\", \"source\": \"Prifti et al. 2010 dynamicPropagation (validated)\", \"n_steps\": {}, \"converged_at\": {}}},\n",
        trace.configs.len() - 1,
        converged_at.map(|s| s.to_string()).unwrap_or_else(|| "null".into())
    )
    .unwrap();

    // Themes.
    let tn: Vec<String> = themes.iter().map(|t| format!("\"{t}\"")).collect();
    let tc: Vec<String> = theme_colors.iter().map(|c| format!("\"{c}\"")).collect();
    write!(out, "  \"themes\": [{}],\n", tn.join(", ")).unwrap();
    write!(out, "  \"theme_colors\": [{}],\n", tc.join(", ")).unwrap();

    // Nodes.
    out.push_str("  \"nodes\": [\n");
    for i in 0..n {
        let nm = &gm.nodes[i];
        let comma = if i + 1 < n { "," } else { "" };
        write!(
            out,
            "    {{\"id\": \"g{}\", \"community\": {}, \"degree\": {}, \"betweenness\": {}, \"atc\": {}, \"atc_last\": {}, \"m0\": {}}}{}\n",
            i + 1,
            nm.community,
            nm.degree,
            num(nm.betweenness),
            num(res.atc[i]),
            num(res.atc_last_config[i]),
            arr(&m0[i]),
            comma
        )
        .unwrap();
    }
    out.push_str("  ],\n");

    // Edges.
    out.push_str("  \"edges\": [\n");
    for (k, &(a, b)) in edges.iter().enumerate() {
        let comma = if k + 1 < edges.len() { "," } else { "" };
        write!(out, "    {{\"source\": \"g{}\", \"target\": \"g{}\"}}{}\n", a + 1, b + 1, comma).unwrap();
    }
    out.push_str("  ],\n");

    // ATC trace: T⁰ .. Tᵖ, each with the node×theme config, plus per-step
    // co-inertia and per-node change.
    out.push_str("  \"atc_trace\": {\n    \"steps\": [\n");
    let n_steps = trace.configs.len();
    for p in 0..n_steps {
        let comma = if p + 1 < n_steps { "," } else { "" };
        // Per-node config rows.
        let rows: Vec<String> = (0..n).map(|i| arr(&trace.configs[p][i])).collect();
        let config = format!("[{}]", rows.join(", "));
        // coinertia/delta align with configs[1..]; step 0 has none.
        let (coin, delta) = if p == 0 {
            ("null".to_string(), "null".to_string())
        } else {
            (num(trace.coinertia[p - 1]), arr(&trace.delta[p - 1]))
        };
        write!(
            out,
            "      {{\"step\": {p}, \"coinertia\": {coin}, \"config\": {config}, \"delta\": {delta}}}{comma}\n"
        )
        .unwrap();
    }
    out.push_str("    ]\n  }\n}\n");

    let path = "web/propagation.graph.json";
    fs::write(path, &out).expect("write example.graph.json");
    println!(
        "wrote {path}: {} nodes, {} steps (converged at step {:?}), {} bytes",
        n,
        n_steps - 1,
        converged_at,
        out.len()
    );
}
