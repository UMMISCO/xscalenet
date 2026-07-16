//! `xscalenet-web` — a small local server that reconstructs a network from a raw
//! matrix in the browser's "Import" panel, reusing the validated Rust core.
//!
//! Serves the static `web/` frontend and exposes `POST /run` (multipart: a matrix
//! file + reconstruction params [+ an annotation file]) → the `graph.json` the
//! network browser loads. Run with:
//!
//!   cargo run --features web --bin xscalenet-web            # serves ./web on :8788
//!
//! Everything runs locally; the heavy reconstruction happens on a blocking thread
//! pool so the async runtime stays responsive.

use std::collections::HashMap;

use axum::{
    extract::{DefaultBodyLimit, Multipart},
    http::StatusCode,
    response::{IntoResponse, Response},
    routing::{get, post},
    Router,
};
use tower_http::services::ServeDir;

use serde_json::{json, Value};

use xscalenet::atc::{self, AtcConfig};
use xscalenet::data::Dataset;
use xscalenet::metrics::{self, Graph};
use xscalenet::network::Orientation;
use xscalenet::pipeline::RunConfig;
use xscalenet::scalenet::{Method, Params};
use xscalenet::{discretize, io};

#[tokio::main]
async fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| "web".to_string());
    let port: u16 = std::env::var("PORT").ok().and_then(|p| p.parse().ok()).unwrap_or(8788);

    let app = Router::new()
        .route("/health", get(|| async { "ok" }))
        .route("/run", post(run_handler))
        // up to 512 MB uploads (SCAPIS matrices are large)
        .layer(DefaultBodyLimit::max(512 * 1024 * 1024))
        .fallback_service(ServeDir::new(&dir));

    let listener = tokio::net::TcpListener::bind(("127.0.0.1", port)).await.expect("bind");
    println!("xscalenet-web → http://localhost:{port}   (serving ./{dir})");
    println!("  open http://localhost:{port}/network.html and use “Import data”.");
    axum::serve(listener, app).await.expect("serve");
}

fn bad(msg: String) -> Response {
    (StatusCode::BAD_REQUEST, msg).into_response()
}

async fn run_handler(mut mp: Multipart) -> Response {
    let mut data: Option<(String, Vec<u8>)> = None;
    let mut annotation: Option<String> = None;
    let mut fields: HashMap<String, String> = HashMap::new();

    loop {
        let field = match mp.next_field().await {
            Ok(Some(f)) => f,
            Ok(None) => break,
            Err(e) => return bad(format!("multipart error: {e}")),
        };
        let name = field.name().unwrap_or("").to_string();
        let fname = field.file_name().unwrap_or("").to_string();
        match name.as_str() {
            "data" => match field.bytes().await {
                Ok(b) => data = Some((fname, b.to_vec())),
                Err(e) => return bad(format!("read data: {e}")),
            },
            "annotation" => match field.bytes().await {
                Ok(b) => annotation = Some(String::from_utf8_lossy(&b).into_owned()),
                Err(e) => return bad(format!("read annotation: {e}")),
            },
            _ => {
                if let Ok(t) = field.text().await {
                    fields.insert(name, t);
                }
            }
        }
    }

    let Some((fname, bytes)) = data else {
        return bad("missing 'data' matrix file".into());
    };

    // Reconstruction is CPU-bound and synchronous (rayon) — keep it off the async runtime.
    match tokio::task::spawn_blocking(move || build_graph_json(fname, bytes, annotation, fields)).await {
        Ok(Ok(json)) => ([("content-type", "application/json")], json).into_response(),
        Ok(Err(e)) => bad(e),
        Err(e) => (StatusCode::INTERNAL_SERVER_ERROR, format!("task failed: {e}")).into_response(),
    }
}

/// Decode uploaded bytes to text, transparently gunzipping `.gz` uploads.
fn decode(bytes: &[u8], name: &str) -> Result<String, String> {
    if name.ends_with(".gz") || bytes.starts_with(&[0x1f, 0x8b]) {
        use std::io::Read;
        let mut s = String::new();
        flate2::read::GzDecoder::new(bytes)
            .read_to_string(&mut s)
            .map_err(|e| format!("gunzip: {e}"))?;
        Ok(s)
    } else {
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }
}

fn build_graph_json(
    fname: String,
    bytes: Vec<u8>,
    annotation: Option<String>,
    f: HashMap<String, String>,
) -> Result<String, String> {
    let flag = |k: &str| f.get(k).map(|v| v == "true" || v == "1" || v == "on").unwrap_or(false);
    let num_f = |k: &str, d: f64| f.get(k).and_then(|v| v.parse().ok()).unwrap_or(d);
    let num_u = |k: &str, d: usize| f.get(k).and_then(|v| v.parse().ok()).unwrap_or(d);

    // Load matrix (+ optional transpose / discretize), exactly as the CLI does.
    let text = decode(&bytes, &fname)?;
    let text = if flag("transpose") { io::transpose_tsv(&text)? } else { text };
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text)?;
    let data = if flag("discretize") {
        discretize::discretize_to_dataset(names, &raw, nrow, num_u("max_clusters", 5))
    } else {
        Dataset::from_raw_columns(names, &raw, nrow)
    };
    if data.ncol() < 2 {
        return Err(format!("need ≥2 variables, got {}", data.ncol()));
    }

    let methods: Vec<Method> = f
        .get("methods")
        .map(|s| {
            s.split(',')
                .filter_map(|m| match m.trim() {
                    "aracne" => Some(Method::Aracne),
                    "bayes_hc" => Some(Method::BayesHc),
                    _ => None,
                })
                .collect()
        })
        .filter(|v: &Vec<Method>| !v.is_empty())
        .unwrap_or_else(|| vec![Method::Aracne, Method::BayesHc]);
    let thresh: Vec<f64> = f
        .get("thresh")
        .map(|s| s.split(',').filter_map(|t| t.trim().parse().ok()).collect())
        .filter(|v: &Vec<f64>| !v.is_empty())
        .unwrap_or_else(|| vec![0.3, 0.8]);

    let params = Params {
        var_perc: num_f("var_perc", 0.2),
        eigen_perc: f.get("eigen_perc").and_then(|v| if v == "auto" || v.is_empty() { None } else { v.parse().ok() }),
        pres_freq_thresh: thresh,
        epsilon: num_f("epsilon", 0.001),
        restart: num_u("restart", 21),
        seed: f.get("seed").and_then(|v| v.parse().ok()).unwrap_or(6196),
        max_iter: num_u("max_iter", 10_000),
        iss: num_f("iss", 1.0),
        maxp: f.get("maxp").and_then(|v| v.parse().ok()).unwrap_or(usize::MAX),
        iterative_eigen: flag("iterative_eigen"),
    };

    let config = RunConfig { methods, params };
    let result = xscalenet::pipeline::run(&data, &config, |_| {})?;
    if result.consensus.is_empty() {
        return Err("no consensus network produced".into());
    }
    let ti = num_u("threshold_index", 0).min(result.consensus.len() - 1);
    let net = &result.consensus[ti].network;
    let n = result.nodes.len();
    let weighted = flag("weighted");

    // Topology (degree / betweenness / community …).
    let g = Graph::from_network(net, n, weighted);
    let gm = metrics::compute(&g);

    // Full node annotation table (all columns), keyed by node id.
    let (an_header, an_map) = annotation.as_deref().map(parse_annotation).unwrap_or_default();
    let attr_cols: Vec<&str> = an_header.iter().skip(1).map(String::as_str).collect();
    let col_index = |c: &str| an_header.iter().position(|h| h == c);
    let node_val = |i: usize, ci: usize| -> Option<&str> {
        an_map.get(&result.nodes[i]).and_then(|row| row.get(ci)).map(String::as_str)
    };

    // Colour-by column: explicit param, else a taxonomy rank, else the first attribute.
    let color_by: Option<String> = f
        .get("color_by")
        .cloned()
        .filter(|c| attr_cols.contains(&c.as_str()))
        .or_else(|| {
            ["phylum_group", "phylum", "genus", "family", "class", "order"]
                .iter()
                .find(|c| attr_cols.contains(c))
                .map(|s| s.to_string())
        })
        .or_else(|| f.get("atc_col").cloned().filter(|c| attr_cols.contains(&c.as_str())))
        .or_else(|| attr_cols.first().map(|s| s.to_string()));

    // colour_groups / colour_map over ALL values of the colour-by column, then
    // themes = the non-missing groups (same order) with colours FROM color_map so
    // legend, node fill and propagation pies all agree; then ATC on those themes.
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
            let col = if is_missing(g0) { "#c9ced6".to_string() } else { PALETTE[i % PALETTE.len()].to_string() };
            color_map.insert(g0.clone(), json!(col));
        }
        themes = color_groups.iter().filter(|g0| !is_missing(g0)).cloned().collect();
        theme_colors = themes.iter().map(|t| color_map[t].as_str().unwrap().to_string()).collect();

        let idx: HashMap<&String, usize> = themes.iter().enumerate().map(|(i, t)| (t, i)).collect();
        let labels: Vec<Option<usize>> = (0..n)
            .map(|i| node_val(i, cb).and_then(|v| idx.get(&v.to_string()).copied()))
            .collect();
        let t0 = atc::membership_from_labels(&labels, themes.len());
        let acfg = AtcConfig { weighted, self_loop: true, coinertia_threshold: 0.999, max_iter: 20, capture_trace: true };
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
            if let Some(a) = &atc { nd.insert("atc".into(), json!(round6(a[i]))); }
            if let Some(a) = &atc_last { nd.insert("atc_last".into(), json!(round6(a[i]))); }
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

    let graph = json!({
        "meta": { "title": format!("Imported — {fname}"), "source": "xscalenet reconstruction",
            "directed": false, "color_by": color_by, "n_nodes": n, "n_edges": net.edges.len() },
        "color_groups": color_groups,
        "color_map": Value::Object(color_map),
        "themes": themes,
        "theme_colors": theme_colors,
        "nodes": nodes_json,
        "edges": edges_json,
        "atc_trace": atc_trace,
    });
    Ok(serde_json::to_string(&graph).unwrap())
}

const PALETTE: [&str; 12] = ["#2f6df6","#3fb950","#db6d28","#a371f7","#e3b341","#e5534b","#39c5cf","#db61a2","#7ee787","#ffa657","#79c0ff","#d2a8ff"];
fn is_missing(s: &str) -> bool {
    matches!(s, "NA" | "Unclassified" | "unclassified" | "" | "unknown" | "Unknown" | "NaN")
}
fn round6(x: f64) -> f64 { if x.is_finite() { (x * 1e6).round() / 1e6 } else { 0.0 } }
// Parse an annotation table (id in column 1) into (header, id -> row). Delimiter
// auto-detected (tab / semicolon / comma).
fn parse_annotation(text: &str) -> (Vec<String>, HashMap<String, Vec<String>>) {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(head_line) = lines.next() else { return (Vec::new(), HashMap::new()) };
    let delim = if head_line.contains('\t') { '\t' } else if head_line.contains(';') { ';' } else { ',' };
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
