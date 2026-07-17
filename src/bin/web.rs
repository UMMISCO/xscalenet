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


use xscalenet::data::Dataset;
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
    let title = format!("Imported — {fname}");
    xscalenet::graphjson::build(
        &title,
        &data,
        &config,
        flag("weighted"),
        num_u("threshold_index", 0),
        annotation.as_deref(),
        f.get("color_by").map(String::as_str),
    )
}
