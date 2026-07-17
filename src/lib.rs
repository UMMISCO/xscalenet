//! xscalenet — efficient Rust reimplementation of the ScaleNet network
//! reconstruction pipeline and its Spectral Consensus Strategy (core path).
//!
//! Library entry point: build a [`data::Dataset`], then call [`run`] with a
//! [`RunConfig`] to get a typed [`PipelineResult`] (per-method and consensus
//! [`Network`]s). The CLI is a thin wrapper that loads the data and formats the
//! result; see [`output`] for the R-compatible writers.

pub mod analysis;
pub mod aracne;
pub mod atc;
pub mod data;
pub mod discretize;
pub mod enrichment;
pub mod hc;
pub mod io;
pub mod lanczos;
pub mod metrics;
pub mod mi;
pub mod multilayer;
pub mod network;
pub mod output;
/// Sequential Rayon shim, compiled only when the `parallel` feature is off.
mod par;
/// Rich `graph.json` builder shared by the server and the wasm build.
#[cfg(feature = "graphjson")]
pub mod graphjson;
/// In-browser (wasm-bindgen) reconstruction entry point, `--features wasm`.
#[cfg(feature = "wasm")]
pub mod wasm;
pub mod pipeline;
pub mod recon;
pub mod scalefree;
pub mod scalenet;
pub mod scs;
pub mod spectral;
pub mod util;

// Public API surface.
pub use network::{
    Consensus, Edge, Network, NetworkKind, Orientation, PipelineResult, RankedEdge, SpectralStats,
};
pub use pipeline::{run, RunConfig};
pub use scalenet::{Method, Params};
