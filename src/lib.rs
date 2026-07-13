//! xscalenet — efficient Rust reimplementation of the ScaleNet network
//! reconstruction pipeline and its Spectral Consensus Strategy (core path).
//!
//! Library entry point: build a [`data::Dataset`], then call [`run`] with a
//! [`RunConfig`] to get a typed [`PipelineResult`] (per-method and consensus
//! [`Network`]s). The CLI is a thin wrapper that loads the data and formats the
//! result; see [`output`] for the R-compatible writers.

pub mod aracne;
pub mod atc;
pub mod data;
pub mod discretize;
pub mod hc;
pub mod metrics;
pub mod mi;
pub mod network;
pub mod output;
pub mod pipeline;
pub mod recon;
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
