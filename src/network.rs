//! Typed domain model for reconstructed networks — the public result of the
//! pipeline. Internally the pipeline computes on the index-keyed [`EdgeList`];
//! this module is the typed boundary the CLI, tests, and any embedder consume,
//! so orientation is an `enum` (not a magic `f64`) and missing values are
//! `Option<f64>` (not `NaN`).

use crate::scalenet::{EdgeList, Method};
use crate::scs::RawEdge;
use crate::util::pair_at;

/// Edge direction on the ordered pair `(source, target)` with `source < target`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation {
    /// `source -> target`.
    Forward,
    /// `target -> source`.
    Backward,
    /// Present but unoriented.
    Undirected,
}

/// Whether a network is one method's reconstruction or the fused consensus.
/// Only affects the numeric orientation code used in the R-compatible output.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum NetworkKind {
    PerMethod,
    Consensus,
}

impl Orientation {
    /// Per-method networks encode orientation as `2 / -2 / 1`; the consensus as
    /// `1 / -1 / 0`.
    fn from_code(code: f64, kind: NetworkKind) -> Orientation {
        match kind {
            NetworkKind::PerMethod => match code as i64 {
                2 => Orientation::Forward,
                -2 => Orientation::Backward,
                _ => Orientation::Undirected,
            },
            NetworkKind::Consensus => {
                if code > 0.0 {
                    Orientation::Forward
                } else if code < 0.0 {
                    Orientation::Backward
                } else {
                    Orientation::Undirected
                }
            }
        }
    }

    /// Numeric code used in the R-compatible `eorient` column for a given kind.
    pub fn code(self, kind: NetworkKind) -> f64 {
        match (kind, self) {
            (NetworkKind::PerMethod, Orientation::Forward) => 2.0,
            (NetworkKind::PerMethod, Orientation::Backward) => -2.0,
            (NetworkKind::PerMethod, Orientation::Undirected) => 1.0,
            (NetworkKind::Consensus, Orientation::Forward) => 1.0,
            (NetworkKind::Consensus, Orientation::Backward) => -1.0,
            (NetworkKind::Consensus, Orientation::Undirected) => 0.0,
        }
    }
}

/// One present edge. `source`/`target` index into [`PipelineResult::nodes`].
#[derive(Clone, Debug)]
pub struct Edge {
    pub source: usize,
    pub target: usize,
    pub orientation: Orientation,
    /// `epresenceScore` — ARACNE presence score, or the consensus average rank.
    pub presence_score: Option<f64>,
    /// `eorientScore` — set only on consensus edges.
    pub orient_score: Option<f64>,
    /// Signed Spearman correlation of the two features.
    pub ecorr: Option<f64>,
}

/// A reconstructed network: the present edges plus whether it is a per-method or
/// consensus network.
#[derive(Clone, Debug)]
pub struct Network {
    pub kind: NetworkKind,
    pub edges: Vec<Edge>,
}

/// A consensus edge in the ranked "rawAvg" ordering (best first).
#[derive(Clone, Debug)]
pub struct RankedEdge {
    pub source: usize,
    pub target: usize,
    pub avg_rank: f64,
    pub avg_ort: f64,
    pub ecorr: Option<f64>,
}

/// The consensus for one threshold: the fused network and its ranked edge list.
#[derive(Clone, Debug)]
pub struct Consensus {
    pub network: Network,
    pub ranked: Vec<RankedEdge>,
}

/// Spectral-embedding summary (surfaced for reporting / UIs).
#[derive(Clone, Copy, Debug)]
pub struct SpectralStats {
    pub k: usize,
    pub subset_m: usize,
    pub n_subgraphs: usize,
}

/// The full typed result of [`crate::run`].
#[derive(Clone, Debug)]
pub struct PipelineResult {
    /// Node (feature) names; edges reference these by index.
    pub nodes: Vec<String>,
    pub methods: Vec<Method>,
    pub thresholds: Vec<f64>,
    /// `per_method[method_i][threshold_i]`.
    pub per_method: Vec<Vec<Network>>,
    /// `consensus[threshold_i]`.
    pub consensus: Vec<Consensus>,
    pub spectral: SpectralStats,
    /// The spectral subgraphs that were reconstructed and fused into the network:
    /// each is a list of node indices (`nodes`), in order `pos, neg` per
    /// eigenvector 1..k. Surfaced so a UI can show "the subnetworks that made it".
    pub subgraphs: Vec<Vec<usize>>,
}

fn opt(x: f64) -> Option<f64> {
    if x.is_nan() {
        None
    } else {
        Some(x)
    }
}

impl Network {
    /// Convert an internal index-keyed edge list into a typed network. Edges are
    /// sorted by `(source, target)` for a deterministic public ordering.
    pub(crate) fn from_edge_list(el: &EdgeList, kind: NetworkKind) -> Network {
        let mut edges: Vec<Edge> = el
            .edges
            .iter()
            .map(|(&idx, rec)| {
                let (source, target) = pair_at(idx, el.p);
                Edge {
                    source,
                    target,
                    orientation: Orientation::from_code(rec.eorient, kind),
                    presence_score: opt(rec.epresence_score),
                    orient_score: opt(rec.eorient_score),
                    ecorr: opt(rec.ecorr),
                }
            })
            .collect();
        edges.sort_by_key(|e| (e.source, e.target));
        Network { kind, edges }
    }
}

impl RankedEdge {
    pub(crate) fn from_raw(raw: &RawEdge, p: usize) -> RankedEdge {
        let (source, target) = pair_at(raw.idx, p);
        RankedEdge {
            source,
            target,
            avg_rank: raw.avg_rank,
            avg_ort: raw.avg_ort,
            ecorr: opt(raw.ecorr),
        }
    }
}
