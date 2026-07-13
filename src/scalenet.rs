//! Per-method ScaleNet pipeline: spectral subsets -> local reconstruction ->
//! gather across subgraphs -> per-threshold edge lists.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::aracne;
use crate::data::Dataset;
use crate::hc;
use crate::recon::SubEdge;
use crate::spectral::{self, Spectral};
use crate::util::{npairs, pair_index};

/// Reconstruction method selection.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Method {
    Aracne,
    BayesHc,
}

impl Method {
    pub fn name(self) -> &'static str {
        match self {
            Method::Aracne => "aracne",
            Method::BayesHc => "bayes_hc",
        }
    }
}

/// Tunables carried through the pipeline.
#[derive(Clone)]
pub struct Params {
    pub var_perc: f64,
    pub eigen_perc: Option<f64>,
    pub pres_freq_thresh: Vec<f64>,
    // aracne
    pub epsilon: f64,
    // bayes_hc
    pub restart: usize,
    pub seed: u64,
}

/// The per-pair values for one present edge (`NaN` = NA). Presence is implicit:
/// a pair is present iff it has a record in the [`EdgeList`].
#[derive(Clone, Default)]
pub struct EdgeRec {
    pub eorient: f64,
    pub epresence_score: f64,
    pub eorient_score: f64, // only set by consensus
    pub ecorr: f64,
}

/// A scaleNet-style edge table, stored **sparsely**: only present edges are
/// materialized, keyed by their `combn` pair index. This keeps memory O(#edges)
/// instead of O(p²) while producing the same full-pair output (absent pairs are
/// emitted as NA on the fly).
pub struct EdgeList {
    pub p: usize,
    pub edges: HashMap<usize, EdgeRec>,
}

impl EdgeList {
    pub fn new(p: usize) -> EdgeList {
        EdgeList {
            p,
            edges: HashMap::new(),
        }
    }
    /// Number of present edges.
    pub fn len(&self) -> usize {
        self.edges.len()
    }
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }
    /// Record for `idx`, if the pair is present.
    pub fn get(&self, idx: usize) -> Option<&EdgeRec> {
        self.edges.get(&idx)
    }
}

/// Shared preprocessing (identical across reconstruction methods).
pub struct Prep {
    pub spectral: Spectral,
    pub subset_m: usize,
    /// Subgraphs: each carries the ascending global column indices it spans.
    pub subgraphs: Vec<Vec<usize>>,
}

/// Compute the spectral embedding and derive the spectral subgraphs.
pub fn prepare(data: &Dataset, params: &Params) -> Result<Prep, String> {
    let ncol = data.ncol();
    let subset_m = (params.var_perc * ncol as f64).ceil() as usize;

    let cols: Vec<&[u32]> = data.cols.iter().map(|c| c.as_slice()).collect();
    let sp = spectral::compute(
        &cols,
        &data.levels,
        &data.names,
        ncol,
        params.eigen_perc,
        subset_m,
    )?;

    // Map surviving node name -> global column index.
    let name_to_global: HashMap<&str, usize> = data
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let q = sp.node_names.len();
    let m = subset_m.min(q);
    let mut subgraphs = Vec::new();
    // Eigenvectors 2..k (skip the trivial smallest one at column 0).
    for c in 1..sp.k {
        // (surviving-row-index, value)
        let mut order: Vec<usize> = (0..q).collect();
        order.sort_by(|&i, &j| {
            sp.eigenvectors[j][c].total_cmp(&sp.eigenvectors[i][c])
        }); // descending by value

        let top: Vec<usize> = order[..m]
            .iter()
            .map(|&r| name_to_global[sp.node_names[r].as_str()])
            .collect();
        let bottom: Vec<usize> = order[q - m..]
            .iter()
            .map(|&r| name_to_global[sp.node_names[r].as_str()])
            .collect();

        let mut pos = top;
        pos.sort_unstable(); // preserve original column order
        let mut neg = bottom;
        neg.sort_unstable();
        subgraphs.push(pos);
        subgraphs.push(neg);
    }

    Ok(Prep {
        spectral: sp,
        subset_m,
        subgraphs,
    })
}

/// Run one reconstruction method over all subgraphs and gather the results.
///
/// Returns one `EdgeList` per presence-frequency threshold (same order as
/// `params.pres_freq_thresh`).
pub fn run_method(
    data: &Dataset,
    prep: &Prep,
    method: Method,
    params: &Params,
) -> Vec<EdgeList> {
    let p = data.ncol();

    // Reconstruct every subgraph (embarrassingly parallel).
    let results: Vec<Vec<SubEdge>> = prep
        .subgraphs
        .par_iter()
        .enumerate()
        .map(|(idx, sub)| match method {
            Method::Aracne => aracne::reconstruct(sub, &prep.spectral.mim, params.epsilon),
            Method::BayesHc => hc::reconstruct(
                sub,
                &data.cols,
                &data.levels,
                params.restart,
                params.seed.wrapping_add(idx as u64),
            ),
        })
        .collect();

    // Gather: presence / forward / backward / possible counts per pair. These
    // are dense u32 arrays because, with overlapping subgraphs, most pairs are
    // "possible" — dense arrays are both faster and smaller than a hash map at
    // this coverage. They are transient; the persistent per-threshold edge lists
    // built below are sparse.
    let n = npairs(p);
    let mut presence = vec![0u32; n];
    let mut forward = vec![0u32; n];
    let mut backward = vec![0u32; n];
    let mut possible = vec![0u32; n];
    let mut scores: HashMap<usize, Vec<f64>> = HashMap::new();

    for (sub, edges) in prep.subgraphs.iter().zip(results.iter()) {
        for i in 0..sub.len() {
            for j in (i + 1)..sub.len() {
                possible[pair_index(sub[i], sub[j], p)] += 1;
            }
        }
        for e in edges {
            let idx = pair_index(e.a, e.b, p);
            presence[idx] += 1;
            match e.eorient {
                2 => forward[idx] += 1,
                -2 => backward[idx] += 1,
                _ => {}
            }
            if !e.score.is_nan() {
                scores.entry(idx).or_default().push(e.score);
            }
        }
    }

    // Averaged presence score (mean over unique observed values), like convert.
    let mut avg_score: HashMap<usize, f64> = HashMap::with_capacity(scores.len());
    for (idx, vals) in scores {
        let mut uniq = vals;
        uniq.sort_by(|x, y| x.total_cmp(y));
        uniq.dedup();
        avg_score.insert(idx, uniq.iter().sum::<f64>() / uniq.len() as f64);
    }

    // Build one sparse edge list per threshold: only pairs with freq >= thr.
    params
        .pres_freq_thresh
        .iter()
        .map(|&thr| {
            let mut el = EdgeList::new(p);
            for idx in 0..n {
                if possible[idx] == 0 {
                    continue;
                }
                let freq = presence[idx] as f64 / possible[idx] as f64;
                if freq < thr {
                    continue;
                }
                let pres_ort = if forward[idx] > backward[idx] {
                    2.0
                } else if backward[idx] > forward[idx] {
                    -2.0
                } else {
                    1.0
                };
                el.edges.insert(
                    idx,
                    EdgeRec {
                        eorient: pres_ort, // 1 / 2 / -2
                        epresence_score: *avg_score.get(&idx).unwrap_or(&f64::NAN),
                        eorient_score: f64::NAN,
                        ecorr: f64::NAN,
                    },
                );
            }
            el
        })
        .collect()
}

/// Fill in the Spearman-based sign (`ecorr`) for every present edge, using
/// precomputed per-column ranks (see [`Dataset::rank_matrix`]). Only present
/// edges are visited and their correlations are computed in parallel.
pub fn fill_ecorr(data: &Dataset, ranks: &[Vec<f64>], el: &mut EdgeList) {
    let p = data.ncol();
    // Pearson on precomputed ranks == Spearman on the codes; parallel over the
    // (sparse) present edges.
    el.edges.par_iter_mut().for_each(|(&idx, rec)| {
        let (a, b) = crate::util::pair_at(idx, p);
        rec.ecorr = crate::data::corr(&ranks[a], &ranks[b]);
    });
}
