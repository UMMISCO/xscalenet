//! Per-method ScaleNet pipeline: spectral subsets -> local reconstruction ->
//! gather across subgraphs -> per-threshold edge lists.

use std::collections::HashMap;

#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(not(feature = "parallel"))]
use crate::par::prelude::*;

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
    pub max_iter: usize,
    pub iss: f64,
    pub maxp: usize,
    /// Use the iterative (Lanczos) eigensolver for the spectral embedding — a
    /// fast path when `eigen_perc` is fixed. Off by default (dense, bit-identical).
    pub iterative_eigen: bool,
    /// Reproduce the pre-fix, **column-order-dependent** pipeline: decompose the
    /// lower triangle of the non-symmetric `Lrw` (as R does), skip eigenvector
    /// sign canonicalisation, and keep the index-order tie-breaks in the subset
    /// selection and in hill-climbing. Only for reproducing results published
    /// before the fix; see [`crate::spectral`].
    pub legacy_lrw: bool,
    /// Number of random variable orderings the Bayesian hill climbing is run
    /// under per subgraph, whose results are then fused (see
    /// [`crate::hc::reconstruct`]). `1` = a single ordering, which makes the
    /// network a function of that arbitrary order. Ignored by ARACNE, which is
    /// order-invariant already.
    pub hc_orders: usize,
    /// Fraction of those orderings an arc must appear in to be kept.
    pub hc_order_freq: f64,
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
    /// How many eigenvectors had their `G+`/`G-` cut fall inside a block of
    /// numerically-indistinguishable components. For those, which variables
    /// enter the subset is arbitrary at the available precision, so the subset —
    /// and the edges derived from it — are not reproducible across input column
    /// orders. Caused by duplicate or near-duplicate input variables.
    pub boundary_ties: usize,
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
        params.iterative_eigen,
        params.seed,
        params.legacy_lrw,
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
    // Eigenvectors whose G+/G- cut falls inside a numerically-tied block.
    let mut boundary_ties = 0usize;
    // Eigenvectors 2..k (skip the trivial smallest one at column 0).
    for c in 1..sp.k {
        // (surviving-row-index, value)
        let mut order: Vec<usize> = (0..q).collect();
        order.sort_by(|&i, &j| {
            // Descending by eigenvector value. An exact tie must not fall back
            // to the surviving-row index, which *is* the input column order, so
            // break it on the node name instead. Values that merely sit close
            // together are left in value order: they are reproducible to ~1e-14
            // and reordering them on a tolerance would replace a stable
            // comparison with a noise-sensitive grouping decision (measured:
            // doing so costs eDNA_own its exact invariance). Where the boundary
            // genuinely falls inside a tied block the selection is ill-posed,
            // and `boundary_ties` reports it rather than hiding it.
            let by_value = sp.eigenvectors[j][c].total_cmp(&sp.eigenvectors[i][c]);
            if params.legacy_lrw {
                by_value
            } else {
                by_value.then_with(|| sp.node_names[i].cmp(&sp.node_names[j]))
            }
        });

        // Is either cut sitting inside a block of numerically-equal values? Then
        // which variables land in the subset is arbitrary at the available
        // precision, and the subset is not reproducible across column orders.
        if !params.legacy_lrw {
            let scale = (0..q)
                .map(|r| sp.eigenvectors[r][c].abs())
                .fold(0.0f64, f64::max);
            let tol = scale * 1e-12;
            let val = |r: usize| sp.eigenvectors[order[r]][c];
            if (val(m - 1) - val(m)).abs() <= tol || (val(q - m - 1) - val(q - m)).abs() <= tol {
                boundary_ties += 1;
            }
        }

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
        boundary_ties,
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
            Method::BayesHc => {
                // Hill climbing scans candidate moves in *local* index order and
                // keeps the first best one, and its restart perturbations draw
                // local indices. With the subgraph in ascending global-column
                // order, both depend on the input column order. Ordering the
                // subgraph by node name makes the whole search a function of the
                // variable *set* alone.
                let mut ordered = sub.clone();
                if !params.legacy_lrw {
                    ordered.sort_by(|&x, &y| data.names[x].cmp(&data.names[y]));
                }
                // The subgraph index is itself order- and name-invariant (it
                // comes from the eigenvector rank and the canonical G+/G- sign),
                // so seeding the orderings from it keeps them reproducible
                // without reintroducing a dependence on the input layout.
                let orders = if params.legacy_lrw { 1 } else { params.hc_orders };
                hc::reconstruct(
                    &ordered,
                    &data.cols,
                    &data.levels,
                    params.restart,
                    params.seed.wrapping_add(idx as u64),
                    params.max_iter,
                    params.iss,
                    params.maxp,
                    orders,
                    params.hc_order_freq,
                )
            }
        })
        .collect();

    // Gather: presence / forward / backward votes and co-selection counts per
    // pair. These are dense arrays because, with overlapping subgraphs, most
    // pairs are "possible" — dense arrays are both faster and smaller than a
    // hash map at this coverage. They are transient; the persistent
    // per-threshold edge lists built below are sparse.
    //
    // The vote arrays are `f64` because a method may vote fractionally (hill
    // climbing under ordering consensus weights an arc by the share of variable
    // orderings that inferred it). A method that votes 1.0 per subgraph, like
    // ARACNE, yields exact integers here and so the same frequencies as before.
    let n = npairs(p);
    let mut presence = vec![0.0f64; n];
    let mut forward = vec![0.0f64; n];
    let mut backward = vec![0.0f64; n];
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
            presence[idx] += e.weight;
            match e.eorient {
                2 => forward[idx] += e.weight,
                -2 => backward[idx] += e.weight,
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
                let freq = presence[idx] / possible[idx] as f64;
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
