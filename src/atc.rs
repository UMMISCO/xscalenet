//! Annotation Transcriptional Centrality (ATC) — propagation of node annotations
//! over the network, following Prifti et al., *Bioinformatics* 26(24):3083 (2010).
//!
//! Given a network and a node×theme membership matrix `T⁰` (`W⁰_ij = 1` if node
//! *i* carries annotating theme *j*), a non-linear dynamical system propagates the
//! annotations:
//!
//! * `Tᵖ = A · Tᵖ⁻¹` — each node's theme weights become the (edge-weighted) sum of
//!   its neighbours' weights from the previous step;
//! * each theme **column of `Tᵖ` is normalized to sum 1** (constant total weight
//!   per theme);
//! * the system is considered convergent when the **co-inertia** of two
//!   consecutive configurations exceeds a threshold (default 0.9).
//!
//! The ATC score of a node aggregates, over the propagation, its local change
//! weighted by the global change:
//!
//! * `D_{i,p} = ‖rowᵢ(Tᵖ) − rowᵢ(Tᵖ⁻¹)‖₂`  (Eq. 1),
//! * `C_p = 1 − co-inertia(Tᵖ⁻¹, Tᵖ)`,
//! * `ATCᵢ = Σ_p D_{i,p} · C_p`.
//!
//! Co-inertia is the **RV coefficient** of the two configuration matrices
//! (`RV(X,Y) = ‖XᵀY‖²_F / (‖XᵀX‖_F · ‖YᵀY‖_F)`), a standard normalized
//! co-inertia measure; the paper's ade4 co-inertia may differ in
//! centering/weighting, so treat these scores as faithful-in-spirit.

use rayon::prelude::*;

use crate::metrics::Graph;

/// ATC parameters.
#[derive(Clone, Copy, Debug)]
pub struct AtcConfig {
    /// Use edge weights (`|ecorr|` etc., from [`Graph::from_network`]) or 1/0.
    pub weighted: bool,
    /// Maximum propagation steps.
    pub max_iter: usize,
    /// Convergence threshold on co-inertia of consecutive configurations.
    pub coinertia_threshold: f64,
}

impl Default for AtcConfig {
    fn default() -> Self {
        AtcConfig {
            weighted: false,
            max_iter: 100,
            coinertia_threshold: 0.9,
        }
    }
}

/// Result of an ATC run.
#[derive(Clone, Debug)]
pub struct AtcResult {
    /// ATC score per node (higher = more influential on theme propagation).
    pub atc: Vec<f64>,
    /// Number of propagation steps performed.
    pub n_iterations: usize,
    /// Co-inertia of each consecutive pair of configurations.
    pub coinertia: Vec<f64>,
}

/// Build a node×theme membership matrix `T⁰` from a single categorical label per
/// node (`None` = unannotated). `W⁰_ij = 1` iff node *i*'s label is theme *j*.
pub fn membership_from_labels(labels: &[Option<usize>], n_themes: usize) -> Vec<Vec<f64>> {
    labels
        .iter()
        .map(|&l| {
            let mut row = vec![0.0f64; n_themes];
            if let Some(t) = l {
                if t < n_themes {
                    row[t] = 1.0;
                }
            }
            row
        })
        .collect()
}

/// Compute ATC scores. `t0` is the node×theme membership matrix (`graph.n` rows).
pub fn compute(graph: &Graph, t0: &[Vec<f64>], cfg: &AtcConfig) -> AtcResult {
    let n = graph.n;
    assert_eq!(t0.len(), n, "t0 must have one row per node");
    let m = t0.first().map(|r| r.len()).unwrap_or(0);

    let mut prev: Vec<Vec<f64>> = t0.to_vec();
    let mut atc = vec![0.0f64; n];
    let mut coinertia = Vec::new();
    let mut iters = 0usize;

    for _ in 0..cfg.max_iter {
        // Propagate: rowᵢ = Σ_{(u,w)∈adj[i]} w · prev[u]   (weighted per config).
        let mut cur: Vec<Vec<f64>> = graph
            .adj
            .par_iter()
            .map(|neigh| {
                let mut row = vec![0.0f64; m];
                for &(u, w) in neigh {
                    let wt = if cfg.weighted { w } else { 1.0 };
                    let src = &prev[u];
                    for j in 0..m {
                        row[j] += wt * src[j];
                    }
                }
                row
            })
            .collect();
        normalize_columns(&mut cur, m);

        let rv = rv_coefficient(&prev, &cur, m);
        let c_p = 1.0 - rv;
        // Accumulate ATCᵢ += ‖rowᵢ(cur) − rowᵢ(prev)‖₂ · C_p.
        atc.par_iter_mut().enumerate().for_each(|(i, a)| {
            let mut d = 0.0;
            for j in 0..m {
                let diff = cur[i][j] - prev[i][j];
                d += diff * diff;
            }
            *a += d.sqrt() * c_p;
        });

        coinertia.push(rv);
        iters += 1;
        prev = cur;
        if rv >= cfg.coinertia_threshold {
            break;
        }
    }

    AtcResult {
        atc,
        n_iterations: iters,
        coinertia,
    }
}

/// Normalize each theme column to sum 1 (columns that sum to 0 are left as 0).
fn normalize_columns(t: &mut [Vec<f64>], m: usize) {
    for j in 0..m {
        let s: f64 = t.iter().map(|row| row[j]).sum();
        if s > 0.0 {
            for row in t.iter_mut() {
                row[j] /= s;
            }
        }
    }
}

/// RV coefficient of two n×m matrices: `‖XᵀY‖²_F / (‖XᵀX‖_F · ‖YᵀY‖_F)`, in [0,1].
fn rv_coefficient(x: &[Vec<f64>], y: &[Vec<f64>], m: usize) -> f64 {
    // Cross/inner product Gram matrices (m×m), computed as Σ_i outer products.
    let mut xy = vec![0.0f64; m * m];
    let mut xx = vec![0.0f64; m * m];
    let mut yy = vec![0.0f64; m * m];
    for i in 0..x.len() {
        let xr = &x[i];
        let yr = &y[i];
        for a in 0..m {
            for b in 0..m {
                xy[a * m + b] += xr[a] * yr[b];
                xx[a * m + b] += xr[a] * xr[b];
                yy[a * m + b] += yr[a] * yr[b];
            }
        }
    }
    let fro_sq = |v: &[f64]| v.iter().map(|z| z * z).sum::<f64>();
    let num = fro_sq(&xy);
    let den = fro_sq(&xx).sqrt() * fro_sq(&yy).sqrt();
    if den == 0.0 {
        0.0
    } else {
        num / den
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::Graph;

    fn from_edges(n: usize, edges: &[(usize, usize)]) -> Graph {
        let mut adj = vec![Vec::new(); n];
        for &(a, b) in edges {
            adj[a].push((b, 1.0));
            adj[b].push((a, 1.0));
        }
        Graph { n, adj }
    }

    #[test]
    fn bridge_node_has_high_atc() {
        // Two theme-homogeneous triangles {0,1,2}=theme 0, {3,4,5}=theme 1,
        // joined by bridge 2-3. The bridge mediates cross-theme propagation.
        let g = from_edges(
            6,
            &[(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)],
        );
        let labels = [Some(0), Some(0), Some(0), Some(1), Some(1), Some(1)];
        let t0 = membership_from_labels(&labels, 2);
        let res = compute(&g, &t0, &AtcConfig::default());
        assert!(res.n_iterations >= 1);
        // Co-inertia should trend upward toward convergence.
        assert!(*res.coinertia.last().unwrap() >= res.coinertia[0] - 1e-9);
        // Bridge nodes (2,3) outscore the peripheral triangle nodes (0,5).
        assert!(res.atc[2] > res.atc[0], "atc {:?}", res.atc);
        assert!(res.atc[3] > res.atc[5], "atc {:?}", res.atc);
        assert!(res.atc.iter().all(|v| v.is_finite()));
    }
}
