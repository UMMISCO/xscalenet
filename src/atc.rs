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

/// ATC parameters (defaults mirror the R `dynamicPropagation`).
#[derive(Clone, Copy, Debug)]
pub struct AtcConfig {
    /// Use edge weights (`|ecorr|` etc., from [`Graph::from_network`]) or 1/0.
    pub weighted: bool,
    /// Include a self-loop in the propagation (R uses `diag(J) <- 1`).
    pub self_loop: bool,
    /// Maximum propagation steps.
    pub max_iter: usize,
    /// Convergence threshold on co-inertia (RV) of consecutive configurations.
    pub coinertia_threshold: f64,
    /// Record every intermediate configuration `Tᵖ` (for step-by-step
    /// visualization). Off by default so the large-graph hot path pays nothing.
    pub capture_trace: bool,
}

impl Default for AtcConfig {
    fn default() -> Self {
        AtcConfig {
            weighted: false,
            self_loop: true,
            max_iter: 100,
            coinertia_threshold: 0.9,
            capture_trace: false,
        }
    }
}

/// A step-by-step record of the propagation, for visualization. Present only when
/// [`AtcConfig::capture_trace`] is set.
#[derive(Clone, Debug)]
pub struct AtcTrace {
    /// The column-normalized configurations `T⁰, T¹, … , Tᵖ` (each is
    /// `n_nodes × n_themes`). Length is `n_iterations + 1`.
    pub configs: Vec<Vec<Vec<f64>>>,
    /// Co-inertia (RV) between config `p-1` and `p`. Length `n_iterations`
    /// (aligned with `configs[1..]`).
    pub coinertia: Vec<f64>,
    /// Per-node L2 change `‖rowᵢ(Tᵖ) − rowᵢ(Tᵖ⁻¹)‖₂` at each step. Length
    /// `n_iterations` (aligned with `configs[1..]`).
    pub delta: Vec<Vec<f64>>,
}

/// Result of an ATC run.
#[derive(Clone, Debug)]
pub struct AtcResult {
    /// ATC score per node: `Σ_p ‖Δrowᵢ‖·(1−coinertiaₚ)`, normalized to sum 1
    /// (matches R `computeMTCDistCoinSum`).
    pub atc: Vec<f64>,
    /// Alternative ATC = row sums of the converged configuration
    /// (R `computeATC_last_config`, which the authors found often preferable).
    pub atc_last_config: Vec<f64>,
    /// Number of propagation steps performed.
    pub n_iterations: usize,
    /// Co-inertia (RV) of each consecutive pair of configurations.
    pub coinertia: Vec<f64>,
    /// Full per-step trace (`None` unless [`AtcConfig::capture_trace`] was set).
    pub trace: Option<AtcTrace>,
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

    // T⁰ is column-normalized before propagation (R normalizes each theme to 1).
    let mut prev: Vec<Vec<f64>> = t0.to_vec();
    normalize_columns(&mut prev, m);

    let mut atc = vec![0.0f64; n];
    let mut coinertia = Vec::new();
    let mut iters = 0usize;

    // Trace buffers (only populated when requested). `configs` starts with T⁰.
    let mut trace_configs: Vec<Vec<Vec<f64>>> = Vec::new();
    let mut trace_delta: Vec<Vec<f64>> = Vec::new();
    if cfg.capture_trace {
        trace_configs.push(prev.clone());
    }

    for _ in 0..cfg.max_iter {
        // Propagate: rowᵢ = w·prev[i] (self-loop) + Σ_{(u,w)∈adj[i]} w·prev[u].
        let mut cur: Vec<Vec<f64>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut row = vec![0.0f64; m];
                if cfg.self_loop {
                    let s = &prev[i];
                    for j in 0..m {
                        row[j] += s[j];
                    }
                }
                for &(u, w) in &graph.adj[i] {
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
        // Per-node L2 change ‖rowᵢ(cur) − rowᵢ(prev)‖₂.
        let dist: Vec<f64> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut d = 0.0;
                for j in 0..m {
                    let diff = cur[i][j] - prev[i][j];
                    d += diff * diff;
                }
                d.sqrt()
            })
            .collect();
        // Accumulate ATCᵢ += ‖Δrowᵢ‖ · C_p.
        for (a, &d) in atc.iter_mut().zip(&dist) {
            *a += d * c_p;
        }

        if cfg.capture_trace {
            trace_configs.push(cur.clone());
            trace_delta.push(dist);
        }

        coinertia.push(rv);
        iters += 1;
        prev = cur;
        if rv >= cfg.coinertia_threshold {
            break;
        }
    }

    // Normalize the dist-coinertia ATC to sum 1 (R `computeMTCDistCoinSum`).
    let s: f64 = atc.iter().sum();
    if s > 0.0 {
        for a in &mut atc {
            *a /= s;
        }
    }
    // Alternative ATC = row sums of the converged configuration.
    let atc_last_config: Vec<f64> = prev.iter().map(|row| row.iter().sum()).collect();

    let trace = if cfg.capture_trace {
        Some(AtcTrace {
            configs: trace_configs,
            coinertia: coinertia.clone(),
            delta: trace_delta,
        })
    } else {
        None
    };

    AtcResult {
        atc,
        atc_last_config,
        n_iterations: iters,
        coinertia,
        trace,
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

/// Column-standardize a config (center to mean 0, scale to unit population sd),
/// matching `ade4::dudi.pca(scale = TRUE)`. Constant columns become 0.
fn standardize(t: &[Vec<f64>], m: usize) -> Vec<Vec<f64>> {
    let n = t.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;
    let mut mean = vec![0.0f64; m];
    for row in t {
        for j in 0..m {
            mean[j] += row[j];
        }
    }
    for v in &mut mean {
        *v /= nf;
    }
    let mut sd = vec![0.0f64; m];
    for row in t {
        for j in 0..m {
            let d = row[j] - mean[j];
            sd[j] += d * d;
        }
    }
    for v in &mut sd {
        *v = (*v / nf).sqrt();
    }
    t.iter()
        .map(|row| {
            (0..m)
                .map(|j| if sd[j] > 0.0 { (row[j] - mean[j]) / sd[j] } else { 0.0 })
                .collect()
        })
        .collect()
}

/// Co-inertia (RV coefficient) of two configurations, on column-standardized
/// tables: `‖XᵀY‖²_F / (‖XᵀX‖_F · ‖YᵀY‖_F)`, in [0,1].
fn rv_coefficient(x_raw: &[Vec<f64>], y_raw: &[Vec<f64>], m: usize) -> f64 {
    let x = standardize(x_raw, m);
    let y = standardize(y_raw, m);
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
    fn matches_r_dynamic_propagation() {
        // The exact 10-node graph J and membership m0 from the R integranet
        // `Propagation_and_ATC_code.R`; atc_last_config must match the R output
        // (co-inertia-independent, so an exact reference).
        let g = from_edges(
            10,
            &[
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
            ],
        );
        // m0: brown, black, yellow
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
        let cfg = AtcConfig {
            weighted: false,
            self_loop: true,
            coinertia_threshold: 2.0, // never converge -> run all steps
            max_iter: 49,             // R stores T0..T49 (49 propagations)
            capture_trace: true,
        };
        let res = compute(&g, &m0, &cfg);
        // Trace has T⁰ plus one config per step; deltas/coinertia align with steps.
        let tr = res.trace.as_ref().expect("trace requested");
        assert_eq!(tr.configs.len(), res.n_iterations + 1);
        assert_eq!(tr.delta.len(), res.n_iterations);
        assert_eq!(tr.coinertia.len(), res.n_iterations);
        // Final config in the trace equals atc_last_config's source rows.
        let last = tr.configs.last().unwrap();
        for i in 0..10 {
            let rowsum: f64 = last[i].iter().sum();
            assert!((rowsum - res.atc_last_config[i]).abs() < 1e-9);
        }
        let r_ref = [
            0.4028144, 0.5090365, 0.6156974, 0.0, 0.4028144, 0.4044723, 0.2421126, 0.2715327,
            0.1115623, 0.0399573,
        ];
        for i in 0..10 {
            assert!(
                (res.atc_last_config[i] - r_ref[i]).abs() < 1e-5,
                "node {i}: rust {} vs R {}",
                res.atc_last_config[i],
                r_ref[i]
            );
        }

        // Co-inertia (RV) vs ade4::coinertia on the same configs.
        let r_coin = [0.513602, 0.918923, 0.972564, 0.989206, 0.991419];
        for (p, &rc) in r_coin.iter().enumerate() {
            assert!(
                (res.coinertia[p] - rc).abs() < 1e-4,
                "coinertia step {p}: rust {} vs R {rc}",
                res.coinertia[p]
            );
        }
        // Dist-coinertia ATC (normalized) vs R computeMTCDistCoinSum.
        let r_atc = [
            0.1127189, 0.1472663, 0.1467608, 0.1469249, 0.1479314, 0.0871726, 0.1080605,
            0.0526891, 0.0095518, 0.0409239,
        ];
        for i in 0..10 {
            assert!(
                (res.atc[i] - r_atc[i]).abs() < 1e-4,
                "atc node {i}: rust {} vs R {}",
                res.atc[i],
                r_atc[i]
            );
        }
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
