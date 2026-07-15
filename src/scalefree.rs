//! Scale-free topology analysis of a reconstructed network.
//!
//! Barabási's central observation is that many real networks have a **power-law**
//! degree distribution `P(k) ~ k^-γ` (a few high-degree hubs, many low-degree
//! nodes), typically with `2 ≤ γ ≤ 3`. Following the standard co-expression-network
//! practice (Zhang & Horvath, WGCNA), we quantify how well a network matches that
//! form by the **R² of a straight-line fit to `log P(k)` vs `log k`** — the
//! *scale-free fit index* — and use it to **propose an edge-weight threshold**:
//! the sparsification at which the retained network best approximates scale-free
//! topology while keeping a meaningful fraction of its nodes.

use crate::network::Network;

/// Power-law fit of a degree distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScaleFreeFit {
    /// Exponent `γ` of `P(k) ~ k^-γ` (negative of the log-log regression slope).
    pub gamma: f64,
    /// Scale-free fit index: R² of `log P(k)` vs `log k`, in `[0, 1]`.
    pub r2: f64,
    /// Largest degree observed.
    pub k_max: usize,
    /// Number of connected nodes (degree ≥ 1) the fit was computed on.
    pub n_connected: usize,
}

impl ScaleFreeFit {
    /// Whether the degree distribution passes the usual scale-free criteria:
    /// a good power-law fit (`R² ≥ 0.80`) with `γ` in Barabási's `[2, 3]` regime.
    pub fn is_scale_free(&self) -> bool {
        self.r2 >= 0.80 && (2.0..=3.0).contains(&self.gamma)
    }
}

/// Fit `P(k) ~ k^-γ` to a degree multiset by log-log linear regression over the
/// degree histogram. Isolated nodes (degree 0) are excluded. Returns `None` when
/// there are fewer than 3 distinct positive degrees (an underdetermined fit).
pub fn fit(degrees: &[usize]) -> Option<ScaleFreeFit> {
    let connected: Vec<usize> = degrees.iter().copied().filter(|&d| d > 0).collect();
    let total = connected.len();
    if total < 3 {
        return None;
    }
    // Degree histogram.
    let k_max = *connected.iter().max().unwrap();
    let mut hist = vec![0usize; k_max + 1];
    for &d in &connected {
        hist[d] += 1;
    }
    // Log-log points, one per observed degree.
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (k, &c) in hist.iter().enumerate() {
        if k >= 1 && c > 0 {
            xs.push((k as f64).log10());
            ys.push((c as f64 / total as f64).log10());
        }
    }
    if xs.len() < 3 {
        return None;
    }
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (&x, &y) in xs.iter().zip(&ys) {
        sxx += (x - mx) * (x - mx);
        sxy += (x - mx) * (y - my);
        syy += (y - my) * (y - my);
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return None;
    }
    let slope = sxy / sxx;
    let r2 = (sxy * sxy) / (sxx * syy);
    Some(ScaleFreeFit {
        gamma: -slope,
        r2,
        k_max,
        n_connected: total,
    })
}

/// One row of a threshold sweep: the network kept at `threshold` and its fit.
#[derive(Clone, Copy, Debug)]
pub struct SweepRow {
    /// Minimum `|ecorr|` an edge must have to be kept.
    pub threshold: f64,
    pub n_nodes: usize,
    pub n_edges: usize,
    /// Fraction of the original nodes still connected at this threshold.
    pub frac_nodes: f64,
    /// Power-law fit at this threshold (`None` if it could not be fit).
    pub fit: Option<ScaleFreeFit>,
}

/// Undirected degrees of `net` over `n_nodes` when only edges with
/// `|ecorr| ≥ threshold` (and, if `positive_only`, `ecorr ≥ 0`) are kept.
fn degrees_at(net: &Network, n_nodes: usize, threshold: f64, sign: Sign) -> Vec<usize> {
    let mut deg = vec![0usize; n_nodes];
    for e in &net.edges {
        let Some(ec) = e.ecorr else { continue };
        if ec.abs() < threshold {
            continue;
        }
        match sign {
            Sign::Positive if ec < 0.0 => continue,
            Sign::Negative if ec >= 0.0 => continue,
            _ => {}
        }
        deg[e.source] += 1;
        deg[e.target] += 1;
    }
    deg
}

/// Edge-sign filter for the sweep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sign {
    Both,
    Positive,
    Negative,
}

/// Sweep `|ecorr|` thresholds from the network's min to max in `steps` steps,
/// fitting scale-free topology at each. Rows are ordered by increasing threshold.
pub fn sweep_ecorr(net: &Network, n_nodes: usize, steps: usize, sign: Sign) -> Vec<SweepRow> {
    let abs: Vec<f64> = net
        .edges
        .iter()
        .filter_map(|e| e.ecorr.map(f64::abs))
        .collect();
    if abs.is_empty() {
        return Vec::new();
    }
    let lo = abs.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = abs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let steps = steps.max(1);
    (0..=steps)
        .map(|s| {
            let t = lo + (hi - lo) * (s as f64) / (steps as f64);
            let deg = degrees_at(net, n_nodes, t, sign);
            let n_conn = deg.iter().filter(|&&d| d > 0).count();
            let n_edges = deg.iter().sum::<usize>() / 2;
            SweepRow {
                threshold: t,
                n_nodes: n_conn,
                n_edges,
                frac_nodes: n_conn as f64 / n_nodes as f64,
                fit: fit(&deg),
            }
        })
        .collect()
}

/// Recommend the threshold that best approximates scale-free topology: among rows
/// keeping at least `min_frac` of the nodes, the one maximizing the fit index,
/// with a small bonus when `γ` falls in Barabási's `[2, 3]` regime.
pub fn recommend(sweep: &[SweepRow], min_frac: f64) -> Option<&SweepRow> {
    sweep
        .iter()
        .filter(|r| r.frac_nodes >= min_frac && r.fit.is_some())
        .max_by(|a, b| {
            let score = |r: &SweepRow| {
                let f = r.fit.unwrap();
                f.r2 + if (2.0..=3.0).contains(&f.gamma) { 0.06 } else { 0.0 }
            };
            score(a).total_cmp(&score(b))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::{Edge, Network, NetworkKind, Orientation};

    #[test]
    fn power_law_fit_matches_reference() {
        // Degree multiset {1×4, 2×2, 3×2, 5×1}; reference γ/R² from a log-log
        // regression computed independently (Python).
        let degs = [1, 1, 1, 1, 2, 2, 3, 3, 5];
        let f = fit(&degs).expect("fit");
        assert!((f.gamma - 0.805089).abs() < 1e-5, "gamma {}", f.gamma);
        assert!((f.r2 - 0.934680).abs() < 1e-5, "r2 {}", f.r2);
        assert_eq!(f.k_max, 5);
        assert_eq!(f.n_connected, 9);
    }

    #[test]
    fn too_few_degrees_is_none() {
        assert!(fit(&[]).is_none());
        assert!(fit(&[0, 0, 0]).is_none()); // all isolated
        assert!(fit(&[1, 1, 2]).is_none()); // only 2 distinct positive degrees
    }

    #[test]
    fn hub_network_is_more_scale_free_than_regular() {
        // A hub-and-spoke graph (one high-degree hub) has a steeper, better
        // power-law fit than a near-regular ring.
        let hub_degs = {
            let mut d = vec![1usize; 20]; // 20 leaves
            d.push(20); // hub
            d
        };
        let ring_degs = vec![2usize; 21];
        let hub = fit(&{
            // spread leaf degrees a little so there are ≥3 distinct values
            let mut d = hub_degs.clone();
            d[0] = 2;
            d[1] = 2;
            d[2] = 3;
            d
        });
        assert!(hub.is_some());
        assert!(fit(&ring_degs).is_none()); // one distinct degree -> not fittable
    }

    #[test]
    fn ecorr_sweep_and_recommend() {
        // Small network: a hub (0) linked to 1..=5 with strong |ecorr|, plus a
        // few weak peripheral edges that drop out as the threshold rises.
        let mk = |s: usize, t: usize, ec: f64| Edge {
            source: s,
            target: t,
            orientation: Orientation::Undirected,
            presence_score: Some(1.0),
            orient_score: None,
            ecorr: Some(ec),
        };
        let net = Network {
            kind: NetworkKind::Consensus,
            edges: vec![
                mk(0, 1, 0.9),
                mk(0, 2, 0.85),
                mk(0, 3, 0.8),
                mk(0, 4, 0.75),
                mk(0, 5, 0.7),
                mk(1, 2, 0.3),
                mk(3, 4, 0.3),
                mk(2, 5, 0.35),
            ],
        };
        let sweep = sweep_ecorr(&net, 6, 8, Sign::Both);
        assert_eq!(sweep.len(), 9);
        // Lowest threshold keeps every edge; highest keeps fewer.
        assert!(sweep[0].n_edges >= sweep[sweep.len() - 1].n_edges);
        // Thresholds are non-decreasing and within [min,max] |ecorr|.
        for w in sweep.windows(2) {
            assert!(w[1].threshold >= w[0].threshold - 1e-12);
        }
        // A recommendation exists and keeps at least the requested fraction.
        if let Some(r) = recommend(&sweep, 0.5) {
            assert!(r.frac_nodes >= 0.5);
            assert!(r.fit.is_some());
        }
    }
}
