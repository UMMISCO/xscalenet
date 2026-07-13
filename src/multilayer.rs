//! Multi-layer **theme-proximity network** — the functional layer of the FunNet /
//! integranet model (Prifti et al. 2008; R `dynamicPropagation`).
//!
//! Node annotations are propagated over the network (weighted neighbour sum with
//! a self-loop, columns L2-normalized each step) until the configuration
//! stabilizes. The converged config gives each theme a profile over the nodes;
//! the **theme-proximity network** is the self-tuning Gaussian affinity between
//! theme profiles:
//!
//! * `d(a,b)` = Euclidean distance between theme columns `a`, `b`;
//! * local scale `σ_i = median{ d(i,·) ≥ median(d(i,·)) }` (Zelnik-Manor–Perona);
//! * `proximity(a,b) = exp(−d(a,b)² / (σ_a·σ_b))`, zero diagonal.
//!
//! Themes are the nodes of this second-layer network; edges are their proximity.

use crate::metrics::Graph;

/// Parameters for the theme-proximity construction.
#[derive(Clone, Copy, Debug)]
pub struct ThemeConfig {
    pub self_loop: bool,
    pub weighted: bool,
    pub max_iter: usize,
    /// Stop when the Frobenius change between configs drops below this.
    pub tol: f64,
    /// Consider only annotated nodes (row sum > 0) for the theme distances.
    pub annotated_only: bool,
}

impl Default for ThemeConfig {
    fn default() -> Self {
        ThemeConfig {
            self_loop: true,
            weighted: false,
            max_iter: 200,
            tol: 1e-9,
            annotated_only: true,
        }
    }
}

/// The theme-proximity (second-layer) network.
#[derive(Clone, Debug)]
pub struct ThemeProximity {
    pub themes: Vec<String>,
    /// `m × m` theme-proximity affinity (zero diagonal).
    pub proximity: Vec<Vec<f64>>,
    /// `m × m` theme-column Euclidean distances.
    pub distance: Vec<Vec<f64>>,
    pub n_iterations: usize,
}

fn l2_normalize_columns(t: &mut [Vec<f64>], m: usize) {
    for j in 0..m {
        let norm = t.iter().map(|r| r[j] * r[j]).sum::<f64>().sqrt();
        if norm > 0.0 {
            for r in t.iter_mut() {
                r[j] /= norm;
            }
        }
    }
}

fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.total_cmp(b));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// Build the theme-proximity network from a node×theme membership `t0`.
pub fn theme_proximity(
    graph: &Graph,
    t0: &[Vec<f64>],
    themes: &[String],
    cfg: &ThemeConfig,
) -> ThemeProximity {
    let n = graph.n;
    let m = t0.first().map(|r| r.len()).unwrap_or(0);

    // Propagate (L2 column-normalized) until stable.
    let mut prev = t0.to_vec();
    l2_normalize_columns(&mut prev, m);
    let mut iters = 0usize;
    for _ in 0..cfg.max_iter {
        let mut cur: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut row = vec![0.0f64; m];
                if cfg.self_loop {
                    for j in 0..m {
                        row[j] += prev[i][j];
                    }
                }
                for &(u, w) in &graph.adj[i] {
                    let wt = if cfg.weighted { w } else { 1.0 };
                    for j in 0..m {
                        row[j] += wt * prev[u][j];
                    }
                }
                row
            })
            .collect();
        l2_normalize_columns(&mut cur, m);
        let change: f64 = (0..n)
            .map(|i| (0..m).map(|j| (cur[i][j] - prev[i][j]).powi(2)).sum::<f64>())
            .sum::<f64>()
            .sqrt();
        prev = cur;
        iters += 1;
        if change < cfg.tol {
            break;
        }
    }

    // Rows (nodes) used for the theme distances.
    let rows: Vec<usize> = if cfg.annotated_only {
        (0..n).filter(|&i| prev[i].iter().sum::<f64>() > 0.0).collect()
    } else {
        (0..n).collect()
    };

    // Theme-column Euclidean distance matrix.
    let mut distance = vec![vec![0.0f64; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let d: f64 = rows
                .iter()
                .map(|&i| (prev[i][a] - prev[i][b]).powi(2))
                .sum::<f64>()
                .sqrt();
            distance[a][b] = d;
            distance[b][a] = d;
        }
    }

    // Self-tuning local sigma per theme, then Gaussian affinity.
    let sigma: Vec<f64> = (0..m)
        .map(|i| {
            let row = &distance[i];
            let med = median(row);
            let upper: Vec<f64> = row.iter().cloned().filter(|&x| x >= med).collect();
            median(&upper)
        })
        .collect();
    let mut proximity = vec![vec![0.0f64; m]; m];
    for a in 0..m {
        for b in 0..m {
            if a != b && sigma[a] > 0.0 && sigma[b] > 0.0 {
                proximity[a][b] = (-(distance[a][b].powi(2)) / (sigma[a] * sigma[b])).exp();
            }
        }
    }

    ThemeProximity {
        themes: themes.to_vec(),
        proximity,
        distance,
        n_iterations: iters,
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
    fn matches_r_theme_proximity() {
        // Same 10-node J graph + m0 as the R dynamicPropagation reference.
        let g = from_edges(
            10,
            &[
                (0, 1), (0, 2), (1, 2), (1, 4), (2, 4), (2, 5),
                (5, 6), (5, 7), (6, 7), (7, 8), (8, 9),
            ],
        );
        let m0: Vec<Vec<f64>> = vec![
            vec![1., 0., 0.], vec![0., 0., 0.], vec![0., 1., 0.], vec![0., 0., 4.],
            vec![1., 1., 0.], vec![0., 0., 0.], vec![1., 0., 0.], vec![0., 0., 1.],
            vec![0., 0., 1.], vec![0., 0., 0.],
        ];
        let themes = vec!["brown".into(), "black".into(), "yellow".into()];
        let cfg = ThemeConfig {
            self_loop: true,
            weighted: false,
            max_iter: 30,
            tol: -1.0, // never stop early -> exactly 30 steps, like the R reference
            annotated_only: true,
        };
        let tp = theme_proximity(&g, &m0, &themes, &cfg);
        // R reference proximities.
        assert!((tp.proximity[0][1] - 0.9316320).abs() < 1e-5, "brown-black {}", tp.proximity[0][1]);
        assert!((tp.proximity[0][2] - 0.2050967).abs() < 1e-5, "brown-yellow {}", tp.proximity[0][2]);
        assert!((tp.proximity[1][2] - 0.1518247).abs() < 1e-5, "black-yellow {}", tp.proximity[1][2]);
        // symmetry + zero diagonal
        assert_eq!(tp.proximity[0][0], 0.0);
        assert!((tp.proximity[1][0] - tp.proximity[0][1]).abs() < 1e-12);
    }
}
