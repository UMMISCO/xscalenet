//! ARACNE reconstruction (`minet::aracne`) on a subset of variables.
//!
//! Slices the subset out of the precomputed global mutual-information matrix
//! (identical to rebuilding it, since `MI(i, j)` is over all rows), then applies
//! the Data Processing Inequality: within every triangle (i, j, k) the weakest
//! edge is removed when its weight is below `min(w_ik, w_jk) - eps`. Surviving
//! edges are reported as present, unoriented, weighted by their MI.

use crate::recon::SubEdge;

/// Run ARACNE on the columns identified by `sub` (ascending global indices),
/// reading edge weights from the global MI matrix `mim_global`.
pub fn reconstruct(sub: &[usize], mim_global: &[Vec<f64>], eps: f64) -> Vec<SubEdge> {
    let m = sub.len();
    if m < 2 {
        return Vec::new();
    }
    // Slice the m×m subset MIM out of the global matrix (global indices).
    let mut mim = vec![vec![0.0f64; m]; m];
    for i in 0..m {
        let row = &mim_global[sub[i]];
        for j in 0..m {
            mim[i][j] = row[sub[j]];
        }
    }

    // Data Processing Inequality. Decide removals reading the original weights,
    // then apply them (single pass, like minet).
    let mut remove = vec![vec![false; m]; m];
    for i in 0..m {
        for j in (i + 1)..m {
            for k in 0..m {
                if k == i || k == j {
                    continue;
                }
                let a = mim[i][j];
                let b = mim[i][k];
                let c = mim[j][k];
                // Only act when (i,j) is the weakest edge of the triangle.
                if a <= b && a <= c {
                    let other_min = b.min(c);
                    if a < other_min - eps {
                        remove[i][j] = true;
                        remove[j][i] = true;
                    }
                }
            }
        }
    }

    let mut edges = Vec::new();
    for i in 0..m {
        for j in (i + 1)..m {
            if remove[i][j] {
                mim[i][j] = 0.0;
            }
            if mim[i][j] > 0.0 {
                edges.push(SubEdge {
                    a: sub[i],
                    b: sub[j],
                    eorient: 1,
                    score: mim[i][j],
                    weight: 1.0,
                });
            }
        }
    }
    edges
}
