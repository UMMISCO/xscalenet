//! Random-walk Laplacian and its spectral decomposition.
//!
//! Follows `LaplacianRW` + `computeEigenVectVal` + `computeEigenVectBestK`:
//! build the affinity (MIM), drop zero-degree nodes, form the random-walk
//! Laplacian `Lrw = I - D^-1 W`, decompose it, sort eigenvalues ascending, and
//! pick the number of eigenvectors either from a fixed percentage or the
//! two-line "elbow".
//!
//! # Why `Lrw` is not decomposed directly
//!
//! `Lrw` is **not symmetric**. R ScaleNet nonetheless calls
//! `eigen(decompose.mat, symmetric = TRUE)`, which silently reads only the
//! **lower triangle** — so entry `(i, j)` is taken as `-W[i,j]/D[i]` with `i`
//! the *larger* index. Permuting the input columns then changes which degree
//! divides each off-diagonal entry, i.e. it changes *the matrix being
//! decomposed*, not merely its row/column order. The eigenvectors change, the
//! spectral subsets in [`crate::scalenet::prepare`] change, and every
//! reconstructed edge changes with them.
//!
//! The fix: decompose the **symmetric** normalised Laplacian
//! `Lsym = I - D^-1/2 W D^-1/2`, which is genuinely symmetric and permutation-
//! equivariant, then map its eigenvectors back to `Lrw`'s. Since
//! `Lrw = D^-1/2 Lsym D^1/2`, the two share eigenvalues and
//! `v_Lrw = D^-1/2 u_Lsym`. Eigenvector signs are then canonicalised so the
//! `G+` / `G-` split of each eigenvector is order-invariant too.
//!
//! The old, order-dependent behaviour is kept behind `Params::legacy_lrw`
//! (`--legacy-lrw`) so results published before the fix stay reproducible.

use nalgebra::DMatrix;

use crate::mi;

/// Convergence tolerance for the dense symmetric eigensolver on the fixed path.
/// See the note in [`compute`]: smaller is *not* better — at machine epsilon the
/// QL sweep exhausts `EIGEN_MAX_ITER` and returns a less accurate result.
const EIGEN_EPS: f64 = 1e-12;
const EIGEN_MAX_ITER: usize = 10_000;

pub struct Spectral {
    /// Names of the nodes that survived the zero-degree filter, row order of
    /// `eigenvectors`.
    pub node_names: Vec<String>,
    /// Eigenvectors, truncated to `k` columns, sorted by ascending eigenvalue.
    /// `eigenvectors[i][c]` is component `i` of eigenvector `c`.
    pub eigenvectors: Vec<Vec<f64>>,
    /// Number of retained eigenvectors (`k`).
    pub k: usize,
    /// All eigenvalues (ascending), before truncation.
    pub eigenvalues: Vec<f64>,
    /// Number of consecutive pairs among the **retained** eigenvalues that are
    /// numerically equal. A degenerate eigenvalue has no unique eigenvector —
    /// any orthonormal basis of its eigenspace is valid and the solver's choice
    /// depends on the row ordering — so the `G+`/`G-` subsets built from it are
    /// not well defined and stay column-order dependent even on the fixed path.
    /// Caused by duplicate or mutual-information-equivalent input variables.
    pub degenerate_pairs: usize,
    /// The global mutual-information affinity matrix (`p × p`, global column
    /// indices, before the zero-degree filter). Retained so ARACNE can slice it
    /// instead of recomputing per subgraph — `MI(i, j)` is over all rows, so the
    /// slice equals a freshly built subset MIM.
    pub mim: Vec<Vec<f64>>,
}

/// Compute the spectral embedding.
///
/// * `cols` / `levels` — level-coded columns and cardinalities for every node.
/// * `names` — node names aligned with `cols`.
/// * `ncol_all` — original number of columns (used to turn `k%` into a count).
/// * `eigen_perc` — fraction of eigenvectors to keep, or `None` for the elbow.
/// * `legacy_lrw` — reproduce the pre-fix, column-order-dependent behaviour
///   (decompose the lower triangle of the non-symmetric `Lrw`, no sign
///   canonicalisation). See the module docs.
#[allow(clippy::too_many_arguments)]
pub fn compute(
    cols: &[&[u32]],
    levels: &[u32],
    names: &[String],
    ncol_all: usize,
    eigen_perc: Option<f64>,
    subset_m: usize,
    iterative: bool,
    seed: u64,
    legacy_lrw: bool,
) -> Result<Spectral, String> {
    // Affinity = mutual information matrix (>=0, zero diagonal).
    let w = mi::build_mim(cols, levels);
    let p = w.len();

    // Degrees = column/row sums.
    let degree: Vec<f64> = (0..p).map(|i| w[i].iter().sum::<f64>()).collect();

    // Identify and drop zero-degree nodes.
    let keep: Vec<usize> = (0..p).filter(|&i| degree[i] > 0.0).collect();
    if keep.is_empty() {
        return Err("all variables have zero degree".into());
    }
    if keep.len() < subset_m {
        return Err("number of connected variables is smaller than subgraph size".into());
    }
    let q = keep.len();
    let node_names: Vec<String> = keep.iter().map(|&i| names[i].clone()).collect();

    // Build the symmetric matrix to decompose. Diagonal is 1 in both cases
    // (W has a zero diagonal).
    let mut m = DMatrix::<f64>::zeros(q, q);
    for a in 0..q {
        let ia = keep[a];
        m[(a, a)] = 1.0;
        for b in 0..a {
            let ib = keep[b];
            let val = if legacy_lrw {
                // Lower triangle of the non-symmetric Lrw = I - D^-1 W, mirrored
                // — what R's `eigen(symmetric = TRUE)` reads. The divisor is the
                // degree of whichever node has the larger index, which is why
                // this is column-order dependent.
                -w[ia][ib] / degree[ia]
            } else {
                // Lsym = I - D^-1/2 W D^-1/2 — symmetric by construction, so
                // permuting the input permutes the matrix and nothing else.
                -w[ia][ib] / (degree[ia] * degree[ib]).sqrt()
            };
            m[(a, b)] = val;
            m[(b, a)] = val;
        }
    }

    // Iterative fast path: only usable when k is fixed (`eigen_perc`), since the
    // elbow needs the full spectrum. Returns `None` if it did not converge, in
    // which case we fall back to the exact dense solver below.
    let iter_result = if iterative && eigen_perc.is_some() {
        let k = ((eigen_perc.unwrap() * ncol_all as f64).ceil() as usize).clamp(2, q);
        crate::lanczos::smallest_eigenpairs(&m, k, seed).map(|(vals, vecs)| (vals, vecs, k))
    } else {
        None
    };

    let (eigenvalues, mut eigenvectors, k) = if let Some(r) = iter_result {
        r
    } else {
        // Full dense symmetric eigendecomposition.
        //
        // nalgebra's default `symmetric_eigen` stops early: on real inputs here
        // it returns eigenvectors with residuals ||Av - \u{3bb}v|| up to ~6e-5, which is
        // the same order as the top-m selection margin (~6e-5 on BRUVS_own), so
        // membership of the G+/G- subsets flips on solver noise alone. Asking
        // for a tighter tolerance brings the worst residual down to ~3e-10.
        // `EIGEN_EPS` is deliberately not machine epsilon: with eps that small
        // the QL sweep hits its iteration cap and returns a *worse* result than
        // the default. Falls back to the default solver if it cannot converge.
        let eig = if legacy_lrw {
            m.clone().symmetric_eigen()
        } else {
            m.clone()
                .try_symmetric_eigen(EIGEN_EPS, EIGEN_MAX_ITER)
                .unwrap_or_else(|| m.clone().symmetric_eigen())
        };
        let vals = eig.eigenvalues;
        let vecs = eig.eigenvectors;
        let mut order: Vec<usize> = (0..q).collect();
        order.sort_by(|&i, &j| vals[i].total_cmp(&vals[j]));
        let eigenvalues: Vec<f64> = order.iter().map(|&i| vals[i]).collect();
        let perc = match eigen_perc {
            Some(pc) => pc,
            None => best_k_perc(&eigenvalues),
        };
        let mut k = (perc * ncol_all as f64).ceil() as usize;
        k = k.clamp(2, q);
        let mut eigenvectors = vec![vec![0.0f64; k]; q];
        for c in 0..k {
            let col = order[c];
            for r in 0..q {
                eigenvectors[r][c] = vecs[(r, col)];
            }
        }
        (eigenvalues, eigenvectors, k)
    };

    if !legacy_lrw {
        // `m` was Lsym, so the eigenvectors are Lsym's. Map them back to Lrw's:
        // Lrw = D^-1/2 Lsym D^1/2, hence v_Lrw = D^-1/2 u_Lsym for the same
        // eigenvalue. Applies to both the dense and the Lanczos path.
        for (r, row) in eigenvectors.iter_mut().enumerate() {
            let s = 1.0 / degree[keep[r]].sqrt();
            for v in row.iter_mut() {
                *v *= s;
            }
        }
        // An eigenvector is only defined up to sign, and the solver's choice is
        // not permutation-stable — a flip swaps an eigenvector's G+ and G-
        // subsets, which also shifts the per-subgraph hill-climbing seeds. Pin
        // the sign to a permutation-invariant rule.
        for c in 0..k {
            let col: Vec<f64> = (0..q).map(|r| eigenvectors[r][c]).collect();
            if canonical_sign(&col, &node_names) < 0.0 {
                for row in eigenvectors.iter_mut() {
                    row[c] = -row[c];
                }
            }
        }
    }

    let degenerate_pairs = (1..k.min(eigenvalues.len()))
        .filter(|&i| eigenvalues[i] - eigenvalues[i - 1] < 1e-9)
        .count();

    Ok(Spectral {
        node_names,
        eigenvectors,
        k,
        eigenvalues,
        degenerate_pairs,
        mim: w,
    })
}

/// The sign to apply to an eigenvector so the choice does not depend on input
/// column order: the sign of the component sum, falling back to the sign of the
/// largest-magnitude component (ties broken by node name) when the sum is ~0.
/// Both rules are invariant under permutation of the components.
fn canonical_sign(col: &[f64], node_names: &[String]) -> f64 {
    let sum: f64 = col.iter().sum();
    if sum.abs() > 1e-12 {
        return if sum < 0.0 { -1.0 } else { 1.0 };
    }
    let mut best = 0usize;
    for i in 1..col.len() {
        let (a, b) = (col[i].abs(), col[best].abs());
        let takes_over = a > b + 1e-15
            || ((a - b).abs() <= 1e-15 && node_names[i] < node_names[best]);
        if takes_over {
            best = i;
        }
    }
    if col.get(best).copied().unwrap_or(0.0) < 0.0 {
        -1.0
    } else {
        1.0
    }
}

/// Simple OLS fit of `y ~ x`, returning (intercept, slope, r_squared).
fn lm(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    let slope = if sxx == 0.0 { 0.0 } else { sxy / sxx };
    let intercept = my - slope * mx;
    let r2 = if sxx == 0.0 || syy == 0.0 {
        0.0
    } else {
        (sxy * sxy) / (sxx * syy)
    };
    (intercept, slope, r2)
}

/// The two-regression elbow heuristic from `computeEigenVectBestK`.
///
/// Returns the retained fraction `best.k / p` (rounded to 4 digits like R).
fn best_k_perc(eigenvalues: &[f64]) -> f64 {
    let step = 10usize;
    let p = eigenvalues.len();
    let xs: Vec<f64> = (1..=p).map(|i| i as f64).collect();

    // First scan: symmetric windows around the middle, growing by `step`.
    let mut r_squared: Vec<f64> = Vec::new();
    let mut n = step;
    while n <= p {
        // R: n.start = ceil(p/2) - ceil(n/2); n.end = ceil(p/2) + ceil(n/2)
        let half_p = ((p as f64) / 2.0).ceil() as isize;
        let half_n = ((n as f64) / 2.0).ceil() as isize;
        let n_start = half_p - half_n;
        let n_end = half_p + half_n;
        if n_start <= 0 || n_end > p as isize {
            break;
        }
        let s = n_start as usize; // 1-based
        let e = n_end as usize;
        let (_, _, r2) = lm(&xs[s - 1..e], &eigenvalues[s - 1..e]);
        r_squared.push(r2);
        n += step;
    }
    if r_squared.is_empty() {
        // Not enough eigenvalues for the heuristic; fall back to a small default.
        return round4((2.0 / p as f64).max(0.0));
    }

    let max_r2 = r_squared.iter().cloned().fold(f64::MIN, f64::max);
    // n1 = max index (1-based, times step) where within 0.1 of the best R^2.
    let mut best_idx = 0usize;
    for (idx, &r2) in r_squared.iter().enumerate() {
        if round4(max_r2 - r2) < 0.1 {
            best_idx = idx + 1; // 1-based
        }
    }
    let n1 = best_idx * step;

    let half_p = ((p as f64) / 2.0).ceil() as isize;
    let half_n1 = ((n1 as f64) / 2.0).ceil() as isize;
    let n1_start = (half_p - half_n1).max(1) as usize;
    let n1_end = (half_p + half_n1).min(p as isize) as usize;

    let (int1, slope1, _) = lm(&xs[n1_start - 1..n1_end], &eigenvalues[n1_start - 1..n1_end]);
    let (int2, slope2, _) = lm(&xs[0..n1_start], &eigenvalues[0..n1_start]);

    // Two near-parallel fitted lines give no usable intersection; fall back to a
    // small default rather than a division by (nearly) zero.
    let denom = slope1 - slope2;
    if denom.abs() < 1e-12 {
        return round4((2.0 / p as f64).max(0.0));
    }
    let best_k = (int2 - int1) / denom;
    if !best_k.is_finite() || best_k <= 0.0 {
        return round4((2.0 / p as f64).max(0.0));
    }
    round4(best_k / p as f64)
}

fn round4(x: f64) -> f64 {
    (x * 10000.0).round() / 10000.0
}
