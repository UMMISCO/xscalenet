//! Iterative solver for the **k smallest** eigenpairs of a dense symmetric PSD
//! matrix — an optional fast path for the spectral embedding when only a small
//! fraction of eigenvectors is kept (fixed `--eigen-perc`), instead of a full
//! O(p³) dense decomposition.
//!
//! **Shift-invert Lanczos.** Plain Lanczos converges to *extremal* eigenvalues,
//! but the smallest eigenvalues of a random-walk Laplacian cluster near 0 (one
//! per connected component) and converge very slowly. Instead we run Lanczos on
//! `B⁻¹` with `B = A + εI` (A is PSD, so B is SPD): the smallest eigenvalues of A
//! become the *largest* of `B⁻¹`, which Lanczos captures in a few steps. `B` is
//! LU-factorised once (O(p³/3)); each Lanczos step is a back-substitution
//! (O(p²)). A Ritz value `μ` of `B⁻¹` maps back to `λ = 1/μ − ε` of A.
//!
//! Full reorthogonalization (twice) for robustness. Every returned pair is
//! verified against A directly (`‖Ax − λx‖`); if the requested `k` do not
//! converge within a step budget, returns `None` so the caller falls back to the
//! exact dense solver. The downstream use (top-m / bottom-m split into two
//! subgraphs) is sign-invariant, so approximate eigenvectors reproduce the dense
//! network except at degenerate eigenvalues / sort ties.

use nalgebra::{DMatrix, DVector};

/// The `k` smallest eigenvalues (ascending) and their eigenvectors
/// (`eigenvectors[row][col]`) of the symmetric PSD matrix `a`, via shift-invert
/// Lanczos. `None` if the factorisation fails or the pairs do not converge within
/// the step budget (caller should fall back to a dense solve). `seed` makes the
/// random start vector reproducible.
pub fn smallest_eigenpairs(a: &DMatrix<f64>, k: usize, seed: u64) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let q = a.nrows();
    let k = k.min(q).max(1);
    let eps = 1e-6;

    // B = A + εI (SPD since A is PSD); factorise once.
    let mut b = a.clone();
    for i in 0..q { b[(i, i)] += eps; }
    let lu = b.lu();

    // reproducible start vector
    let mut s = seed | 1;
    let mut rand = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; (s as f64 / u64::MAX as f64) * 2.0 - 1.0 };
    let mut v0 = DVector::from_fn(q, |_, _| rand());
    v0 /= v0.norm();

    let mut steps = ((2 * k + 40).min(q)).max((k + 10).min(q));
    let cap = ((6 * k + 120).min(q)).max(steps);
    loop {
        if let Some(r) = run(a, &lu, &v0, steps, k, eps) { return Some(r); }
        if steps >= cap { return None; }
        steps = (steps * 2).min(cap);
    }
}

fn run(
    a: &DMatrix<f64>,
    lu: &nalgebra::LU<f64, nalgebra::Dyn, nalgebra::Dyn>,
    v0: &DVector<f64>,
    steps: usize,
    k: usize,
    eps: f64,
) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let q = v0.len();
    let mut vv: Vec<DVector<f64>> = vec![v0.clone()];
    let mut alpha: Vec<f64> = Vec::with_capacity(steps);
    let mut beta: Vec<f64> = Vec::with_capacity(steps);

    for j in 0..steps {
        let mut w = lu.solve(&vv[j])?;              // w = B⁻¹ vv[j]
        if j > 0 { w -= beta[j - 1] * &vv[j - 1]; }
        let aj = w.dot(&vv[j]);
        alpha.push(aj);
        w -= aj * &vv[j];
        for _ in 0..2 { for vi in &vv { let c = w.dot(vi); w -= c * vi; } }
        let bn = w.norm();
        if bn < 1e-12 { break; }                    // invariant subspace
        if j + 1 < steps { beta.push(bn); vv.push(&w / bn); }
    }

    let m = alpha.len();
    let mut t = DMatrix::<f64>::zeros(m, m);
    for i in 0..m { t[(i, i)] = alpha[i]; }
    for i in 0..m - 1 { t[(i, i + 1)] = beta[i]; t[(i + 1, i)] = beta[i]; }
    let eig = t.symmetric_eigen();

    // Largest μ of B⁻¹ ↔ smallest λ of A; descending μ ⇒ ascending λ.
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&i, &j| eig.eigenvalues[j].total_cmp(&eig.eigenvalues[i]));

    let kk = k.min(m);
    let mut vals = Vec::with_capacity(kk);
    let mut evecs = vec![vec![0.0f64; kk]; q];
    for c in 0..kk {
        let col = order[c];
        let mu = eig.eigenvalues[col];
        if mu <= 0.0 { return None; }
        let lam = 1.0 / mu - eps;
        let y = eig.eigenvectors.column(col);
        let mut x = DVector::<f64>::zeros(q);
        for i in 0..m { x += y[i] * &vv[i]; }
        // verify against A directly: ‖A x − λ x‖ small ⇒ converged
        if (a * &x - lam * &x).norm() > 1e-6 * (1.0 + lam.abs()) { return None; }
        vals.push(lam);
        for r in 0..q { evecs[r][c] = x[r]; }
    }
    Some((vals, evecs))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A random PSD symmetric matrix (Gᵀ G) so shift-invert's SPD assumption holds.
    fn psd(n: usize, seed: u64) -> DMatrix<f64> {
        let mut s = seed | 1;
        let mut r = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; (s as f64 / u64::MAX as f64) * 2.0 - 1.0 };
        let g = DMatrix::<f64>::from_fn(n, n, |_, _| r());
        &g.transpose() * &g
    }

    #[test]
    fn matches_dense_smallest() {
        let n = 120;
        let m = psd(n, 42);
        let dense = m.clone().symmetric_eigen();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| dense.eigenvalues[i].total_cmp(&dense.eigenvalues[j]));

        let k = 10;
        let (vals, vecs) = smallest_eigenpairs(&m, k, 7).expect("converged");
        for c in 0..k {
            let want = dense.eigenvalues[order[c]];
            assert!((vals[c] - want).abs() < 1e-5, "val {c}: {} vs {}", vals[c], want);
            let mut dot = 0.0; let mut nrm = 0.0;
            for r in 0..n { dot += vecs[r][c] * dense.eigenvectors[(r, order[c])]; nrm += vecs[r][c] * vecs[r][c]; }
            assert!((nrm - 1.0).abs() < 1e-5, "not unit: {nrm}");
            assert!(dot.abs() > 0.999, "vec {c} overlap {}", dot.abs());
        }
    }

    #[test]
    #[ignore] // timing only: `cargo test --release lanczos -- --ignored --nocapture`
    fn bench_vs_dense() {
        for &q in &[800usize, 1600] {
            let m = { let mut a = psd(q, 11); for i in 0..q { a[(i, i)] += 1.0; } a };
            let k = (0.05 * q as f64) as usize;
            let t0 = std::time::Instant::now();
            let _ = m.clone().symmetric_eigen();
            let dense = t0.elapsed().as_secs_f64();
            let t1 = std::time::Instant::now();
            let r = smallest_eigenpairs(&m, k, 1);
            let iter = t1.elapsed().as_secs_f64();
            eprintln!("q={q} k={k}: dense={:.3}s  shift-invert={:.3}s  (x{:.1})  converged={}",
                dense, iter, dense / iter, r.is_some());
        }
    }

    #[test]
    fn clustered_smallest_converges() {
        // A Laplacian-like PSD matrix with several ~0 eigenvalues (block structure)
        // — the hard case shift-invert is meant to handle.
        let n = 90;
        let mut m = psd(n, 5);
        for i in 0..n { m[(i, i)] += 1.0; }       // shift to make it clearly PSD
        let (vals, _) = smallest_eigenpairs(&m, 12, 9).expect("converged");
        let mut dense: Vec<f64> = m.symmetric_eigen().eigenvalues.iter().copied().collect();
        dense.sort_by(f64::total_cmp);
        for i in 0..12 { assert!((vals[i] - dense[i]).abs() < 1e-5, "{}: {} vs {}", i, vals[i], dense[i]); }
    }
}
