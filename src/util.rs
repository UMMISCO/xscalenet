//! Small helpers shared across the pipeline.

/// Number of unordered pairs among `p` items.
pub fn npairs(p: usize) -> usize {
    p * (p.saturating_sub(1)) / 2
}

/// Index of the ordered pair `(a, b)` with `a < b` in `combn` order over `p`
/// items (row-major over the upper triangle), matching R's `combn(names, 2)`.
pub fn pair_index(a: usize, b: usize, p: usize) -> usize {
    debug_assert!(a < b && b < p);
    let base = a * (p - 1) - a * (a.wrapping_sub(1)) / 2;
    base + (b - a - 1)
}

/// Inverse of [`pair_index`]: the ordered pair `(a, b)` with `a < b` at linear
/// `combn` index `idx` over `p` items. O(1).
pub fn pair_at(idx: usize, p: usize) -> (usize, usize) {
    // Row `a` starts at cumulative count s(a) = a*(p-1) - a*(a-1)/2.
    let s = |a: usize| a * (p - 1) - a * (a.wrapping_sub(1)) / 2;
    let pf = p as f64 - 0.5;
    let mut a = (pf - (pf * pf - 2.0 * idx as f64).max(0.0).sqrt()).floor() as usize;
    // Correct any floating-point drift.
    while a > 0 && s(a) > idx {
        a -= 1;
    }
    while s(a + 1) <= idx {
        a += 1;
    }
    let b = idx - s(a) + a + 1;
    (a, b)
}

/// Format a value the way the R outputs do: `NA` for missing, otherwise a
/// full-precision round-trippable number (integers without a decimal point).
pub fn fmt(x: f64) -> String {
    if x.is_nan() {
        "NA".to_string()
    } else if x == x.trunc() && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        format!("{}", x)
    }
}
