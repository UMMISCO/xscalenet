//! Annotation enrichment by the one-sided Fisher exact test, following the R
//! `integranet` `testEnrichment` / `enrichmentAnalysis`.
//!
//! For each annotating term, test whether it is **over-represented** in an
//! interest set (e.g. the nodes of a network module) relative to a reference
//! background (e.g. all annotated nodes), via the 2×2 table
//!
//! ```text
//!                 annotated (term)   other
//!   interest       intr_a            intr_t − intr_a
//!   non-interest   ref_a − intr_a    …
//! ```
//!
//! The `greater` p-value is the hypergeometric right tail `P(X ≥ intr_a)` with
//! `N = ref_t`, `K = ref_a`, `n = intr_t`. p-values are corrected with
//! Benjamini–Hochberg (`p.adjust(method = "BH")`).

/// One enriched-term result.
#[derive(Clone, Debug)]
pub struct EnrichmentTerm {
    pub term: String,
    pub interest_count: usize,
    pub reference_count: usize,
    pub pvalue: f64,
    /// Sample odds ratio `(a·d)/(b·c)` (R reports the conditional-MLE estimate;
    /// this is the simpler sample OR).
    pub odds: f64,
    /// Benjamini–Hochberg adjusted p-value.
    pub qvalue: f64,
}

fn lchoose(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    libm::lgamma(n as f64 + 1.0) - libm::lgamma(k as f64 + 1.0) - libm::lgamma((n - k) as f64 + 1.0)
}

/// One-sided (over-representation) Fisher exact p-value: hypergeometric right
/// tail `P(X ≥ intr_a)` for `N = ref_t`, `K = ref_a`, `n = intr_t`.
pub fn fisher_greater(intr_a: usize, ref_a: usize, intr_t: usize, ref_t: usize) -> f64 {
    let (nn, kk, n) = (ref_t as u64, ref_a as u64, intr_t as u64);
    let a = intr_a as u64;
    if kk > nn || n > nn {
        return f64::NAN;
    }
    let denom = lchoose(nn, n);
    let kmax = n.min(kk);
    let mut p = 0.0;
    for x in a..=kmax {
        // need 0 <= n - x <= N - K
        if x > n {
            break;
        }
        let nk = n - x;
        if nk > nn - kk {
            continue;
        }
        p += (lchoose(kk, x) + lchoose(nn - kk, nk) - denom).exp();
    }
    p.clamp(0.0, 1.0)
}

/// Benjamini–Hochberg adjusted p-values (matches R `p.adjust(method = "BH")`).
pub fn bh_adjust(p: &[f64]) -> Vec<f64> {
    let m = p.len();
    if m == 0 {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..m).collect();
    idx.sort_by(|&i, &j| p[i].total_cmp(&p[j]));
    let mut q = vec![0.0f64; m];
    let mut prev = 1.0f64;
    // Walk ranks high → low, enforcing monotonicity.
    for r in (0..m).rev() {
        let i = idx[r];
        let rank = (r + 1) as f64;
        let val = (p[i] * m as f64 / rank).min(prev);
        q[i] = val.min(1.0);
        prev = q[i];
    }
    q
}

/// Enrichment of the distinct terms in `interest` against the `reference`
/// background (which should contain the interest set). Returns terms sorted by
/// ascending p-value.
pub fn enrichment_analysis(interest: &[String], reference: &[String]) -> Vec<EnrichmentTerm> {
    let intr_t = interest.len();
    let ref_t = reference.len();
    let mut terms: Vec<String> = interest.to_vec();
    terms.sort();
    terms.dedup();

    let count = |list: &[String], t: &str| list.iter().filter(|x| x.as_str() == t).count();

    let mut rows: Vec<EnrichmentTerm> = terms
        .into_iter()
        .map(|t| {
            let intr_a = count(interest, &t);
            let ref_a = count(reference, &t);
            let pvalue = fisher_greater(intr_a, ref_a, intr_t, ref_t);
            // 2×2: a=intr_a, b=intr_o, c=other_a, d=other_o
            let b = intr_t.saturating_sub(intr_a) as f64;
            let c = ref_a.saturating_sub(intr_a) as f64;
            let d = (ref_t.saturating_sub(intr_t)) as f64 - c;
            let a = intr_a as f64;
            let odds = if b * c != 0.0 { (a * d) / (b * c) } else { f64::INFINITY };
            EnrichmentTerm {
                term: t,
                interest_count: intr_a,
                reference_count: ref_a,
                pvalue,
                odds,
                qvalue: 0.0,
            }
        })
        .collect();

    let q = bh_adjust(&rows.iter().map(|r| r.pvalue).collect::<Vec<_>>());
    for (r, &qv) in rows.iter_mut().zip(q.iter()) {
        r.qvalue = qv;
    }
    rows.sort_by(|a, b| a.pvalue.total_cmp(&b.pvalue));
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_r_fisher_and_bh() {
        // interest = a a a b c ; reference (background) = 4a 2b 3c 1d
        let interest: Vec<String> = ["a", "a", "a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let reference: Vec<String> = ["a", "a", "a", "a", "b", "b", "c", "c", "c", "d"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let res = enrichment_analysis(&interest, &reference);
        // Reference from R: for each term, fisher.test(..., alternative="greater")$p.value
        // a: P(X>=3 | N=10,K=4,n=5) ; b: X>=1,K=2 ; c: X>=1,K=3
        let expect_p = |term: &str| res.iter().find(|r| r.term == term).unwrap().pvalue;
        // R fisher.test(alternative="greater")$p.value:
        assert!((expect_p("a") - 0.2619048).abs() < 1e-6, "a p={}", expect_p("a"));
        assert!((expect_p("b") - 0.7777778).abs() < 1e-6, "b p={}", expect_p("b"));
        assert!((expect_p("c") - 0.9166667).abs() < 1e-6, "c p={}", expect_p("c"));
        // R p.adjust(c(.2619,.7778,.9167),"BH") = (0.7857143, 0.9166667, 0.9166667)
        let q = |term: &str| res.iter().find(|r| r.term == term).unwrap().qvalue;
        assert!((q("a") - 0.7857143).abs() < 1e-6, "a q={}", q("a"));
        assert!((q("b") - 0.9166667).abs() < 1e-6, "b q={}", q("b"));
        assert!((q("c") - 0.9166667).abs() < 1e-6, "c q={}", q("c"));
        // sorted ascending by p
        assert_eq!(res[0].term, "a");
    }
}
