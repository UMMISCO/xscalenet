//! Spectral Consensus Strategy: rank-fuse the per-method edge lists into one
//! consensus network (`scs`).

use std::collections::{HashMap, HashSet};

use crate::data::Dataset;
use crate::scalenet::{EdgeList, EdgeRec, Method};

/// Per-method metadata used by the consensus.
pub struct MethodInfo {
    pub method: Method,
    /// Does the method provide edge orientation?
    pub ort: bool,
    /// Edge weight used for ranking: `true` = `ecorr`, `false` = `epresenceScore`.
    pub weight_ecorr: bool,
}

/// A consensus edge in the ranked "rawAvg" output.
pub struct RawEdge {
    pub idx: usize,
    pub avg_rank: f64,
    pub avg_ort: f64,
    pub ecorr: f64,
}

/// Build the consensus for one threshold from `(info, edgelist)` pairs, in the
/// method order given (later methods win ties for `ecorr`).
///
/// Returns `(consensus_edge_list, ranked_raw_edges)`.
pub fn consensus(
    data: &Dataset,
    methods: &[(MethodInfo, &EdgeList)],
) -> (EdgeList, Vec<RawEdge>) {
    let p = data.ncol();

    // Per-method rank (1-based) of each inferred edge, and its edge count.
    let mut ranks: Vec<HashMap<usize, usize>> = Vec::with_capacity(methods.len());
    let mut counts: Vec<usize> = Vec::with_capacity(methods.len());

    for (info, el) in methods {
        let mut inferred: Vec<usize> = el.edges.keys().copied().collect();
        // weight for ranking
        let weight = |idx: usize| -> f64 {
            let rec = &el.edges[&idx];
            let w = if info.weight_ecorr {
                rec.ecorr
            } else {
                rec.epresence_score
            };
            w.abs()
        };
        // Order by |weight| descending, NaN last (idx tiebreak keeps determinism).
        inferred.sort_by(|&a, &b| {
            let wa = weight(a);
            let wb = weight(b);
            match (wa.is_nan(), wb.is_nan()) {
                (true, true) => a.cmp(&b),
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                (false, false) => wb.total_cmp(&wa).then(a.cmp(&b)),
            }
        });
        let mut rank = HashMap::with_capacity(inferred.len());
        for (r, &idx) in inferred.iter().enumerate() {
            rank.insert(idx, r + 1);
        }
        counts.push(inferred.len());
        ranks.push(rank);
    }

    // Union of inferred edges, in ascending index order (deterministic).
    let mut union_set: HashSet<usize> = HashSet::new();
    for (_, el) in methods {
        union_set.extend(el.edges.keys().copied());
    }
    let mut union: Vec<usize> = union_set.into_iter().collect();
    union.sort_unstable();

    let mut out = EdgeList::new(p);

    let mut raw: Vec<RawEdge> = Vec::with_capacity(union.len());

    for &idx in &union {
        // Normalised score per method with data; lower is better.
        let mut norm_sum = 0.0;
        let mut norm_cnt = 0.0;
        let mut num = 0.0; // weighted orientation numerator
        let mut den = 0.0; // weighted orientation denominator
        let mut ecorr = f64::NAN;

        for (mi, (info, el)) in methods.iter().enumerate() {
            let nm = counts[mi];
            if nm == 0 {
                continue; // method inferred nothing -> skip (its column is NA)
            }
            let denom = (nm + 1) as f64;
            let this_rank = ranks[mi].get(&idx).copied();
            let norm = match this_rank {
                Some(r) => r as f64 / denom,
                None => 1.0, // (nm+1)/(nm+1)
            };
            norm_sum += norm;
            norm_cnt += 1.0;

            if info.ort {
                if this_rank.is_some() {
                    // Orientation from this method's eorient: 1->0, then /2.
                    let mut o = el.edges[&idx].eorient;
                    if o == 1.0 {
                        o = 0.0;
                    }
                    let ortv = o / 2.0; // 2->1, -2->-1, 0->0
                    let w = 1.0 - norm;
                    num += w * ortv;
                    den += w;
                }
            }

            // ecorr: last method (in order) that inferred this edge wins.
            if let Some(rec) = el.get(idx) {
                if !rec.ecorr.is_nan() {
                    ecorr = rec.ecorr;
                }
            }
        }

        let score_mean = if norm_cnt > 0.0 {
            norm_sum / norm_cnt
        } else {
            f64::NAN
        };
        let ort_wmean = if den != 0.0 { num / den } else { 0.0 };

        let eorient = if ort_wmean > 0.0 {
            1.0
        } else if ort_wmean < 0.0 {
            -1.0
        } else {
            0.0
        };
        out.edges.insert(
            idx,
            EdgeRec {
                eorient,
                epresence_score: score_mean,
                eorient_score: ort_wmean,
                ecorr,
            },
        );

        raw.push(RawEdge {
            idx,
            avg_rank: score_mean,
            avg_ort: ort_wmean,
            ecorr,
        });
    }

    // Ranked raw output: ascending consensus score (best first).
    raw.sort_by(|a, b| a.avg_rank.total_cmp(&b.avg_rank));

    (out, raw)
}
