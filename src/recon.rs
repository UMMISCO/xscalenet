//! Shared types for local network reconstruction.

/// An inferred edge from a subgraph reconstruction.
///
/// `a` and `b` are global column indices with `a < b`. `eorient` follows the
/// ScaleNet convention on the ordered pair `(a, b)`:
///   * `2`  forward  (a -> b)
///   * `-2` backward (b -> a)
///   * `1`  present but unoriented
///
/// `score` is the edge weight (ARACNE MI weight; `NaN` for hill-climbing).
///
/// `weight` is how strongly this subgraph votes for the edge when the results
/// are gathered: `1.0` for a method that simply infers an edge or not (ARACNE),
/// and the fraction of variable orderings that inferred it for hill-climbing
/// ordering consensus. The gather step sums these votes, so the
/// presence frequency keeps its meaning — the share of local reconstructions
/// that found the edge — while a borderline arc contributes partially instead of
/// flipping a whole vote on a tie-break.
#[derive(Clone, Copy, Debug)]
pub struct SubEdge {
    pub a: usize,
    pub b: usize,
    pub eorient: i32,
    pub score: f64,
    pub weight: f64,
}
