# The ScaleNet algorithm

This document describes the method as implemented in `xscalenet`, following the
2016 paper. The goal is to make the *why* and the *formulas* legible next to the
code that implements them.

## The problem

Reconstructing a network over `p` variables from `n` observations is hard when
`p` is large: global structure-learning methods (Bayesian networks, information-
theoretic methods) either do not scale or lose accuracy as dimensionality grows,
and no single method is best across data types. In the paper's microbiome case,
hill-climbing alone did not converge in 48 hours on ~2,000 variables.

ScaleNet's idea: **don't learn one big network — learn many small, overlapping
ones on spectrally-chosen variable subsets, then aggregate.** This (a) keeps each
local problem in a regime where the embedded method is fast and accurate, and
(b) lets several methods be fused into a consensus that is more robust than any
one of them.

## Overview

```
        data (n obs × p features, discrete levels)
                    │
                    ▼
   ┌────────────────────────────────────────────┐
   │ SCS-spectral : find path-related subsets    │
   │   W   = mutual-information affinity (p×p)    │
   │   Lrw = I − D⁻¹W  (random-walk Laplacian)    │
   │   eigendecompose Lrw, keep k eigenvectors    │
   │   per eigenvector → two subsets G⁺,m / G⁻,m  │
   └────────────────────────────────────────────┘
                    │  2·(k−1) subgraphs, m ≪ p vars each
                    ▼
   ┌────────────────────────────────────────────┐
   │ SCS-learn : reconstruct each subgraph       │
   │   ARACNE and/or Bayesian hill-climbing      │
   │   (embarrassingly parallel)                 │
   └────────────────────────────────────────────┘
                    │  per-method subgraph edge lists
                    ▼
   ┌────────────────────────────────────────────┐
   │ gather : presence-frequency + orientation   │
   │   freq = #inferred / #co-selected           │
   │   threshold → per-method global network     │
   └────────────────────────────────────────────┘
                    │  one network per method
                    ▼
   ┌────────────────────────────────────────────┐
   │ SCS-consensus : rank-fusion across methods  │
   └────────────────────────────────────────────┘
                    ▼
              consensus network
```

## 1. Affinity: mutual-information matrix `W`

Vertex similarity is the pairwise **mutual information** between variables,
because MI is a general (non-linear) measure well suited to high-dimensional
biological data. Entropy uses the **Miller–Madow** bias-corrected estimator in
nats:

```
H(X) = − Σ_x p(x) ln p(x) + (m − 1) / (2N)
```

where `m` is the number of non-empty cells and `N` the sample size, and

```
MI(X,Y) = H(X) + H(Y) − H(X,Y),   clamped to ≥ 0, zero on the diagonal.
```

`W` is symmetric `p×p` (implemented in [`mi.rs`](../src/mi.rs), upper triangle
computed in parallel).

## 2. Spectral embedding: the random-walk Laplacian

With degree matrix `D = diag(Σ_j W_ij)`, the **random-walk normalized Laplacian**
is

```
Lrw = I − D⁻¹ W.
```

Zero-degree vertices (no MI with anything) are dropped. The first eigenvectors of
`Lrw` minimize the normalized-cut criterion of the relaxed clustering problem, so
their entries encode each vertex's membership strength to path-related
components. `Lrw` is eigendecomposed and eigenvalues sorted ascending
([`spectral.rs`](../src/spectral.rs)).

*Implementation note:* like the R original, the solver treats `Lrw` as symmetric
(reading the lower triangle), which is a modelling choice baked into ScaleNet —
`Lrw` is not literally symmetric. This is reproduced deliberately for fidelity.

### Choosing `k` (number of eigenvectors)

Either a fixed fraction of `p`, or an **elbow heuristic**: sweep growing windows
of the sorted eigenvalue curve, fit a line to each, find where two fitted lines
(the flat head and the rising elbow) intersect, and take that index as `k`.

## 3. Spectral variable subsets `G⁺,m` and `G⁻,m`

For each retained eigenvector `v` (skipping the trivial smallest one), sort its
entries. The `m` variables with the **most positive** entries form `G⁺,m` and the
`m` **most negative** form `G⁻,m`, where `m = ⌈var_perc · p⌉`. This yields
`2·(k−1)` overlapping subgraphs, each over `m ≪ p` variables
([`scalenet.rs::prepare`](../src/scalenet.rs)).

## 4. Local reconstruction (embedded methods)

Each subgraph is reconstructed independently and in parallel. Two methods are
embedded:

- **ARACNE** ([`aracne.rs`](../src/aracne.rs)): builds the subset MI matrix and
  applies the Data Processing Inequality — for every triangle, remove the weakest
  edge if it is weaker than the other two by more than `epsilon`. Edges are
  undirected; the score is the retained MI weight.
- **Bayesian hill-climbing** ([`hc.rs`](../src/hc.rs)): greedy add/delete/reverse
  search maximizing the **BDeu** score (decomposable, `lgamma`-based, imaginary
  sample size 1), from an empty graph, with random restarts. Edges are directed.

## 5. Gather: presence frequency and orientation

Across all subgraphs a pair `(x_i, x_j)` may be *co-selected* several times. Its
reliability is the **presence frequency**

```
freq(i,j) = (# subgraphs where the edge was inferred) / (# subgraphs where both i,j were selected).
```

Orientation is resolved by majority vote (forward / backward / unoriented on ties).
An edge is retained in the global network when `freq ≥ threshold`
([`scalenet.rs::run_method`](../src/scalenet.rs)). A signed **Spearman**
correlation (`ecorr`) is attached to each retained edge.

## 6. Spectral Consensus Strategy (rank fusion)

Given one global network per method, the consensus
([`scs.rs`](../src/scs.rs)) is a **normalized-rank aggregation** (Borda-style):

- Within each method, rank its inferred edges by `|weight|` (ARACNE → presence
  score, hill-climbing → `|ecorr|`); normalize the rank to `(0,1]`. An edge a
  method did not infer gets the worst rank.
- The consensus score of an edge is the mean normalized rank across methods
  (lower = stronger, more agreed). Edges are output sorted by this score.
- Orientation is a rank-weighted vote across the orienting methods (better-ranked
  edges weigh more).

The result is a single network that is typically more precise than any embedded
method alone.

## Parameters (defaults)

| Symbol | CLI | Default | Meaning |
|---|---|---|---|
| methods | `--methods` | `aracne,bayes_hc` | embedded reconstruction methods |
| `var_perc` (→ m) | `--var-perc` | `0.2` | fraction of variables per subgraph |
| `k` | `--eigen-perc` | `auto` | eigenvectors (auto = elbow heuristic) |
| thresholds | `--thresh` | `0.3,0.8` | presence-frequency cut-offs |
| `epsilon` | `--epsilon` | `0.001` | ARACNE DPI tolerance |
| restarts | `--restart` | `21` | hill-climbing random restarts |
| seed | `--seed` | `6196` | PRNG seed |
| discretize | `--discretize` | off | PAM-discretize numeric columns first |

## 7. Downstream analysis (metrics + ATC)

Once a network is reconstructed, two analyses run on it (opt-in via `--metrics` /
`--atc`, in [`metrics.rs`](../src/metrics.rs) and [`atc.rs`](../src/atc.rs)):

- **Topological metrics** — degree, strength (weighted degree), betweenness and
  closeness centrality, local/global clustering, connected components, and
  **Louvain communities** with modularity.
- **Enrichment** ([`enrichment.rs`](../src/enrichment.rs)) — per-module annotation
  over-representation via the one-sided **Fisher exact** test + **Benjamini–Hochberg**
  correction, validated against R `testEnrichment`/`p.adjust` to <1e-6.
- **Modular centrality** — degree and betweenness computed *within* each node's
  community (vs the total network), the modular-vs-total distinction of Prifti 2010.
- **Multi-layer theme-proximity network** ([`multilayer.rs`](../src/multilayer.rs),
  FunNet 2008 / R `dynamicPropagation`) — annotations are propagated (L2-normalized)
  to convergence, then themes are linked by a self-tuning Gaussian kernel over the
  Euclidean distances of their node-profiles: `proximity(a,b) = exp(−d²/(σ_a·σ_b))`.
  This is the second network layer. **Validated against R to <1e-5.**
- **Annotation Transcriptional Centrality (ATC)** — Prifti et al.,
  *Bioinformatics* 2010. Node annotations (e.g. taxonomy) are propagated over the
  network as a non-linear dynamical system: from a column-normalized node×theme
  membership `T⁰`, iterate `Tᵖ = A·Tᵖ⁻¹` (with a self-loop, per-theme column
  normalization) until the **co-inertia** (RV of column-standardized configs) of
  consecutive configurations converges. Two scores are reported: the
  dist-coinertia sum `Σ_p ‖Δrowᵢ‖·(1−co-inertiaₚ)` (normalized), and
  `atc_last = rowSums(final config)`. This is the "functional layer" of the FunNet
  two-layer model (Prifti et al., *Bioinformatics* 2008). The implementation is
  **numerically validated against the R `integranet::dynamicPropagation`**:
  propagation to <1e-5, co-inertia (vs `ade4::coinertia`) and ATC to <1e-4.
