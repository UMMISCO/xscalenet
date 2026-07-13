# Performance & scaling

This document describes the cost model, the measured behavior, the known
bottlenecks, and the optimization direction. Notation: `p` = number of features,
`n` = number of observations, `m = ⌈var_perc·p⌉` = variables per subgraph,
`k` = retained eigenvectors.

## Cost model per stage

| Stage | Cost | Notes |
|---|---|---|
| Discretization (PAM, optional) | O(k²·n²) per column × p, **serial** | dominant; see below |
| MI affinity `W` | O(p²·n) | parallel over pairs |
| Eigendecomposition | O(p³) | dense, full spectrum, though only `k ≪ p` used |
| Subgraph reconstruction | O((k−1) · cost(m)) | parallel; ARACNE O(m³), hill-climbing O(m²·n·restarts·iters) |
| Gather + `fill_ecorr` | O(p²) + O(#edges · n log n) | serial |
| Output (edge lists) | O(p²) rows written | full-pair, every pair emitted |

Two structural facts dominate: the pipeline is **O(p²) in memory and I/O**
(dense pair space) and **O(p³) in compute** (dense eigendecomposition), while the
input data and the resulting networks are **sparse**.

## Measured behavior

On a dense metagenomics matrix (2417 features × 9817 observations, ~87% missing),
pre-binned to ≤5 levels, ARACNE-only, release build on an 8-core laptop:

| features `p` | total (after buffering fix) | peak RSS |
|---:|---:|---:|
| 250  | ~1.5 s | 190 MB |
| 500  | 2.8 s | 319 MB |
| 1000 | 9.4 s | 770 MB |
| 2417 | 152 s (compute-bound) | 1.75 GB |

The reconstruct/gather/output block scales ≈ O(p²). Adding `bayes_hc` multiplies
cost by ~6× at `p=250` and more as `p` (hence `m`) grows.

## Known bottlenecks

### 1. Discretization does not scale (blocker)
`--discretize` runs 1-D PAM (`pam` in [`discretize.rs`](../src/discretize.rs)) on
the **full per-row** value vector — BUILD + exhaustive SWAP, i.e. O(k²·n²) per
SWAP iteration — for every `k` in `pamk`, for every column, **serially**. With
`n ≈ 10⁴` non-missing values per column this does not complete on realistic data.
The fix is a proper 1-D discretizer (exact 1-D k-means DP, or weighted k-medoids
on unique/gridded values) plus `par_iter` across columns.

### 2. Dense O(p²) pair space
`EdgeList` allocates parallel arrays over all `p(p−1)/2` pairs; gather, `ecorr`,
consensus, and output all sweep that space; and every pair is written to disk
(≈2.9M rows per file at `p=2417`, mostly "absent"). A sparse edge map keyed by
present pairs makes these O(#edges).

### 3. Dense eigendecomposition
Only `k ≈ 66` of 2417 eigenpairs are used, yet the full dense spectrum is
computed (O(p³)). An iterative solver (Lanczos / implicitly-restarted) on the
sparsified Laplacian returns just the `k` smallest pairs in ~O(k·nnz).

## Fixes applied

On the reference matrix (2417 features × 9817 obs), ARACNE, the full run went from
effectively unusable to **25 s**:

| Change | Effect |
|---|---|
| **Scalable discretizer** — weighted-representative fast path for large columns + parallel columns | `--discretize` **∞ → 8.3 s** for the whole matrix; small inputs stay bit-identical |
| **Buffered edge-list writers** (`BufWriter`) — was one `write()` syscall per row | removed the I/O stall; **>266 s → 152 s** |
| **ecorr precompute + parallel** — rank each column once, then O(n) Pearson on ranks over present edges | **152 s → 53 s**, output bit-identical |
| **ARACNE reuses the global MI matrix** — slice instead of rebuild per subgraph | reconstruct block **36 s → 8.5 s**; **53 s → 25 s**, bit-identical |
| **Sparse `EdgeList`** — present edges only (was dense O(p²) arrays) | peak RSS **~2.2 GB → ~1.9 GB**, bit-identical; foundation for present-only output |

After these, the **MI + Laplacian + eigendecomposition stage (~14.5 s)** is the
dominant *time* cost (O(p³) eigen + O(p²·n) MI).

### Memory profile (p=2417)

Peak RSS ~1.9 GB is **not** in the edge lists. It is:
1. **~1.5 GB during input parsing** — all p·n ≈ 24 M cells are held as `String`s.
2. **~190 MB** for the precomputed rank matrix (p×n f64), plus reconstruct transients.

## Remaining optimization direction

1. **Compact input parsing** — parse per column into typed values (numeric →
   `f64`, else interned categoricals) instead of 24 M `String`s (~8× less memory).
   This is the current #1 memory target.
2. **Iterative eigensolver** for the `k` smallest eigenpairs (removes the O(p³)
   eigen). Note: the auto-`k` elbow heuristic needs the full spectrum, so this
   requires either a fixed `k` or a separate cheap spectrum estimate.
3. **Present-only output** (optional) — now trivial on the sparse `EdgeList`; cuts
   output file size ~6×.
4. Exploit the ~87% sparsity in the MI pass; parallelize the remaining serial tails;
   bound hill-climbing with a max-parents cap; fix bayes_hc run-to-run nondeterminism.
