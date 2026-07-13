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

- **Buffered edge-list writers** — the writers wrote each of the millions of rows
  with an unbuffered `writeln!` straight to the `File` (one `write()` syscall per
  row), making large runs I/O-bound. Wrapping them in `BufWriter` (output bytes
  unchanged) gave a **5–6× speedup at `p=500–1000`** and removed the I/O stall at
  `p=2417`, leaving the run compute-bound.

## Optimization direction (summary)

1. Replace the discretizer and parallelize it across columns (unblocks `--discretize`).
2. Sparse edge model + present-only output (removes the O(p²) memory/I/O).
3. Iterative eigensolver for the `k` smallest eigenpairs (removes the O(p³) compute).
4. Exploit the ~87% sparsity in the MI pass; parallelize the serial tails.
5. Bound hill-climbing with a max-parents cap.

See the (internal) roadmap for sequencing and status.
