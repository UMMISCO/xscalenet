# Changelog

## 0.2.0 — unreleased

### Fixed: the reconstruction depended on input column order

Shuffling the columns of an input matrix changed the reconstructed network
(edge Jaccard **0.23–0.49** across 10 permutations on four real genus
presence/absence datasets; renaming columns changed nothing, so the networks
were a function of column *position*). It is now exactly invariant
(**Jaccard 1.0000**).

Three independent causes, all fixed:

1. **The decomposed matrix was order-dependent** (`src/spectral.rs`).
   `Lrw = I − D⁻¹W` is not symmetric, but it was decomposed as if it were, from
   its lower triangle — mirroring R ScaleNet's `eigen(decompose.mat,
   symmetric = TRUE)`. That divides entry `(i, j)` by the degree of whichever
   endpoint has the larger index, so a permutation changes the matrix itself.
   The default path now decomposes the symmetric normalised Laplacian
   `Lsym = I − D^-1/2 W D^-1/2` and maps its eigenvectors back to `Lrw`'s via
   `v = D^-1/2 u` (identical eigenvalues, since `Lrw = D^-1/2 Lsym D^1/2`).
   Eigenvector signs are canonicalised, since a sign flip swaps an
   eigenvector's `G+`/`G-` subsets and shifts the per-subgraph hill-climbing
   seeds.

2. **The eigensolver was not accurate enough for the selection it fed**
   (`src/spectral.rs`). nalgebra's default `symmetric_eigen` returned
   eigenvectors with residuals `‖Av − λv‖` up to **6e-5** on these inputs —
   the same order as the top-`m` selection margin (6e-5 on one dataset), so
   subset membership flipped on solver noise alone. The fixed path requests a
   tighter tolerance (`1e-12`, bringing the worst residual to ~3e-10) and falls
   back to the default solver if it cannot converge. Note that asking for
   machine epsilon is *worse*: the QL sweep then exhausts its iteration cap.

3. **Hill climbing scanned variables in input column order** (`src/scalenet.rs`,
   `src/hc.rs`). `climb` keeps the first best move in local index order, and
   `perturb` draws local indices, so both the tie-breaking and the restart
   perturbations followed the input order. Each subgraph's variables are now
   ordered by node name before the search, making it a function of the variable
   *set* alone.

### Added

* `--legacy-lrw` — reproduce the pre-0.2 behaviour exactly (all three points
  above), so analyses published with 0.1.x can be regenerated.
* Warnings for the two cases where the subsets are **intrinsically** ambiguous
  and no tie-break rule can help: degenerate retained eigenvalues, and a
  `G+`/`G-` cut falling inside a block of numerically-tied components. Both are
  caused by duplicate or mutual-information-equivalent input variables (the
  bundled `pop2mat` fixture has 19 duplicate columns and triggers both);
  collapse such variables before reconstruction.
* `Spectral::degenerate_pairs` and `Prep::boundary_ties` expose those counts to
  library callers.
* `tests/order_invariance.rs` — permutation invariance of the spectral subsets
  and of both methods' edge sets, verification that the back-mapped eigenvectors
  really satisfy `Lrw·v = λv` against a direct multiplication by the asymmetric
  `Lrw`, agreement between the dense and Lanczos paths, and a guard that the
  legacy path stays order-dependent.
* `--scale-free` and `--iterative-eigen` are now documented in the README (they
  already existed in the CLI).

### Changed

* **Results differ from 0.1.x.** The elbow heuristic now sees the correct
  eigenvalues, so `k` grows and with it the edge count (on the four reference
  datasets: `k` 5→9, 8→11, 4→4, 6→10; consensus edges at presFreq 0.3
  166→294, 679→846, 111→100, 263→400). Use `--legacy-lrw` for the old numbers.
* The R-fidelity regression tests now pin `legacy_lrw = true`: the golden
  fixtures were produced by R's own order-dependent spectral step, so matching
  them requires reproducing it.

### Note on R ScaleNet

The same flaw is present in the R package and was confirmed there, not just
inferred: `LaplacianRW.R` builds `decompose.mat = diag(n) − D⁻¹ %*% W`
(asymmetric) and `computeEigenVectVal.R` calls `eigen(..., symmetric = TRUE)`;
`eigen(matrix(c(1,99,7,1),2,2,byrow=TRUE), symmetric=TRUE)` returns the
eigenvalues of the lower-triangle symmetrisation, and permuting a toy affinity
matrix through the R code path changes the subsets (per-eigenvector subset
Jaccard 0.00–0.60). Applying the `Lsym` fix in R restores Jaccard 1.00.
Analyses produced with R ScaleNet are affected in the same way.
