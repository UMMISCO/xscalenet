# Changelog

## 0.3.0 — unreleased

### Fixed: the Bayesian step depended on variable *names*

0.2 made the reconstruction invariant to column order, partly by scanning
hill-climbing moves in node-name order. That is reproducible but **not
invariant**: the network stayed a function of an arbitrary ordering, now
alphabetical, so renaming variables while keeping the data identical moved edges
(rename Jaccard 0.65–0.93 on the four reference matrices). In a study where a
taxonomic update renames a genus, that path changes results.

The Bayesian step now runs each subgraph under `--hc-orders` **random** variable
orderings and fuses them. Hill climbing keeps the first best-scoring move it
finds, so the visit order decides every tie; averaging over random orders removes
the dependence on any particular one. The orderings are drawn from `--seed` and
the subgraph index — both invariant to names and column positions — so the result
is invariant to renaming and reordering up to Monte-Carlo error.

**How the votes are fused.** Each arc votes with the fraction of orderings that
inferred it, and those fractional votes go into the *existing* presence-frequency
consensus (`SubEdge::weight`, summed by `gather`). The alternative — thresholding
within each subgraph and casting a binary vote — was measured and rejected: it
puts a knife edge at the threshold where the ordering-frequency distribution is
densest, and it was consistently less reproducible. Fractional voting also keeps
`--thresh` meaning exactly what it did ("the share of local reconstructions that
found this edge"), only estimated from a soft vote rather than from one
tie-break, and it degrades to the old binary vote at `--hc-orders 1`. ARACNE
votes 1.0 and so is bit-identical to before, which the golden fixtures confirm.

### Added

* `--hc-orders` (default **400**) and `--hc-order-freq` (default `0.0`, a noise
  floor). 400 clears a 2 % edge difference under renaming and reordering on the
  reference matrices; the README tabulates reproducibility and runtime against
  the ordering count. The wasm build defaults to 24, since it is single-threaded.
* `SubEdge::weight` — how strongly a subgraph votes for an edge. `gather`
  accumulates `f64` votes instead of integer counts.
* Tests: rename-invariance for both methods at the shipped defaults, and a guard
  that single-ordering hill climbing is reproducible but *not* rename-invariant,
  so the distinction the README draws cannot go stale.

### Changed

* **`--restart` now defaults to 0 when `--hc-orders > 1`.** Restarts and
  orderings both diversify the greedy search; under consensus the restarts are
  redundant and actively harmful — 21 restarts alongside 400 orderings *lowered*
  rename reproducibility (0.974–0.978 vs 0.988–0.995) and doubled the runtime,
  because each restart adds variation that is never averaged. An explicit
  `--restart` is honoured, and `--hc-orders 1` or `--legacy-lrw` keeps 21.
* **Results differ from 0.2.** `--legacy-lrw` still reproduces 0.1.x exactly: it
  forces a single ordering and the published restart count.

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
