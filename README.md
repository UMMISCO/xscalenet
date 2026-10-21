# xscalenet

An efficient **Rust** reimplementation of the [ScaleNet](scalenet/) R package —
network reconstruction from high-dimensional data via spectral decomposition,
plus the **Spectral Consensus Strategy (SCS)** that fuses several embedded
reconstruction methods into one robust network.

This is the **core path** of the original package: spectral variable subsets →
ARACNE and/or Bayesian hill-climbing local reconstructions → gather → consensus,
with optional PAM discretization. It keeps everything in memory and reconstructs
subgraphs in parallel, avoiding the R version's per-subgraph disk round-trips.

📖 **Documentation:** see [`docs/`](docs/) — the
[algorithm](docs/algorithm.md), [architecture](docs/architecture.md),
[R mapping](docs/r-mapping.md), and [performance](docs/performance.md).

## Algorithm (and where it maps to the R code)

| Stage | R source | Rust module |
|-------|----------|-------------|
| Load + factor-encode input | `setEnvironment.R` | `data.rs` |
| (optional) PAM discretization | `discretize.R` | `discretize.rs` |
| Affinity = mutual-information matrix (Miller–Madow, `mi.mm`) | `LaplacianRW.R` / `minet::build.mim` | `mi.rs` |
| Random-walk Laplacian `Lrw = I − D⁻¹W`, drop zero-degree nodes | `LaplacianRW.R` | `spectral.rs` |
| Symmetric eigendecomposition + two-line "elbow" for `k` | `computeEigenVectVal.R`, `computeEigenVectBestK.R` | `spectral.rs` |
| Spectral subsets: top-`m` / bottom-`m` per eigenvector | `createSubInputData.spectral.R` | `scalenet.rs::prepare` |
| ARACNE (DPI) local reconstruction | `rMethod.aracne.R` / `minet::aracne` | `aracne.rs` |
| Bayesian hill-climbing (BDeu score) | `rMethod.hc.R` / `bnlearn::hc` | `hc.rs` |
| Gather subgraphs → presence / orientation / frequency | `scaleNet.gatherSubGraphs.R` | `scalenet.rs::run_method` |
| Threshold + Spearman sign → edge list | `scaleNet.convertToScaleNetFormat.R` | `scalenet.rs`, `output.rs` |
| SCS rank-fusion consensus | `scs.R` | `scs.rs` |

## Build

```sh
cargo build --release
```

## Usage

```sh
xscalenet <input(.tsv|.gz)> <output_dir> [options]
```

`input`: tab-separated, **features in columns, observations in rows**, with a
header of feature names. Plain or gzip (`.gz`) is auto-detected. Use `--transpose`
for a features-in-rows matrix (feature id in column 1, header = observation ids).

| Option | Default | Meaning |
|--------|---------|---------|
| `--methods` | `aracne,bayes_hc` | reconstruction methods to embed |
| `--var-perc` | `0.2` | fraction of variables per subgraph (`m`) |
| `--eigen-perc` | `auto` | fraction of eigenvectors (`auto` = elbow heuristic) |
| `--thresh` | `0.3,0.8` | presence-frequency thresholds |
| `--epsilon` | `0.001` | ARACNE DPI tolerance |
| `--restart` | `21` | hill-climbing random restarts |
| `--seed` | `6196` | PRNG seed |
| `--discretize` | off | PAM-discretize numeric columns first |
| `--transpose` | off | input has features in rows |
| `--edges-only` | off | write only present edges (skip absent pairs) |
| `--max-clusters` | `5` | max PAM clusters per discretized column |
| `--maxp` | unlimited | max parents per node (hill-climbing) |
| `--max-iter` | `10000` | hill-climbing iteration cap |
| `--iss` | `1.0` | BDeu imaginary sample size |
| `--metrics` | off | per-node topological metrics on the consensus network |
| `--metrics-weighted` | off | use `|ecorr|` edge weights for metrics/ATC |
| `--atc <annot.tsv>` | — | propagate a node annotation (ATC centrality) |
| `--atc-col <name>` | col 2 | annotation column to propagate |
| `--theme-network` | off | build the multi-layer theme-proximity network |
| `--json` | off | also export `graph.json` (nodes+metrics+edges+theme layer) |
| `--scale-free` | off | sweep `|ecorr|` thresholds, fit `P(k) ~ k^-γ`, recommend one |
| `--iterative-eigen` | off | Lanczos eigensolver (needs a fixed `--eigen-perc`; not bit-identical) |
| `--legacy-lrw` | off | reproduce the pre-0.2 **column-order-dependent** spectral step (see below) |
| `--hc-orders` | `400` | random variable orderings fused per subgraph by the Bayesian step |
| `--hc-order-freq` | `0.0` | drop arcs seen in fewer than this fraction of orderings |

Example (matches the `scs` call in the R package):

```sh
xscalenet ../scalenet/inst/extdata/pop2mat.txt out --discretize
```

Or run directly on a raw gzipped, features-in-rows matrix:

```sh
xscalenet all_samples.20M.nsat.gz out --transpose --discretize
```

### Reproducibility: column order

Up to version 0.1 the reconstruction **depended on the order of the input
columns**: renaming variables changed nothing, but shuffling them changed the
network (edge Jaccard 0.23–0.49 across 10 permutations on four real datasets).

`Lrw = I − D⁻¹W` is not symmetric, yet R ScaleNet decomposes it with
`eigen(symmetric = TRUE)`, which reads only the **lower triangle** — so each
off-diagonal entry is divided by the degree of whichever endpoint has the larger
index. Permuting the columns changes *the matrix being decomposed*, hence the
eigenvectors, the spectral subsets, and every edge. xscalenet reproduced this
faithfully, because it was porting it.

Since 0.2 the default path decomposes the symmetric normalised Laplacian
`Lsym = I − D^-1/2 W D^-1/2` and maps its eigenvectors back to `Lrw`'s
(`v = D^-1/2 u`, same eigenvalues), canonicalises eigenvector signs, and orders
each subgraph's variables by name before hill climbing. The result is **exactly
order-invariant** (Jaccard 1.0000 over 10 permutations on all four datasets).
Two side effects are expected and correct: the elbow heuristic now sees the real
eigenvalues, so `k` and the edge count both rise; and results differ from
pre-0.2 runs. Pass `--legacy-lrw` to reproduce those.

Two inputs are *intrinsically* ambiguous, and the run now warns instead of
pretending otherwise:

* **duplicate or MI-equivalent variables** make eigenvalues exactly degenerate,
  so no unique eigenvector exists;
* a **`G+`/`G-` cut inside a block of numerically-tied components** makes subset
  membership arbitrary at the available precision.

Collapse duplicate variables before reconstruction to remove both.

### Reproducibility: variable names

Order-invariance is not enough. Hill climbing keeps the *first* best-scoring move
it finds, so the order in which it visits candidate moves decides every tie and
steers the greedy path. Version 0.2 fixed that order to the node names: that made
runs **reproducible but not invariant** — the network was still a function of an
arbitrary ordering, now alphabetical, so a taxonomic rename
(`Chromis` → `Pycnochromis`) moved edges (Jaccard 0.65–0.93 on the reference
matrices). ARACNE was unaffected throughout, since its DPI decisions read the
unmodified MI matrix.

Since 0.3 the Bayesian step runs under `--hc-orders` **random** variable
orderings per subgraph and fuses them. Each arc votes with the fraction of
orderings that inferred it, and those fractional votes feed the existing
presence-frequency consensus — so `--thresh` keeps its published meaning, "the
share of local reconstructions that found this edge", now estimated from a soft
vote instead of one tie-break. `--hc-order-freq` is only a noise floor below
which an arc contributes nothing. The orderings are drawn from `--seed` and the
subgraph index, neither of which depends on names or column positions, so the
result is invariant to renaming and reordering **up to Monte-Carlo error**, which
shrinks as the ordering count grows:

| `--hc-orders` | rename Jaccard | s/run (135 vars) |
|---|---|---|
| 1 (= 0.2 behaviour) | 0.65–0.93 | 0.1 |
| 24 | 0.91–0.96 | 2.0 |
| 100 | 0.95–0.98 | 4.7 |
| 200 | 0.98–0.99 | 8.9 |
| **400 (default)** | **0.99** | **18** |

The default is chosen to clear a 2 % edge difference. Lower it knowingly on large
matrices, or use `--methods aracne`, which is exactly invariant at no cost.

**`--restart` defaults to 0 when ordering consensus is on.** Restarts and
orderings both diversify the greedy search, and under consensus the restarts are
redundant: keeping 21 of them alongside 400 orderings *lowered* rename
reproducibility (0.974–0.978 vs 0.988–0.995) while doubling the runtime, because
each restart adds variation that is never averaged. An explicit `--restart` is
always honoured, and `--hc-orders 1` keeps the published default of 21.

Hill climbing remains **seed**-dependent, by design — it is a stochastic search.
Record `--seed`. At the default ordering count the seed and rename sensitivities
coincide (both ≈0.99), which is the expected floor: renaming resamples the
orderings, so it is statistically a seed change.

### Outputs (same layout as the R package)

```
out/globalGraph/globalNet_presFreq<t>/edgesList.<method>.txt   # per method, per threshold
out/consensusGraph/globalNet_presFreq<t>/edgesList.txt          # fused consensus network
out/consensusGraph/globalNet_presFreq<t>/edgesList.consensus.rawAvg.txt
```

Edge-list columns: `x  y  epresenceScore  epresence  eorientScore  eorient  ecorr`.

With `--metrics` / `--atc`, a per-node table is also written:

```
out/consensusGraph/globalNet_presFreq<t>/nodeAnalysis.txt
# node  degree  strength  betweenness  closeness  clustering  community  atc
```

`community` is a Louvain module id; `atc`/`atc_last` are the Annotation
Transcriptional Centrality (Prifti et al. 2010) of the propagated annotation.
With an annotation, per-module **Fisher enrichment** is also written:

```
out/consensusGraph/globalNet_presFreq<t>/enrichment.txt      # per-module Fisher enrichment
out/consensusGraph/globalNet_presFreq<t>/themeNetwork.txt    # theme-proximity network (--theme-network)
out/consensusGraph/globalNet_presFreq<t>/graph.json          # self-contained graph for web/IntegraNet (--json)
```

See [docs/integration.md](docs/integration.md) for the library API and JSON schema.

Example:

```sh
xscalenet net.tsv out --discretize --metrics --atc taxonomy.tsv --atc-col genus
```

## Validation

Run against the bundled R fixtures (`scalenet/tests/scalent_results/`, generated
by the original package on `pop2mat.txt`):

| Output | Result vs R fixture |
|--------|---------------------|
| **ARACNE** edges @ 0.3 / 0.8 | **exact match** — 205/205 and 168/168, zero differences, *and* `epresenceScore`/`ecorr` agree to 15 significant digits |
| Eigenvector count / subgraph structure | exact (`k=13`, `m=9`, 24 subgraphs) |
| Consensus format, ranking, `ecorr` | matches (top edges same nodes, same order, identical `ecorr`) |
| **bayes_hc** edges | ~80 % edge overlap (121/153 @ 0.3) |

The pipeline surrounding reconstruction (MI, Laplacian, eigen, gather, convert,
Spearman sign, consensus math) is **numerically exact** against R. The bayes_hc
differences are the expected consequence of *algorithmically-faithful* (not
bit-identical) structure learning: our hill-climber differs from `bnlearn` in
move tie-breaking, restart RNG, and search order, and BN structure search is
sensitive to all three. Consensus disagreements trace entirely to bayes_hc —
ARACNE contributes none.

## Performance

The whole `pop2mat` run (both methods, both thresholds, discretization) finishes
in well under a second. The two structurally expensive stages — the O(p²) MI
matrix and the 2·(k−1) subgraph reconstructions — are parallelized with `rayon`.
The dense symmetric eigendecomposition (via `nalgebra`, pure Rust, no system
LAPACK dependency) dominates for large feature counts.

## Not implemented (out of the chosen "core path" scope)

- Alternative subset strategies: `random`, `spectralKmeans`,
  `spectralFuzzyCmeans{Order,Sample}`, `spectralBipartition`.
- Plotting / PDF outputs (eigenvalue plots, best-k plot, `discretize.plot`).

These are comparison-only variants in the original package; the default
`spectral` path is the one exercised by the package's own README example.
