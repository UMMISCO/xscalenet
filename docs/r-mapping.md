# Mapping to the original R package

`xscalenet` ports the **default "spectral + SCS" core path** of the
[scalenet](https://github.com/UMMISCO/scalenet) R package. Each Rust module is a
close transliteration of one or more R source files.

## Stage-by-stage mapping

| Stage | R source | Rust |
|---|---|---|
| Load + factor-encode input | `setEnvironment.R` | `data.rs` |
| PAM discretization (optional) | `discretize.R` | `discretize.rs` |
| MI affinity (Miller–Madow `mi.mm`) | `LaplacianRW.R` / `minet::build.mim` | `mi.rs` |
| Random-walk Laplacian `Lrw = I − D⁻¹W` | `LaplacianRW.R` | `spectral.rs` |
| Eigendecomposition + elbow `k` | `computeEigenVectVal.R`, `computeEigenVectBestK.R` | `spectral.rs` |
| Spectral subsets `G⁺,m` / `G⁻,m` | `createSubInputData.spectral.R` | `scalenet.rs::prepare` |
| ARACNE (DPI) | `rMethod.aracne.R` / `minet::aracne` | `aracne.rs` |
| Bayesian hill-climbing (BDeu) | `rMethod.hc.R` / `bnlearn::hc` | `hc.rs` |
| Gather → presence/orientation/frequency | `scaleNet.gatherSubGraphs.R` | `scalenet.rs::run_method` |
| Threshold + Spearman sign → edge list | `scaleNet.convertToScaleNetFormat.R` | `scalenet.rs`, `output.rs` |
| SCS rank-fusion consensus | `scs.R` | `scs.rs` |

## What the Rust port improves

- **No disk round-trips.** The R package writes every subset and every subgraph
  result to disk and re-reads them to aggregate. `xscalenet` keeps everything in
  memory and passes matrices directly.
- **Parallel reconstruction.** The R core `spectral` path reconstructs subgraphs
  sequentially (its `foreach` loop is commented out); `xscalenet` reconstructs
  them with `rayon`.
- **Self-contained.** No `minet`/`bnlearn`/`fpc` dependency stack and no system
  LAPACK — ARACNE, hill-climbing, MI, and the eigensolver are all in-crate.

## Fidelity (validated against the R fixtures on `pop2mat`)

| Output | Result vs R |
|---|---|
| Eigenvector count / subgraph structure | exact (`k=13`, `m=9`, 24 subgraphs) |
| ARACNE present-edge set | **exact** — 205/205 @0.3 and 168/168 @0.8 |
| ARACNE `epresenceScore` | matches to ~13 significant digits |
| ARACNE `ecorr` magnitude | matches to 15 significant digits |
| Consensus format / ranking / top edges | matches (same nodes, same order) |
| bayes_hc present-edge set | ~80% overlap |

Known differences:

- **`ecorr` sign** flips on ~25% of ARACNE edges. The discretization produces the
  *same partition* (hence identical MI and ARACNE structure) but assigns cluster
  labels in a different value-order for some columns; because `ecorr` is a
  rank-based Spearman on those labels, its sign reverses. Tracked for a fix.
- **bayes_hc structure** (~20% of edges) differs because hand-rolled hill-climbing
  diverges from `bnlearn` in restart RNG, move tie-breaking, and search order —
  Bayesian-network structure search is sensitive to all three. This is expected of
  an *algorithmically faithful* (not bit-identical) port.

## Intentionally out of scope (the "core path" only)

- Alternative subset strategies from the R package (kept there "for comparison"):
  `random`, `spectralKmeans`, `spectralFuzzyCmeansOrder`,
  `spectralFuzzyCmeansSample`, `spectralBipartition`.
- The **3off2** reconstruction method (benchmarked in the paper; already dropped
  from the R package's exported methods).
- Non-`pamk` discretizers (`equalwidth`, `equalfreq`).
- Diagnostic plots (eigenvalue / eigenvector / best-k PDFs).

These would only be needed to reproduce the paper's *method-comparison* figures;
the default `spectral` path is the one exercised by the package's own example.
