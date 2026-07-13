# Architecture

`xscalenet` is a small (~2 kLOC), self-contained Rust crate with four
dependencies: **nalgebra** (dense symmetric eigendecomposition), **rayon**
(data-parallelism), **libm** (`lgamma` for BDeu), and **flate2** (pure-Rust gzip
input). No LAPACK, no serde, no I/O framework — TSV parsing and the CLI are
hand-rolled.

## End-to-end pipeline

```
 main.rs::run
   │
   ├─ data::parse_tsv_raw            read TSV → column-major raw strings
   ├─ discretize::discretize_raw     (optional) PAM-discretize numeric columns
   ├─ data::from_raw_columns         factor-encode columns → integer levels
   │
   ├─ scalenet::prepare              ── shared spectral preprocessing ──
   │     ├─ mi::build_mim                MI affinity W               (rayon)
   │     ├─ spectral::compute            Lrw, eigendecompose, pick k
   │     └─ (subset selection)           2·(k−1) spectral subgraphs
   │
   ├─ for each method:
   │     ├─ scalenet::run_method      reconstruct every subgraph      (rayon)
   │     │     ├─ aracne::reconstruct
   │     │     └─ hc::reconstruct
   │     │   then gather → presence/orientation/frequency → threshold
   │     ├─ scalenet::fill_ecorr      Spearman sign per present edge
   │     └─ output::write_edge_list   per-method global edge lists
   │
   └─ for each threshold:
         ├─ scs::consensus            rank-fusion across methods
         └─ output::write_edge_list / write_raw_consensus
```

## Modules

| Module | Responsibility |
|---|---|
| [`main.rs`](../src/main.rs) | CLI parsing, orchestration, per-method/threshold output wiring |
| [`data.rs`](../src/data.rs) | TSV load, factor-encoding to integer levels, Spearman on codes |
| [`discretize.rs`](../src/discretize.rs) | Optional PAM/`pamk` discretization of numeric columns |
| [`mi.rs`](../src/mi.rs) | Miller–Madow mutual-information matrix (`build_mim`) |
| [`spectral.rs`](../src/spectral.rs) | Random-walk Laplacian, eigendecomposition, elbow `k` heuristic |
| [`scalenet.rs`](../src/scalenet.rs) | `prepare` (subgraphs), `run_method` (reconstruct + gather), `fill_ecorr` |
| [`aracne.rs`](../src/aracne.rs) | ARACNE with Data Processing Inequality |
| [`hc.rs`](../src/hc.rs) | Bayesian hill-climbing, BDeu score, xorshift64\* PRNG |
| [`recon.rs`](../src/recon.rs) | `SubEdge` type + orientation convention |
| [`scs.rs`](../src/scs.rs) | Spectral Consensus Strategy (normalized-rank fusion) |
| [`output.rs`](../src/output.rs) | Edge-list writers (buffered) |
| [`util.rs`](../src/util.rs) | pair indexing (`combn` order), R-compatible number formatting |
| [`lib.rs`](../src/lib.rs) | module declarations |

## Data model

The pipeline is **column-major and index-based**, mirroring the R original:

- `Dataset { names, cols: Vec<Vec<u32>>, levels, nrow }` — each column is a vector
  of integer level codes, ordered by sorted string value (matching R factor
  ordering).
- Variable pairs use a fixed **`combn` upper-triangle order**; a pair maps to a
  linear index via a closed-form formula (`util::pair_index`).
- `EdgeList` stores **parallel arrays over all `p(p−1)/2` pairs**
  (`epresence`, `eorient`, `epresence_score`, `eorient_score`, `ecorr`).
- Orientation is encoded as `2` (forward), `-2` (backward), `1` (unoriented),
  with `NaN` meaning absent/NA.

> This dense pair-space model is the main structural inheritance from R and the
> source of the pipeline's O(p²) memory and time. See
> [performance.md](performance.md) and the roadmap for the sparse replacement.

## Parallelism

Two `rayon` fan-outs, both deterministic:

1. The MI matrix upper triangle (`mi::build_mim`) — parallel over pairs.
2. Subgraph reconstruction (`scalenet::run_method`) — parallel over the
   `2·(k−1)` subgraphs. Hill-climbing seeds are derived per-subgraph-index
   (`seed + idx`) so results are reproducible regardless of scheduling.

Gather, `fill_ecorr`, consensus, and output are currently serial.

## CLI

```
xscalenet <input(.tsv|.gz)> <output_dir> [options]
```

Input is auto-detected as plain or gzip. `--transpose` accepts features-in-rows.

| Option | Default | Meaning |
|---|---|---|
| `--methods` | `aracne,bayes_hc` | reconstruction methods to embed |
| `--var-perc` | `0.2` | fraction of variables per subgraph (`m`) |
| `--eigen-perc` | `auto` | fraction of eigenvectors (`auto` = elbow) |
| `--thresh` | `0.3,0.8` | presence-frequency thresholds |
| `--epsilon` | `0.001` | ARACNE DPI tolerance |
| `--restart` | `21` | hill-climbing random restarts |
| `--seed` | `6196` | PRNG seed |
| `--discretize` | off | PAM-discretize numeric columns first |
| `--transpose` | off | input has features in rows (feature id in col 1) |
| `--edges-only` | off | write only present edges |
| `--max-clusters` | `5` | max PAM clusters per discretized column |
| `--maxp` | unlimited | max parents per node (hill-climbing) |
| `--max-iter` | `10000` | hill-climbing iteration cap |
| `--iss` | `1.0` | BDeu imaginary sample size |

### Output layout (same as the R package)

```
<out>/globalGraph/globalNet_presFreq<t>/edgesList.<method>.txt
<out>/consensusGraph/globalNet_presFreq<t>/edgesList.txt
<out>/consensusGraph/globalNet_presFreq<t>/edgesList.consensus.rawAvg.txt
```

Edge-list columns: `x  y  epresenceScore  epresence  eorientScore  eorient  ecorr`.
