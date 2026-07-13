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
| [`lib.rs`](../src/lib.rs) | crate root; re-exports the public API |
| [`pipeline.rs`](../src/pipeline.rs) | `run(data, config, progress) -> PipelineResult` — the library entry point (no I/O) |
| [`network.rs`](../src/network.rs) | typed domain model: `Network`, `Edge`, `Orientation`, `Consensus`, `PipelineResult` |
| [`main.rs`](../src/main.rs) | thin CLI: parse args, load/discretize, call `pipeline::run`, format |
| [`data.rs`](../src/data.rs) | TSV load, factor-encoding to integer levels, Spearman on codes |
| [`discretize.rs`](../src/discretize.rs) | Optional PAM/`pamk` discretization of numeric columns |
| [`mi.rs`](../src/mi.rs) | Miller–Madow mutual-information matrix (`build_mim`) |
| [`spectral.rs`](../src/spectral.rs) | Random-walk Laplacian, eigendecomposition, elbow `k` heuristic |
| [`scalenet.rs`](../src/scalenet.rs) | `prepare` (subgraphs), `run_method` (reconstruct + gather), `fill_ecorr`; internal `EdgeList` |
| [`aracne.rs`](../src/aracne.rs) | ARACNE with Data Processing Inequality |
| [`hc.rs`](../src/hc.rs) | Bayesian hill-climbing, BDeu score, xorshift64\* PRNG |
| [`recon.rs`](../src/recon.rs) | `SubEdge` type + orientation convention |
| [`scs.rs`](../src/scs.rs) | Spectral Consensus Strategy (normalized-rank fusion) |
| [`metrics.rs`](../src/metrics.rs) | topological metrics: degree, betweenness/closeness (Brandes/BFS), clustering, Louvain communities + modularity |
| [`atc.rs`](../src/atc.rs) | Annotation Transcriptional Centrality (Prifti 2010) — annotation propagation over the network |
| [`enrichment.rs`](../src/enrichment.rs) | Fisher over-representation test + Benjamini–Hochberg (per-module annotation enrichment) |
| [`multilayer.rs`](../src/multilayer.rs) | multi-layer theme-proximity network (FunNet functional layer) — annotation propagation + self-tuning Gaussian theme affinity |
| [`analysis.rs`](../src/analysis.rs) | `analyze()` — one typed call running metrics + ATC + enrichment + theme network; the functional-layer library contract |
| [`output.rs`](../src/output.rs) | Edge-list + node-analysis writers (buffered), formatting the typed `Network` |
| [`util.rs`](../src/util.rs) | pair indexing (`combn` order + O(1) inverse), R-compatible number formatting |

## Library API

The pipeline is usable as a library, decoupled from I/O:

```rust
let data = /* build a Dataset from a TSV / discretize */;
let config = RunConfig { methods, params };
let result: PipelineResult = xscalenet::run(&data, &config, |msg| eprintln!("{msg}"))?;
// result.per_method[method][threshold] : Network
// result.consensus[threshold]          : Consensus { network, ranked }
```

The CLI (`main.rs`) is a thin wrapper: it loads/discretizes the input, calls
`run`, and formats the result via `output`. Tests and a future web frontend
consume the same `PipelineResult`.

## Data model

Two layers: a **typed public boundary** and an **index-based internal core**.

- **Typed (public):** `Network { kind, edges: Vec<Edge> }`, where
  `Edge { source, target, orientation: Orientation, presence_score: Option<f64>,
  orient_score, ecorr }`. `Orientation` is an `enum` (`Forward` / `Backward` /
  `Undirected`) and missing values are `Option<f64>` — the magic `2/-2/1` codes
  and `NaN`-as-NA never leak past the boundary. `NetworkKind` records whether a
  network is per-method or consensus (only affects the output orientation code).
- **Internal (compute):** `Dataset { names, cols: Vec<Vec<u32>>, levels, nrow }`
  (column-major level codes) and a **sparse** `EdgeList` — a
  `HashMap<pair_index, EdgeRec>` holding only present edges. Pairs use `combn`
  upper-triangle order with a closed-form `util::pair_index` and its O(1) inverse
  `util::pair_at`. The typed `Network` is built from the `EdgeList` at the pipeline
  boundary.

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
