# Track B — Extend xscalenet

**Audience:** contributors to the engine.
**Prerequisite:** Track A modules [A0](track-user.md#a0--setup-15-min)–[A2](track-user.md#a2--your-first-reconstruction-20-min) — you should have run the CLI once and
seen `k=13, 24 subgraphs, m=9` before you read the code that produces it.
**Outcome:** you can add a reconstruction method, ship it to all four build
targets, and prove it did not break the R fidelity guarantees.

Answers: [solutions.md](solutions.md).

---

## B0 — The map (20 min)

~5.7 kLOC of Rust, four dependencies (**nalgebra** for the dense symmetric
eigendecomposition, **rayon** for data-parallelism, **libm** for `lgamma`,
**flate2** for gzip). No LAPACK, no serde in the core, no I/O framework — the TSV
parser and the CLI are hand-rolled.

```
core pipeline      data · discretize · mi · spectral · scalenet · aracne · hc · scs · pipeline
domain boundary    network (typed Network/Edge/Orientation) · recon · util
analysis layer     metrics · atc · enrichment · multilayer · scalefree · analysis
edges of the crate main · bin/web · wasm · graphjson · output · io · lanczos · par
```

Start here, in this order — it is the shortest path to a mental model:

1. [`pipeline.rs`](../../src/pipeline.rs) (93 lines) — the whole algorithm, readable
   in one screen.
2. [`network.rs`](../../src/network.rs) — the typed boundary everything else crosses.
3. [`scalenet.rs`](../../src/scalenet.rs) — `prepare` / `run_method` / `fill_ecorr`,
   where the real work is orchestrated.
4. [`analysis.rs`](../../src/analysis.rs) — one call that runs the whole functional layer.

Full table: [architecture.md](../architecture.md#modules).

> **Exercise B0.1** — Without opening any other file, read `pipeline.rs` and list
> the five stages of [A1](track-user.md#a1--the-idea-in-ten-minutes-15-min) with the
> function that implements each.

---

## B1 — Trace one run through the code (30 min)

```
main::run
 ├─ io::read_input                 gzip auto-detect
 ├─ io::transpose_tsv              --transpose
 ├─ Dataset::parse_tsv_raw         zero-copy &str slices, column-major
 ├─ discretize::discretize_to_dataset | Dataset::from_raw_columns
 │                                 → factor-encoded u32 level codes
 ├─ pipeline::run
 │   ├─ scalenet::prepare
 │   │   ├─ mi::build_mim          MI affinity W                     (rayon, over pairs)
 │   │   └─ spectral::compute      Lrw = I − D⁻¹W, eigendecompose, elbow k
 │   │       (+ subset selection)  → 2·(k−1) subgraphs of m vars
 │   ├─ data::rank_matrix          ranks once, reused by every ecorr
 │   ├─ per method: scalenet::run_method                             (rayon, over subgraphs)
 │   │   ├─ aracne::reconstruct  |  hc::reconstruct
 │   │   └─ gather → presence frequency + orientation vote → threshold
 │   ├─ scalenet::fill_ecorr       Spearman on the precomputed ranks
 │   └─ scs::consensus             normalized-rank fusion, per threshold
 ├─ output::write_edge_list / write_raw_consensus
 └─ analysis::analyze              metrics · ATC · enrichment · theme layer
```

Two invariants hold this together and are worth stating before you change
anything:

- **Determinism.** Both rayon fan-outs are order-independent, and hill-climbing
  seeds are derived per subgraph index (`seed + idx`), so results do not depend on
  thread scheduling. Two identical runs produce byte-identical edge lists —
  verify it yourself, it is [exercise A4.2](track-user.md#a4--the-parameters-that-matter-30-min).
- **The typed boundary.** Inside the pipeline, edges live in an index-keyed
  `EdgeList` with magic orientation codes and `NaN`-as-NA. Those never escape:
  `Network::from_edge_list` converts them to `Orientation` enums and `Option<f64>`
  at the boundary. Anything you add downstream consumes the typed form.

> **Exercise B1.1** — `fill_ecorr` runs *after* thresholding, not before. What
> would the cost be if it ran on all pairs, at `p = 2417`?
>
> **Exercise B1.2** — `data.rank_matrix()` is computed once in `pipeline::run`
> rather than inside `fill_ecorr`. Find the entry in
> [performance.md](../performance.md#fixes-applied) that quantifies why.

---

## B2 — The data model (25 min)

### `Dataset` — column-major level codes

```rust
Dataset { names: Vec<String>, cols: Vec<Vec<u32>>, levels: …, nrow: usize }
```

Column-major because every consumer (MI, ranks, discretization) works one variable
at a time. Values are **factor-encoded integer levels**, not floats: MI needs
discrete cells, and the discretizer orders cluster labels ascending by medoid
value so the codes stay value-monotonic — which is exactly why xscalenet's `ecorr`
sign is correct where the R fixture's is not
([r-mapping.md](../r-mapping.md#fidelity-validated-against-the-r-fixtures-on-pop2mat)).

### `EdgeList` — sparse pair space

```rust
pub struct EdgeList { pub p: usize, pub edges: HashMap<usize, EdgeRec> }
```

**Presence is implicit: a pair is present iff it has a record.** The key is the
`combn` upper-triangle pair index, with a closed-form `util::pair_index(i, j, p)`
and an O(1) inverse `util::pair_at(idx, p)`. That is what keeps memory O(#edges)
instead of O(p²) while still emitting the full-pair R-compatible output (absent
pairs are written as `NA` on the fly).

### The typed boundary

`Network { kind, edges: Vec<Edge> }` with `Orientation` as an enum and
`Option<f64>` for missing values. `NetworkKind` exists only because the R output
encodes orientation differently for per-method (`2 / −2 / 1`) and consensus
(`1 / −1 / 0`) files — a formatting concern, kept at the edge of the system.

> **Exercise B2.1** — `pair_index` and `pair_at` must round-trip for every pair.
> Write a `#[test]` that checks it for `p = 45` and `p = 1000`.
>
> **Exercise B2.2** — `EdgeRec` still uses `NaN` for NA while `Edge` uses
> `Option<f64>`. Why is that split deliberate rather than sloppy?

---

## B3 — The library API (30 min)

Two entry points, both I/O-free: `run` (reconstruct) and `analyze` (functional
layer). Everything else — CLI, server, wasm, Python bindings — is a wrapper.

Save as `examples/topnodes.rs` and run with
`cargo run --release --example topnodes -- tests/data/pop2mat.txt`:

```rust
use xscalenet::analysis::{analyze, AnalysisConfig};
use xscalenet::data::Dataset;
use xscalenet::pipeline::{run, RunConfig};
use xscalenet::scalenet::{Method, Params};
use xscalenet::{discretize, io};

fn main() -> Result<(), String> {
    let path = std::env::args().nth(1).ok_or("usage: <matrix.tsv>")?;
    let text = io::read_input(&path)?;
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text)?;
    let data = discretize::discretize_to_dataset(names, &raw, nrow, 5);

    let config = RunConfig {
        methods: vec![Method::Aracne],
        params: Params {
            var_perc: 0.2, eigen_perc: None, pres_freq_thresh: vec![0.3],
            epsilon: 0.001, restart: 21, seed: 6196, max_iter: 10_000,
            iss: 1.0, maxp: usize::MAX, iterative_eigen: false, legacy_lrw: false,
            hc_orders: 400, hc_order_freq: 0.0,
        },
    };
    let result = run(&data, &config, |msg| eprintln!("{msg}"))?;

    let net = &result.consensus[0].network;
    let an = analyze(net, result.nodes.len(), None, &AnalysisConfig::default());
    println!("{} edges, {} communities, Q={:.3}",
             net.edges.len(), an.metrics.n_communities, an.metrics.modularity);

    let mut top: Vec<_> = an.metrics.nodes.iter().enumerate().collect();
    top.sort_by(|a, b| b.1.betweenness.partial_cmp(&a.1.betweenness).unwrap());
    for (i, m) in top.iter().take(3) {
        println!("{:<12} degree={:<3} betweenness={:.4}",
                 result.nodes[*i], m.degree, m.betweenness);
    }
    Ok(())
}
```

Expected output (ARACNE only, threshold 0.3):

```
# Kept k=13 eigenvectors, 24 subgraphs, m=9 vars/subgraph
# Reconstructing with aracne...
# Consensus @ presFreq 0.3: 205 edges
205 edges, 3 communities, Q=0.397
GU:1771      degree=20  betweenness=0.1572
GU:124       degree=16  betweenness=0.1368
GU:860       degree=12  betweenness=0.0942
```

Note the consensus of a *single* method is that method's own network — 205 edges,
matching the ARACNE golden fixture. Useful sanity check when you add a method.

> **Exercise B3.1** — Extend it: pass an annotation through `analyze` (see
> `io::read_annotation_labels`) and print the top-3 by ATC instead. Do the two
> rankings agree?
>
> **Exercise B3.2** — `analyze` takes `n_nodes` separately from `network`. Why
> can't it infer the node count from the edges?

---

## B4 — Add a reconstruction method (45 min) — the core exercise

This is the exercise that transfers ownership of the codebase. `Method` is a
two-variant enum, and a third variant has to be threaded through **seven** places.
Finding them all is the point.

The touch points, in dependency order:

| # | File | What to add |
|---|---|---|
| 1 | `src/<yourmethod>.rs` | `reconstruct(sub, …) -> Vec<SubEdge>` — global column indices `a < b`, `eorient` in `{2, −2, 1}`, `score` = weight or `NaN` |
| 2 | `src/lib.rs` | `pub mod <yourmethod>;` |
| 3 | `src/scalenet.rs` | the `Method` enum variant + its `name()` |
| 4 | `src/scalenet.rs::run_method` | the dispatch arm (~line 178) |
| 5 | `src/scs.rs::MethodInfo::for_method` | does it orient? is it ranked by `ecorr` or by presence score? |
| 6 | `src/main.rs`, `src/bin/web.rs`, `src/wasm.rs` | the three independent `"name" => Method::…` parsers |
| 7 | `docs/` + `README.md` | the `--methods` table |

Step 5 is the one people get wrong. `MethodInfo` tells the consensus how to treat
your method:

```rust
Method::Aracne  => MethodInfo { ort: false, weight_ecorr: false }, // unoriented, ranked by presence score
Method::BayesHc => MethodInfo { ort: true,  weight_ecorr: true  }, // oriented,   ranked by |ecorr|
```

Get `weight_ecorr` wrong and your method is ranked on a scale it never produced —
the fusion silently degrades instead of failing.

Step 6 is the design smell to notice: **the same string→enum mapping is duplicated
in three files**, and `wasm.rs` even accepts aliases (`"hc"`, `"bayeshc"`) the
others reject. A `FromStr for Method` in `scalenet.rs` would collapse all three.
That is a legitimate PR.

### Suggested exercise: `correlation`

A deliberately simple method — keep pairs whose `|Spearman|` on the subset exceeds
a cutoff, undirected, score = `|rho|`. It is not a good reconstruction method; it
is a good scaffold, and it lets you see the consensus behave when one method
disagrees loudly with the others.

> **Exercise B4.1** — Implement it, then run
> `--methods aracne,correlation` on `pop2mat`. How does the consensus edge count
> compare with `--methods aracne` alone, and does the ranking in
> `edgesList.consensus.rawAvg.txt` still look sane?
>
> **Exercise B4.2** — Set `weight_ecorr: true` for it, rerun, and diff the rawAvg
> file. Explain what the fusion did.

---

## B5 — Build targets and features (25 min)

One crate, four artefacts, driven by Cargo features:

| Feature | Default | Pulls in | Produces |
|---|---|---|---|
| `parallel` | **on** | rayon | multi-threaded core |
| `graphjson` | off | serde_json | the rich `graph.json` builder |
| `web` | off | axum, tokio, tower-http (+ parallel, graphjson) | `xscalenet-web` server |
| `wasm` | off | wasm-bindgen (+ graphjson) | in-browser reconstruction |

```sh
cargo build --release                                     # CLI + lib
cargo build --release --features web --bin xscalenet-web  # server
./scripts/build-wasm.sh                                   # → web/wasm/
```

The subtlety worth internalizing: **wasm32 is single-threaded**, so the `wasm`
feature turns `parallel` off — and [`src/par.rs`](../../src/par.rs) provides a
sequential shim mirroring the exact slice of the rayon API this crate uses
(`into_par_iter`, `par_iter`, `map_init`). Only *scheduling* changes; the numerics
stay byte-identical, which is why the browser build gives the same network as the
CLI.

`scripts/build-wasm.sh` documents two real toolchain traps in its header (a
Homebrew `cargo` shadowing rustup's, and a nightly `rust-lld` that cannot find
`libLLVM.dylib`). Read them before you debug a link failure from scratch.

> **Exercise B5.1** — Build the wasm module and confirm the in-browser
> reconstruction produces the same edge count as the CLI on `pop2mat`. If it
> doesn't, which of the two invariants in B1 has been broken?
>
> **Exercise B5.2** — Why does `graphjson` exist as a separate feature instead of
> living inside `web`?

---

## B6 — The `graph.json` contract (20 min)

`graph.json` is the suite's interchange format — the CLI, the Python bindings,
the server, the wasm build and the external converters all produce it, and the web
apps consume it. Full spec: [suite.md](../suite.md#graphjson--the-interchange-format).

**There are two producers, and they do not emit the same thing.** Know which you
are looking at:

| | `output::graph_json` (CLI `--json`) | `graphjson::build` (server + wasm) |
|---|---|---|
| `nodes` (metrics, ATC) | ✅ | ✅ |
| `edges` | ✅ | ✅ |
| `theme_network` | ✅ | — |
| `meta` (title, `k`, `subset_m`, `color_by`) | — | ✅ |
| `color_map` / `color_groups` | — | ✅ |
| `themes` / `theme_colors` | — | ✅ |
| `atc_trace` (propagation animation) | — | ✅ |
| `subgraphs` (the fused spectral subsets) | — | ✅ |
| all annotation columns per node | one (`annotation`) | all |

Verify it yourself:

```sh
# minimal (3 top-level keys: nodes, edges, theme_network)
./target/release/xscalenet tests/data/pop2mat.txt out --discretize --metrics \
  --atc annot.tsv --theme-network --json
python3 -c "import json;print(list(json.load(open('out/consensusGraph/globalNet_presFreq0.3/graph.json'))))"

# rich (9 top-level keys)
curl -s -X POST localhost:8788/run -F data=@tests/data/pop2mat.txt \
  -F annotation=@annot.tsv -F color_by=group -F methods=aracne -F thresh=0.3 \
  -F discretize=true | python3 -c "import json,sys;print(sorted(json.load(sys.stdin)))"
```

So a CLI `--json` graph loads in the browser but arrives without colours, without
the propagation timeline, and without the spectral-subgraph view. That asymmetry
is a real gap, not a documented design decision — **converging the two producers
is a good first substantial PR.**

> **Exercise B6.1** — Both producers already depend on the same `analyze` result.
> Sketch the smallest change that lets the CLI emit the rich form. What does it
> need that `analyze` does not currently return?

---

## B7 — Tests and validation against R (30 min)

```sh
cargo test --release      # 12 unit + 4 regression, 1 ignored (a Lanczos benchmark)
```

Three layers, each answering a different question:

1. **Regression vs R goldens** ([`tests/regression.rs`](../../tests/regression.rs)) —
   `spectral_structure_matches_r` (k=13, m=9, 24 subgraphs), `aracne_matches_golden`
   (the present-edge set and `epresenceScore`), `run_api_returns_typed_networks`,
   `bayes_hc_is_deterministic`.
2. **Numerical validation against R inside the modules** — `atc::matches_r_dynamic_propagation`
   (propagation <1e-5, co-inertia and ATC <1e-4 vs `integranet::dynamicPropagation`),
   `multilayer::matches_r_theme_proximity` (<1e-5), `enrichment::matches_r_fisher_and_bh` (<1e-6).
3. **Property/behaviour tests** — `two_triangles_bridge`,
   `bridge_node_has_high_atc`, `hub_network_is_more_scale_free_than_regular`.

The fidelity guarantees these lock in are tabulated in
[r-mapping.md](../r-mapping.md#fidelity-validated-against-the-r-fixtures-on-pop2mat):
ARACNE is exact (205/205 @0.3, 168/168 @0.8), `epresenceScore` matches to ~13
significant digits, and `bayes_hc` overlaps R's `bnlearn` by ~80 % **by design** —
hand-rolled hill-climbing diverges in restart RNG, tie-breaking and search order.
Do not "fix" that to 100 %; do not let it drift below 80 % either.

### The open determinism bug

Reconstruction is deterministic. **Louvain is not.** Four identical runs on
`pop2mat` @0.3:

```
254 edges, 4 communities, modularity 0.343
254 edges, 3 communities, modularity 0.329
254 edges, 3 communities, modularity 0.337
254 edges, 4 communities, modularity 0.343
```

Cause: `louvain` in [`metrics.rs`](../../src/metrics.rs) picks a node's best
community by iterating `w_to: HashMap<usize, f64>`, and coarsens the supergraph
through another `HashMap`. Rust's default hasher is **randomly seeded per
process**, so near-ties (the comparison uses a `1e-12` margin) and the neighbour
visit order in the next pass resolve differently every run.

Fixes, cheapest first: iterate the candidate communities in a deterministic order
(collect into a `Vec` and sort by community id before scanning), and make the
coarsened adjacency deterministic (`BTreeMap`, or sort after collecting). Then add
a regression test asserting a stable `(n_communities, modularity)` on a fixture —
without it the bug comes back silently.

This matters beyond tidiness: `community` feeds `modular_degree`,
`modular_betweenness`, and every enrichment `qvalue` in `enrichment.txt`.

> **Exercise B7.1** — Fix it, and prove the fix with a test that fails on `main`.
>
> **Exercise B7.2** — `bayes_hc_is_deterministic` passes today. Why doesn't it
> catch the Louvain problem?

---

## B8 — Performance (25 min)

Cost model, with `p` features, `n` observations, `m = ⌈var_perc·p⌉`, `k`
eigenvectors:

| Stage | Cost | Parallel? |
|---|---|---|
| PAM discretization | O(k²·n²) per column × p | yes (per column) |
| MI affinity `W` | O(p²·n) | yes (per pair) |
| Eigendecomposition | O(p³) | **no** |
| Subgraph reconstruction | O((k−1)·cost(m)); ARACNE O(m³) | yes (per subgraph) |
| Gather, `fill_ecorr`, consensus, output | O(p²) + O(#edges·n log n) | no |

Measured on 2417 features × 9817 observations, ARACNE, 8-core laptop: the full run
went from effectively unusable to **25 s** through seven bit-identical
optimizations — buffered writers, precomputed ranks, ARACNE reusing the global MI
matrix, sparse `EdgeList`, zero-copy parsing, MI marginal precompute, and a
scalable discretizer. All of them are tabulated in
[performance.md](../performance.md#fixes-applied), and reading that table is the
best available lesson in *how this codebase optimizes*: measure, change one thing,
assert bit-identical output.

What is left is not free:

1. **Eigendecomposition (~10 s, single-threaded, O(p³)).** Only `k ≈ 66` of 2417
   eigenpairs are used. [`lanczos.rs`](../../src/lanczos.rs) exists and is wired to
   `--iterative-eigen`, but it is **off by default** for two reasons: the auto-`k`
   elbow heuristic needs the full spectrum, and a different solver yields different
   (still valid) eigenvectors — so the subgraphs, and hence the network, change.
2. **PAM discretization (~6–8 s).** `CLUSTER_CAP` is a speed/fidelity knob. The
   principled fix (exact 1-D k-means DP) is both faster and more accurate, and
   changes results.
3. **Memory.** Peak RSS at `p=2417` is ~1.9 GB, and **~1.5 GB of it is input
   parsing** — 24 M cells held as strings — not the edge lists. Optimize where the
   bytes actually are.

The rule the codebase follows, and that you should follow: **a change that alters
results is a scientific decision, not a performance decision.** It needs a flag, a
default, and a line in the docs — never a silent speedup.

> **Exercise B8.1** — Time `--iterative-eigen` against the default on `pop2mat`
> with a fixed `--eigen-perc`, and diff the consensus edge lists. Report both
> numbers: how much faster, and how different.
>
> **Exercise B8.2** — At `p = 5000` with `var_perc = 0.2`, which stage dominates,
> and what is the first thing you would measure before optimizing it?

---

## Capstone

Ship a complete change:

1. Implement the `correlation` method from B4 across all seven touch points.
2. Add a `FromStr for Method` and use it in all three parsers (the B4 design smell).
3. Fix the Louvain determinism bug from B7, with a test that fails without the fix.
4. Verify: `cargo test --release` green, `pop2mat` consensus byte-identical to
   `main` for `--methods aracne,bayes_hc`, and the browser still loads a graph
   reconstructed through `POST /run`.
5. Write the PR description: what changed, what stayed bit-identical, what did not
   and why that is acceptable.
