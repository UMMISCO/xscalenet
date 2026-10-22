# Solutions

Numbers below were produced on `tests/data/pop2mat.txt` with the release build.
Where an exercise depends on your own data or judgement, the answer gives the
reasoning to look for rather than a value.

---

## Track A

**A0.1 — why `file://` fails.** The app fetches `projects.json` and the
`*.graph.json` files at runtime, and browsers block `fetch` on the `file://`
origin. The page also probes `/health` to decide whether the reconstruction
button is live. Serve over HTTP (`xscalenet-web`, or `python3 -m http.server`
for the static-only features).

**A1.1 — `m` and the subgraph count.** `m = ⌈0.2 · 45⌉ = 9`. With `k = 13`
eigenvectors, `2·(k−1) = 24` subgraphs — the trivial smallest eigenvector is
skipped, and each of the remaining 12 contributes a `G⁺,m` and a `G⁻,m`. Both
numbers are printed by the run and asserted by `spectral_structure_matches_r`.

**A2.1 — extra thresholds.** Three folders appear
(`globalNet_presFreq0.3 / 0.5 / 0.8`), and the 0.3 edge list is **byte-identical**
to the two-threshold run. Thresholds are independent cuts applied to the *same*
gathered presence frequencies — adding one costs a comparison sweep, not a
re-reconstruction. Consensus sizes: 254 @0.3, 250 @0.5, 207 @0.8. Passing several
thresholds is therefore nearly free, and you should.

**A2.2 — `--discretize` does nothing here.** `head -3 tests/data/pop2mat.txt`
shows a **binary 0/1** matrix — presence/absence data, already discrete. PAM on a
two-level column returns those two levels, so the factor codes are unchanged and
every downstream stage is identical. On relative abundances it would matter
enormously: without discretization the loader factor-encodes by **string order**,
so MI would be computed over ~`n` singleton levels (no shared cells, meaningless
affinity) and the `ecorr` sign can invert on multi-digit numerics. The stage
responsible is A1 step 1, the MI affinity — it needs discrete cells to estimate
`p(x, y)` at all.

**A3.1 — consensus vs per-method counts.** At threshold 0.8:

```
aracne 168 · bayes_hc 119 · consensus 207 · union 207 · intersection 80
```

The consensus is exactly the **union**, and it must be: rank fusion assigns an
edge a method did not infer the *worst* rank rather than dropping it, so any edge
inferred by at least one method survives into the consensus network. The
consensus's contribution is not the edge *set*, it is the **ordering** in
`edgesList.consensus.rawAvg.txt` — the 80 edges both methods found rank far above
the 127 only one of them found.

**A3.2 — strong negative `ecorr`.** In a microbiome network: two taxa whose
abundances move in opposite directions across samples — competitive exclusion,
mutually exclusive niches, or two states of the same community type. In the Amazon
co-purchase network: an edge exists because the products share a co-purchase
context, while the negative sign says the two are rarely bought by the same
customers — substitutes rather than complements. In both cases: **the edge comes
from MI (non-linear, unsigned), the sign comes from Spearman.** They answer
different questions and can legitimately disagree with intuition.

**A4.1 — the sweep.** ARACNE only:

| `var_perc` | `m` | subgraphs | edges @0.3 | edges @0.8 |
|---|---|---|---|---|
| 0.1 | 5 | 24 | 117 | 111 |
| 0.2 | 9 | 24 | 205 | 168 |
| 0.4 | 18 | 24 | 181 | 139 |

**The edge count is not monotone.** From `m=5` to `m=9`, more pairs get
co-selected at all, so more pairs can clear the threshold — the count rises. From
`m=9` to `m=18`, the presence-frequency *denominator* (`#subgraphs where both
variables were selected`) grows faster than the numerator: bigger subgraphs make
co-selection common, while ARACNE's DPI prunes more aggressively inside a denser
subset, so frequencies fall and edges drop out. There is an interior optimum, and
`0.2` is a reasonable default rather than a floor to raise blindly.

The subgraph count stays at 24 because it is `2·(k−1)`, and `k` comes from the
elbow of the eigenvalue curve — which depends on `W` and the Laplacian, not on
`var_perc`. `--var-perc` changes subgraph *size*; `--eigen-perc` changes subgraph
*count*.

**A4.2 — why reruns are byte-identical.** Two reasons, both deliberate. The rayon
fan-outs are order-independent: MI writes into a fixed pair-indexed structure, and
subgraph results are keyed by pair index and merged associatively, so thread
interleaving cannot change the result. And hill-climbing's random restarts are
seeded **per subgraph index** (`seed + idx`), not from a shared RNG, so each
subgraph draws the same restarts whichever thread runs it and whenever it runs.
Verified: `diff` on the consensus and both per-method edge lists across two runs
reports no difference.

**A5.1 — top ATC vs top betweenness.** They usually differ, and the reason is the
definition. Betweenness asks "how many shortest paths cross this node?" — a purely
topological question. ATC asks "how much did this node's *annotation mix* keep
changing while the annotation diffused?" — high for nodes sitting at the frontier
between themes, low for nodes buried inside one homogeneous theme *even if they
are structural hubs*. A hub whose neighbours all share its theme has high
betweenness and low ATC. That decoupling is the whole point of an
annotation-relative centrality.

**A5.2 — reporting modules honestly.** Repeat the run (5–10 times), report the
*range* of module counts and Q, and only name a module if its membership is stable
across repeats. Better: report enrichment signatures rather than module ids — "a
module enriched for X" survives a repartition in a way that "module 3" does not.
Until [B7](track-dev.md#b7--tests-and-validation-against-r-30-min) is fixed,
treat any single partition as one sample, not as the answer.

**A6.1 — scale-free threshold on Amazon.** Report γ, R² and the node fraction the
tool prints. Accept the threshold when γ lands in 2–3 **and** R² ≥ 0.80 **and** it
keeps a defensible share of nodes; a fit that is beautiful because it kept 12 % of
the network is not a threshold, it is a subsample. When no threshold satisfies the
criteria, say so — many real networks are not scale-free, and the fit failing is
information.

**A6.2 — Louvain on the filtered view vs the full graph.** They answer different
questions. On the **full graph**, modules are a property of the reconstructed
network and are comparable with anyone else's run on the same graph. On the
**filtered view**, modules are a property of *your* analysis — the structure that
survives your `|ecorr|` cut and exclusions. Report the filtered-view modules when
your claim is about the finalized topology, and always report the filter alongside
them. Never compare a filtered partition with someone else's unfiltered one.

**A6.3 — ATC before and after recompute.** The stored `atc` was computed once on
the raw full network, so a broad, weakly connected hub (host DNA, a dietary
contaminant, a catch-all class) touches every theme and dominates the score. Once
you exclude it and recompute on the visible sub-network, the propagation runs on
the topology you are actually interpreting, and nodes that genuinely bridge two
biological themes rise to the top. Two sentences to a colleague: *"The published
ATC ranks a contaminant first because it connects to everything. Recomputed on the
cleaned network, the top nodes are the taxa that actually sit between phyla."*

**A7.1 — ATC on a different annotation.** The ranking changes, sometimes
completely, while the network does not change at all. That is the correct
behaviour and the thing to internalize: **ATC is not a property of a node, it is a
property of a node *with respect to an annotation*.** Reporting "node X has high
ATC" without naming the annotation is meaningless.

---

## Track B

**B0.1 — the five stages in `pipeline.rs`.** `scalenet::prepare` covers stages 1–3
(`mi::build_mim` → affinity, `spectral::compute` → embedding and `k`, subset
selection → the `2·(k−1)` subgraphs); `scalenet::run_method` covers stage 4 (local
reconstruction) and stage 5's gather/threshold; `scs::consensus` covers the
fusion. `scalenet::fill_ecorr` between them attaches the Spearman sign, and
`data.rank_matrix()` is hoisted above the method loop so the ranks are computed
once.

**B1.1 — `fill_ecorr` after thresholding.** At `p = 2417` there are ~2.92 M pairs
but only tens of thousands of present edges. Running it before the cut would do
~2.9 M × O(n) Pearson evaluations on `n = 9817`-long rank vectors — ~10¹⁰
operations — to compute a value that is discarded for every absent pair.

**B1.2 — the rank precompute.** [performance.md](../performance.md#fixes-applied):
*"ecorr precompute + parallel — rank each column once, then O(n) Pearson on ranks
over present edges: 152 s → 53 s, output bit-identical."* Ranking is O(n log n)
per column; doing it once per column instead of twice per edge is the win.

**B2.1 — round-trip test.**

```rust
#[test]
fn pair_index_round_trips() {
    for p in [45usize, 1000] {
        for i in 0..p {
            for j in (i + 1)..p {
                assert_eq!(pair_at(pair_index(i, j, p), p), (i, j));
            }
        }
    }
}
```

**B2.2 — `NaN` inside, `Option` outside.** `EdgeRec` is the hot internal record,
touched in the innermost loops of gather, `fill_ecorr` and consensus. `f64` with
`NaN`-as-NA keeps it `Copy` and branch-free, and it matches the R output
convention directly. `Edge` is the public type, where "missing" must be
impossible to confuse with a number — a consumer that forgets to check `is_nan()`
gets silent garbage, while one that forgets to unwrap an `Option` does not
compile. Cheap inside, safe outside; the conversion happens once, in
`Network::from_edge_list`.

**B3.1 — ATC vs betweenness ranking.** Same reasoning as A5.1: they rank different
things. Pass `Some((&labels, &themes))` from `io::read_annotation_labels` into
`analyze` and sort by `an.atc`.

**B3.2 — why `analyze` needs `n_nodes`.** A `Network` carries only its **edges**.
Isolated nodes (degree 0 after thresholding) appear in no edge and would vanish,
which would silently shrink every denominator — density, closeness normalization,
modularity, and the ATC membership matrix. The node count has to come from the
dataset, not be inferred from the topology.

**B4.1 / B4.2 — the `correlation` method.** Expect the consensus edge set to grow
to the **union** of ARACNE's edges and yours (see A3.1), while the *ranking* in
`edgesList.consensus.rawAvg.txt` stays sensible only if `MethodInfo` describes
your method truthfully. With `weight_ecorr: true` on a method whose `score` is
already `|rho|`, you rank it by the Spearman correlation attached at
`fill_ecorr` time rather than by the score it produced — for this particular
method the two are nearly the same quantity, which is exactly why it is a good
teaching case: the fusion degrades quietly, and only a diff of the rawAvg file
shows it. On a method whose score is an unrelated scale (a BIC, a p-value), the
same mistake silently destroys the ranking.

**B5.1 — wasm vs CLI edge counts.** They must match. `par.rs` changes only
*scheduling*, never arithmetic or iteration order over the data, and the pipeline
is order-independent by construction. If the counts differ, one of the two B1
invariants is broken — most likely something in your new code accumulates in a
container whose iteration order differs between the two builds.

**B5.2 — why `graphjson` is its own feature.** The wasm build needs the rich
`graph.json` builder but must not pull in axum/tokio/tower-http — a server stack
that does not compile to `wasm32-unknown-unknown` and would bloat the module.
Splitting `graphjson` out lets `web = [parallel, graphjson, axum, tokio, …]` and
`wasm = [graphjson, wasm-bindgen]` share the builder with no common baggage.

**B6.1 — converging the two `graph.json` producers.** `graphjson::build` currently
takes `(&Dataset, &RunConfig, …)` and runs the pipeline **itself**, so the CLI
cannot call it without reconstructing twice. The smallest honest change is to
split it into `build_from(result: &PipelineResult, an: &NetworkAnalysis, …)` plus a
thin wrapper that keeps today's signature. Then the CLI can call `build_from`
with the result it already has. Two things `analyze` does not currently return and
that the rich form needs: the **full annotation table** (`analyze` reduces it to a
single label vector, while the rich form carries every column per node) and the
**`atc_trace`** (captured by `AtcConfig::capture_trace`, currently not requested by
`analysis::analyze`). `subgraphs` is already on `PipelineResult`, and the colour
map is derivable from the annotation column.

**B7.1 — fixing Louvain.** Replace the `HashMap` scan that picks the best
community with a deterministic order — collect `w_to` into a `Vec`, sort by
community id, then scan — and make the coarsened adjacency deterministic
(`BTreeMap`, or sort after collecting). Then assert stability:

```rust
#[test]
fn louvain_is_deterministic() {
    let g = /* fixture graph */;
    let (c0, q0) = /* metrics::compute(&g) */;
    for _ in 0..5 {
        let (c, q) = /* metrics::compute(&g) */;
        assert_eq!(c, c0);
        assert!((q - q0).abs() < 1e-12);
    }
}
```

On `main` this fails only *across processes*, because Rust's hasher is seeded per
process — so make the test compare against a **committed expected partition**, not
just against another call in the same process.

**B7.2 — why the existing test misses it.** `bayes_hc_is_deterministic` calls
`scalenet::run_method` twice and compares the two `EdgeList`s. It never touches
`metrics::compute`, so Louvain is outside its reach — and both calls happen inside
**one process**, sharing one hasher seed, so even a `HashMap`-order bug in the
reconstruction path could survive it. Cross-process determinism needs a committed
golden, which is exactly what `aracne_matches_golden` does and what the metrics
layer lacks.

**B8.1 — `--iterative-eigen`.** On `pop2mat` (`p = 45`) the whole run is under a
second, so expect no measurable speedup — the point of the exercise is the second
number: the consensus edge lists **differ**, because a different solver returns
different (still valid) eigenvectors, hence different subsets, hence a different
network. Report both, and note that this is why the flag is off by default and why
it requires a fixed `--eigen-perc` (the auto-`k` elbow needs the full spectrum).

**B8.2 — the dominant stage at `p = 5000`.** The dense eigendecomposition, O(p³) ≈
1.25 × 10¹¹ and single-threaded, while MI is O(p²·n) but parallel and reconstruction
is parallel over subgraphs. Before optimizing it, **measure the actual per-stage
wall time on your matrix** — the `n` of your data decides whether MI or the
eigensolver wins, and at `p = 2417` peak memory turned out to be dominated by
input parsing, not by anything the cost model highlights. Measure first; the
codebase's own optimization history is a list of surprises.
