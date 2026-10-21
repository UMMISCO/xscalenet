# Track A — Use xscalenet

**Audience:** analysts and biologists. No Rust needed; you will type `cargo build`
once and never look at Rust again.
**Outcome:** you can take a raw matrix to a reconstructed, annotated, filtered,
module-labelled network and defend every choice you made along the way.

Modules are cumulative. Answers: [solutions.md](solutions.md).

---

## A0 — Setup (15 min)

```sh
git clone <repo> && cd xscalenet
cargo build --release
cargo test  --release
```

`cargo test --release` must end with `4 passed` (regression) and `12 passed`
(unit). Those tests compare the engine against fixtures verified against the
original R package — a green suite means your build reproduces published results.

Then the local server, which you will need from module A6:

```sh
cargo build --release --features web --bin xscalenet-web
./target/release/xscalenet-web web
# xscalenet-web → http://localhost:8788   (serving ./web)
```

**Check:** `curl -s localhost:8788/health` → `ok`, and
<http://localhost:8788/network.html> opens with a network selector.

> **Exercise A0.1** — The server prints the folder it serves. Why does opening
> `web/network.html` directly from the file manager *not* work?

---

## A1 — The idea, in ten minutes (15 min)

The problem: reconstructing a network over `p` variables from `n` observations
breaks down when `p` grows. In the 2016 paper's microbiome case, Bayesian
hill-climbing on ~2 000 variables **did not converge in 48 hours**.

ScaleNet's answer is not a faster solver. It is a change of problem:

> Don't learn one big network. Learn many small, overlapping ones on
> **spectrally chosen** variable subsets, then aggregate.

Five stages, and you will see each of them in the output:

1. **Affinity.** Every variable pair gets a mutual-information score
   (Miller–Madow bias-corrected). MI, not correlation, because it catches
   non-linear dependence.
2. **Spectral embedding.** Build the random-walk Laplacian `Lrw = I − D⁻¹W` and
   eigendecompose it. The low eigenvectors encode "which variables belong to the
   same path-related component".
3. **Subsets.** For each retained eigenvector, take the `m` most positive entries
   (`G⁺,m`) and the `m` most negative (`G⁻,m`). With `k` eigenvectors kept you get
   `2·(k−1)` overlapping subgraphs of `m ≪ p` variables each.
4. **Local reconstruction.** Each subgraph is reconstructed independently and in
   parallel, by ARACNE (information-theoretic, undirected) and/or Bayesian
   hill-climbing (BDeu score, directed).
5. **Gather + consensus.** A pair's **presence frequency** is
   `#subgraphs where the edge was inferred / #subgraphs where both variables were selected`.
   Edges above the threshold survive. Across methods, a normalized-rank (Borda)
   fusion produces one consensus network.

The trick to internalize: **step 3 is what makes it scale, step 5 is what makes it
trustworthy.** Everything you tune in module A4 shifts one of those two.

Reference: [algorithm.md](../algorithm.md).

> **Exercise A1.1** — With `p = 45` and `var_perc = 0.2`, what is `m`? If the
> elbow heuristic keeps `k = 13`, how many subgraphs get reconstructed?

---

## A2 — Your first reconstruction (20 min)

```sh
./target/release/xscalenet tests/data/pop2mat.txt out --discretize
```

`pop2mat` is the fixture the R package ships with: 45 variables, 25 observations.
`--discretize` PAM-bins numeric columns, because MI needs discrete levels — it is
in the command because the R package's own `scs` example has it. Keep reading; on
*this* fixture it does nothing, and finding out why is exercise A2.2.

The run narrates itself on stderr:

```
# Discretizing numeric columns (pamk)...
# Loaded 45 variables x 25 observations
# Kept k=13 eigenvectors, 24 subgraphs, m=9 vars/subgraph
# Reconstructing with aracne...
# Reconstructing with bayes_hc...
# Consensus @ presFreq 0.3: 254 edges
# Consensus @ presFreq 0.8: 207 edges
# Done. Outputs under out
```

Read that block line by line — it is the algorithm of A1, in order. `k=13`,
`m=9`, `24 = 2·(13−1)` subgraphs: exactly the structure the R package produces,
and exactly what the regression test locks in.

The layout on disk mirrors the R package:

```
out/globalGraph/globalNet_presFreq0.3/edgesList.aracne.txt     # one method
out/globalGraph/globalNet_presFreq0.3/edgesList.bayes_hc.txt
out/consensusGraph/globalNet_presFreq0.3/edgesList.txt         # the fused network
out/consensusGraph/globalNet_presFreq0.3/edgesList.consensus.rawAvg.txt
                                        …presFreq0.8/…         # one folder per threshold
```

> **Exercise A2.1** — Run it again with `--thresh 0.3,0.5,0.8`. How many folders
> appear under `consensusGraph/`, and does the 0.3 network change?
>
> **Exercise A2.2** — Drop `--discretize` and rerun. The two consensus edge lists
> are **byte-identical**. Look at `head -3 tests/data/pop2mat.txt` and explain why.
> Then say what would have happened on a matrix of relative abundances, and which
> stage of A1 is responsible.

---

## A3 — Reading the outputs (25 min)

### The edge list

```
x	y	epresenceScore	epresence	eorientScore	eorient	ecorr
GU:185_<<_GU:178	GU:185	GU:178	0.10115329911570095	1	-1	-1	0.9185586535436918
GU:185_<<_GU:1020	GU:185	GU:1020	NA	0	NA	NA	NA
```

R convention: the first field is an unnamed row name (`x_<<_y`), so the header
has one fewer name than the rows have fields. Then:

| Column | Meaning |
|---|---|
| `x`, `y` | the two variables, `x < y` in the internal ordering |
| `epresenceScore` | ARACNE presence score, or the consensus average rank |
| `epresence` | `1` present, `0` absent |
| `eorientScore` | orientation confidence (consensus only) |
| `eorient` | consensus: `1` forward · `−1` backward · `0` unoriented. Per-method files use `2 / −2 / 1` |
| `ecorr` | **signed Spearman correlation** of the two variables |

Two traps worth ten minutes of discussion:

- **Every pair is written, not just the edges.** `edgesList.txt` at `p = 45` has
  991 lines — `45·44/2 = 990` pairs plus the header — most of them `NA`. Use
  `--edges-only` to write present edges only. At `p = 2417` the full-pair form is
  ~2.9 M rows *per file*.
- **`ecorr` is not the edge's strength, it is its sign and its direction of
  co-variation.** Edge existence comes from MI (non-linear, unsigned); `ecorr`
  says whether the two co-occur (+) or exclude each other (−). Filtering by
  `|ecorr|` — which you will do constantly in the browser — is therefore a
  *linear-monotone* filter layered on a *non-linear* graph. That is deliberate and
  useful, but say it out loud when you present results.

### The consensus vs the per-method files

`edgesList.consensus.rawAvg.txt` is the same edges ranked best-first by average
normalized rank. When two methods disagree, this file is where you look: an edge
ARACNE ranked first and hill-climbing never inferred lands mid-table.

> **Exercise A3.1** — In `out/consensusGraph/globalNet_presFreq0.8/edgesList.txt`,
> count present edges (`epresence == 1`) with `awk`. Compare with
> `edgesList.aracne.txt` at the same threshold. Which is larger, and why must it be?
>
> **Exercise A3.2** — Find an edge with `ecorr < −0.8`. What does a strong negative
> `ecorr` mean biologically in a microbiome network, and what would it mean in the
> Amazon co-purchase network (module A6)?

---

## A4 — The parameters that matter (30 min)

Nine flags exist; four of them change your science. Learn those, leave the rest.

| Flag | Default | What it really controls |
|---|---|---|
| `--var-perc` | `0.2` | `m = ⌈var_perc·p⌉`, the size of each local problem. **Bigger = richer local structure, superlinearly more compute** (ARACNE is O(m³)). |
| `--eigen-perc` | `auto` | `k`, how many spectral directions are explored. `auto` = elbow heuristic on the eigenvalue curve. More `k` = more subgraphs = more coverage. |
| `--thresh` | `0.3,0.8` | The presence-frequency cut. **This is your precision/recall dial**, and it is free — pass several and compare. |
| `--methods` | `aracne,bayes_hc` | ARACNE alone is fast and R-exact. Adding `bayes_hc` costs ~6× at `p=250` and grows with `m`, and buys you orientation. |

The rest — `--epsilon`, `--restart`, `--maxp`, `--max-iter`, `--iss`,
`--max-clusters`, `--seed` — are method internals. Touch them when you have a
reason, and record the value when you do.

Two flags are in the CLI but not yet in the README:

- `--scale-free` — sweep `|ecorr|` thresholds, fit a power law at each, and
  recommend one (module A6 does the same interactively).
- `--iterative-eigen` — Lanczos eigensolver instead of the dense one. Fast path
  for large `p` with a **fixed** `--eigen-perc`; not bit-identical to the default.

### The experiment

Run the sweep and fill the table. This is the module's whole point: **numbers you
produced yourself beat any claim in this document.**

```sh
for vp in 0.1 0.2 0.4; do
  echo "== var_perc=$vp"
  ./target/release/xscalenet tests/data/pop2mat.txt out_vp$vp \
    --discretize --methods aracne --var-perc $vp 2>&1 | grep -E "Kept|Consensus"
done
```

| `var_perc` | `m` | subgraphs | edges @0.3 | edges @0.8 |
|---|---|---|---|---|
| 0.1 | | 24 | | |
| 0.2 | 9 | 24 | | |
| 0.4 | | 24 | | |

(`pop2mat` runs in well under a second at every setting, so don't bother timing it
here — timing is module B8's job, on a matrix big enough to measure.)

> **Exercise A4.1** — Fill the table. The edge count is **not monotone** in
> `var_perc`: it rises from 0.1 to 0.2 and falls again at 0.4. Explain both
> directions using the presence-frequency *denominator* (`#subgraphs where both
> variables were selected`). Why does the subgraph count stay at 24?
>
> **Exercise A4.2** — Run twice with identical arguments and `diff` the consensus
> edge lists. They are identical. Which design choice guarantees that even though
> subgraphs are reconstructed in parallel and hill-climbing uses random restarts?

---

## A5 — The functional layer (30 min)

Reconstruction gives you a topology. The functional layer tells you what it
*means*. It is opt-in, and it needs an annotation to do its best work.

```sh
# topology only
./target/release/xscalenet tests/data/pop2mat.txt out --discretize --metrics

# with an annotation: ATC, per-module enrichment, the theme layer, and graph.json
./target/release/xscalenet tests/data/pop2mat.txt out --discretize --metrics \
  --atc annot.tsv --atc-col group --theme-network --json
```

The annotation file is a TSV whose **first column is the node id** and whose other
columns are attributes; `--atc-col` picks one (default: column 2). Ids that do not
match a node are ignored; nodes with no row get no theme.

Four artefacts appear next to the edge list:

**`nodeAnalysis.txt`** — one row per node:
`degree · strength · betweenness · closeness · clustering · community ·
modular_degree · modular_betweenness · atc · atc_last`.

The **modular** columns are the Prifti 2010 distinction that most tools miss:
`betweenness` is a node's brokerage in the *whole* network, `modular_betweenness`
is its brokerage *inside its own community*. A node can be globally peripheral and
locally central — a module's internal hub — and only the pair of numbers shows it.

**`enrichment.txt`** — for each Louvain community, which annotation values are
over-represented: one-sided Fisher exact test + Benjamini–Hochberg `qvalue`. This
answers "is this module actually *about* something?".

**`themeNetwork.txt`** — the FunNet second layer. Annotations are propagated over
the network to convergence, each theme gets a node-profile, and themes are linked
by a self-tuning Gaussian kernel over profile distances. **Two themes are close
when their footprints diffuse to the same regions of the network** — an indirect,
higher-order relation, not a count of direct contacts between their members.

**ATC** (`atc`, `atc_last`) — Annotation Transcriptional Centrality. The
annotation is propagated as a dynamical system until the co-inertia of successive
configurations converges; a node's ATC accumulates how much its theme mix kept
changing along the way. Nodes that sit *between* themes score high; nodes deep
inside one homogeneous theme score low. It is a centrality **relative to an
annotation**, which is why the same network gives different ATC for taxonomy vs
function.

### ⚠ The caveat to state before anyone builds on modules

Reconstruction is deterministic; **Louvain community detection is not**. Four
identical runs on `pop2mat` @0.3:

```
254 edges, 4 communities, modularity 0.343
254 edges, 3 communities, modularity 0.329
254 edges, 3 communities, modularity 0.337
254 edges, 4 communities, modularity 0.343
```

The edge set never moves; the partition does. Consequences: `community`,
`modular_degree`, `modular_betweenness` and every enrichment `qvalue` can shift
between runs. **Do not report a module count as a finding without repeating the
run.** (Cause and fix: [B7](track-dev.md#b7--tests-and-validation-against-r-30-min).)

> **Exercise A5.1** — Build a 2-column annotation from `nodeAnalysis.txt` (node id
> + any grouping you like), rerun with `--atc --theme-network`, and find the node
> with the highest ATC. Is it also the highest-betweenness node? Why usually not?
>
> **Exercise A5.2** — Run the same command three times and diff the `community`
> column of `nodeAnalysis.txt`. How would you report modules honestly given what
> you see?

---

## A6 — The network browser (45 min)

This is where the analysis actually happens. Server running from A0:

<http://localhost:8788/network.html>

Pick a bundled network from **Networks ▾** — start with **ATC 10-node reference**
(small enough that every node's pie chart is legible), then move to **Amazon
co-purchase** or **Flickr** for scale.

### Import your own matrix

**⬆ Import data** → pick a matrix (variables in columns, samples in rows, or tick
transpose) → choose methods/threshold → optional annotation table + colour-by
column → **Reconstruct**. It runs the same validated Rust engine, locally,
multi-threaded. The equivalent from a terminal, which is worth doing once so you
know what the button does:

```sh
curl -X POST http://localhost:8788/run \
  -F "data=@tests/data/pop2mat.txt" \
  -F "annotation=@annot.tsv" -F "color_by=group" \
  -F "methods=aracne" -F "thresh=0.3" -F "discretize=true" \
  -o imported.graph.json
```

The response is a **rich** `graph.json`: `meta` (with `k` and `subset_m`),
`color_map`, `themes`, per-node annotation + metrics + ATC, the consensus `edges`,
the full `atc_trace`, and `subgraphs` — the actual spectral subsets that were
fused. An imported network is first-class: every tool below works on it.

### The workflow that matters

Do these in order. Each step changes what the next one means.

1. **Colour by** an annotation, and toggle legend entries off to drop what you do
   not want (contaminants, host DNA, an uninformative catch-all class).
2. **Edge filters** — min `|ecorr|`, sign, min presence score, top-N, min degree.
   Live, and they compose with everything else.
3. **Scale-free topology** — fit `P(k) ~ k^−γ` on the current view, read **γ, R²,
   k_max**, and let it **recommend a min-|ecorr| threshold**. One click applies it.
   This is the principled way to pick a filter instead of eyeballing density. The
   CLI does the same sweep with `--scale-free`; on `pop2mat` @0.3 it reports:

   ```
   -> recommended min|ecorr| ~ 0.901: gamma=0.44, R2=0.617, keeps 71% of nodes
   ```

   Read that critically: **R² = 0.617 is below WGCNA's 0.80 and γ = 0.44 is far
   outside Barabási's 2 ≤ γ ≤ 3.** A 45-node network is simply too small to be
   scale-free, and the tool says so rather than pretending. On a real network you
   expect the recommendation to land in the criteria; when it doesn't, that is a
   finding about your data, not a number to quote anyway.
4. **Modularity** — run Louvain on the full graph *or on the filtered view only*,
   optionally weighted by `|ecorr|`. Read **Q** and the module count, then read
   each module's **taxonomy composition bar + homogeneity %** and, with enrichment
   on, its **signature** (the top Fisher-enriched feature, e.g.
   "▲ Lachnospiraceae · family · q=3e-17").
5. **Recompute ATC on the current view** (↻ in the ATC panel). The stored `atc`
   was computed once on the raw full network. After you have excluded
   contaminants and set your filters, the stored value describes a network you are
   no longer looking at. Recomputing re-runs the propagation in-browser on exactly
   the visible sub-network. **This is the single most important habit in the
   browser** — a broad weak hub otherwise dominates ATC and hides the biology.
6. **ATC propagation** — play/scrub the timeline: node size grows with accrued
   ATC, node colour shows the evolving theme mix (blend or pie), with the per-step
   co-inertia. This is A5's dynamical system, animated on your network.
7. **Taxonomic / functional layers** — collapse to any rank as the FunNet
   functional layer, laid out by MDS of the profile distances, so spatial
   proximity ≈ functional proximity.
8. **Spectral subgraphs** — for a network you reconstructed here, scrub the
   `2·(k−1)` subsets that were fused, or hit **Overlap heatmap** to tint every
   node/edge by how many subgraphs cover it. The bright core is where the spectral
   embedding concentrated its attention. This is A1 step 3, made visible.
9. **Export** — PNG · Nodes CSV · Edges CSV (of the *filtered* view) · `graph.json`
   carrying your layout and modules.

### State is the deliverable

Filters, colouring, size-by, layout, modules, focus, propagation step, selection,
pan/zoom and node positions are saved to the **URL** and to `localStorage`. A
refresh restores your session; the address bar is a shareable, reproducible view:

```
network.html?ec=0.4&sign=pos&colorby=phylum&sizeby=atc_last&prop=4&atcrec=1
```

Copy the URL into your lab notebook alongside the command that produced the graph
and your analysis is reproducible end to end.

> **Exercise A6.1** — Load the Amazon network, colour by `group`, and use the
> scale-free recommendation to pick a threshold. Report γ, R², and the fraction of
> nodes kept. Would you accept this threshold?
>
> **Exercise A6.2** — On the same network, run Louvain on the *filtered view*, then
> on the *full graph*. Do the modules mean the same thing? Which do you report?
>
> **Exercise A6.3** — Find the top-ATC node before and after **↻ Recompute on
> current view** with a hub excluded. Explain the change to a colleague in two
> sentences.

---

## A7 — Reproducibility and sharing (20 min)

Three artefacts, and you should be able to say which one you are handing over.

- **The command.** Every parameter that changed your result, including `--seed`.
- **`graph.json`.** One self-contained object: nodes with annotations and metrics,
  edges with `ecorr`/`score`/`orient`, colours, themes, `atc_trace`, layout. It is
  the suite's lingua franca — the CLI, the Python bindings, the server, the wasm
  build and the converters all produce it, and the web apps consume it. Format spec:
  [suite.md](../suite.md).
- **The view URL.** What you were looking at when you drew your conclusion.

To organize several datasets, the `project/` workspace holds one folder per
network (`project.json` manifest + `make_graph.py` converter + `notes.md`), and
`project/make_index.py` regenerates `web/projects.json` — which is what fills the
**Networks ▾** selector.

```sh
python3 project/<id>/make_graph.py                              # → web/<id>.graph.json
cargo run --example graph_atc -- web/<id>.graph.json <color_by> # add ATC, no re-run
python3 project/make_index.py                                   # refresh the selector
```

Note `graph_atc`: it adds ATC to an **existing** graph using the validated engine,
without re-running reconstruction. Use it when you change your mind about which
annotation to propagate.

> **Exercise A7.1** — Take the `imported.graph.json` from A6, add ATC on a
> different annotation column with `graph_atc`, and load both into the browser.
> Which nodes change rank, and what does that tell you about ATC?

---

## Capstone

On a matrix of your choice (or `pop2mat`):

1. Reconstruct with two thresholds, justify the one you keep.
2. Produce `nodeAnalysis.txt` with an annotation; name the three highest-ATC nodes
   and say why they are high.
3. In the browser: pick an `|ecorr|` threshold via the scale-free fit, detect
   modules on the filtered view, and give each module a one-line signature from
   its enrichment.
4. Recompute ATC on the finalized topology and report what moved.
5. Hand in: the command, the `graph.json`, the view URL, and five sentences of
   findings — including one sentence on what the Louvain caveat means for your
   module claims.
