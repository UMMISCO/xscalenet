# xscalenet web

Self-contained, dependency-free frontends for xscalenet output. No build step —
just static files. Two **independent** tools (they don't share a dataset):

- **`network.html`** — browse/analyse a real reconstructed network (e.g. SCAPIS),
  with in-context ATC propagation. This is the main app.
- **`index.html`** — a standalone, dataset-agnostic player for the propagation
  *algorithm* on **small** annotated graphs (default: the 10-node R-validated
  reference). Best when per-node pie charts are legible.

## ATC propagation player — `index.html`

A standalone, dataset-agnostic player for the annotation-propagation dynamics
(ATC, Prifti et al. *Bioinformatics* 2010), for **small** annotated graphs.
Animates it **step by step**: each node is a pie chart of its theme weights
`Tᵖ`, the timeline scrubs through the propagation, the co-inertia curve tracks
convergence, and the bars rank nodes by ATC. It loads **any** graph.json with an
`atc_trace`; the bundled default is the 10-node reference. For large real
networks, use the network browser instead (it animates the same propagation).

```sh
# 1. generate the bundled demo (the R-validated 10-node test network)
cargo run --example propagation_demo        # -> web/propagation.graph.json

# 2. serve this folder and open it (fetch is blocked on file://)
cd web && python3 -m http.server
#   http://localhost:8000/index.html
```

Play/scrub the timeline, hover a node for its per-step weights, or drag /
scroll-zoom the network. Use **Load graph.json** to drop in any graph that
carries an `atc_trace` block (see below). `?step=N` in the URL jumps to a step.

The demo is the exact graph `J` + membership `m0` that xscalenet validates
against R (`src/atc.rs::matches_r_dynamic_propagation`): 10 nodes, 3 themes
(brown / black / yellow), converging at step 2.

### Exports

**⬇ Export** (header) → **PNG** of the current view · **Nodes CSV** (every
annotation column + degree/betweenness/community/module/ATC) · **Edges CSV** of
the *filtered* view · **graph.json** carrying the current layout + module
assignment (re-loadable to restore both).

## reconstruction server — `xscalenet-web`

To reconstruct a network **from a raw matrix** inside the browser, run the local
server (it embeds the validated Rust engine and serves this folder):

```sh
cargo run --features web --bin xscalenet-web        # serves ./web on http://localhost:8788
```

Open <http://localhost:8788/network.html>, click **⬆ Import data**, pick a matrix
(variables in columns, samples in rows — or tick transpose), choose params and an
optional annotation table, and hit **Reconstruct**. Everything is local; large
matrices use the full multi-threaded engine.

Imported data is **first-class**: the server returns a *rich* `graph.json` — every
annotation column carried per node, `color_map` + `themes`, per-node metrics
(degree/betweenness/community) and ATC, and the full `atc_trace`. So an imported
network gets the whole toolset (taxonomy colouring, taxonomic layers, module ×
taxonomy, ATC propagation with blend/pie), exactly like the bundled SCAPIS graph.

`POST /run` accepts a multipart form (`data` file + `methods`, `thresh`,
`transpose`, `discretize`, `weighted`, … fields; optional `annotation` file +
`color_by` column — defaults to a taxonomy rank).

When the page is served by a plain static server instead (`python3 -m
http.server`), the reconstruction button is inert (it detects `/health`); use
**Load graph.json** with a pre-built graph in that case.

## network browser — `network.html`

Canvas view of a reconstructed network (`graph.json`) for exploring larger
graphs:

- **Colour by taxonomy / annotation** — a dropdown over any categorical node
  field (phylum, family, …) with a click-to-toggle legend.
- **Select a node & act on it** — click a node to *pin* it (Shift/⌘-click for
  multi-select; click empty / Esc to clear). The Node panel then shows its full
  taxonomy, degree, ATC, and a **signed partner list** (neighbours sorted by
  `ecorr`, + co-occurring / − mutually-excluding), each clickable. Actions:
  **Focus ego-network** (isolate the node + its neighbours — a floating chip
  shows/clears it; **Expand +1** grows it), **Select <phylum>** (all same-taxon
  nodes), **Copy** (info to clipboard), **Export ⬇** (the ego-network as a
  `graph.json`). Select **two** nodes → the **shortest path** between them (fewest
  hops over the visible graph) is highlighted on the canvas and listed with
  taxonomy. Deep-link: `?select=<id>[,<id2>]&focus=1`.
- **Edge filters** — min `|ecorr|`, sign (both ± / positive / negative), min
  presence score, keep top-N by `|ecorr|`, and min node degree. All live, and
  compose with the ego-network focus.
- **Scale-free topology** — fits the degree distribution to a power law
  `P(k) ~ k^−γ` (log-log linear regression), reports **γ, R², k_max**, judges it
  against the scale-free criteria (Barabási's 2 ≤ γ ≤ 3; WGCNA's R² ≥ 0.80), and
  **recommends a min-|ecorr| threshold** that maximises the power-law fit while
  keeping a meaningful fraction of nodes. One click applies it.
- **Modularity** — detect communities with **Louvain** (multi-level, optionally
  weighted by `|ecorr|`, on the full graph or **the filtered view only**);
  reports the modularity **Q** and #modules, colours the network by module, and
  lists modules with a **taxonomy composition bar + homogeneity %** (cross
  modules × taxonomy — are the communities taxonomically pure?), and — with
  **enrichment** on — each module's **signature**: the most over-represented
  annotation feature (any taxonomy rank or functional column) by a one-sided
  **Fisher exact test + Benjamini-Hochberg** (a phylum-pure Bacillota module
  resolves to e.g. "▲ Lachnospiraceae · family · q=3e-17"). Click a module to
  isolate it. `?modules=1` deep-links a coloured partition.
- **Taxonomic layers** — collapse the species network to any rank (genus →
  phylum): each node becomes a taxon at its members' layout centroid, sized by
  number of species, with between-taxa edges (width ∝ number of species-edges).
- **Layout** — a selector: *original (igraph)* when nodes carry `x`/`y` (e.g.
  exported from the R igraph object's layout matrix), *force-directed*, *circular
  (by phylum)*, *concentric (by degree)*, and *grid (by phylum)*.
- **Simulate** (⚛) — a **live**, **Barnes-Hut** (O(n·log n)) force simulation on
  the *visible* subgraph — stays smooth on large imported networks (≈12× faster
  than exact repulsion at 15k nodes):
  turn it on once and the layout re-flows smoothly as you filter (no re-layout
  button per change, no teleporting); dragging a node nudges its neighbours;
  **Freeze** when happy. This is the ergonomic way to "re-lay-out after
  filtering." (`?sim=1` / `?layout=` for shareable views.)
- **Size by** degree, taxon size, or ATC (dist·co-inertia / last-config), when
  available.
- **ATC propagation** — when the graph has an `atc_trace` (see below), a
  play/scrub timeline animates the annotation propagation *on the real network*:
  node **size** grows with the ATC accrued so far, and node **colour** shows the
  evolving theme mix `Tᵖ` — either a **colour blend** or a **pie chart** per node
  (toggle) — with the per-step co-inertia. Only applies the theme colouring when
  you're colouring by the annotation; any other colouring (e.g. module) is kept.
- Hover for taxonomy + degree + ATC; drag nodes, pan the background, and zoom by
  scroll, **two-finger pinch** (trackpad or touch), or Safari pinch gesture.

Shareable views via URL query:
`network.html?ec=0.4&sign=pos&colorby=phylum&sizeby=atc_last&prop=4`.

### Adding ATC to a network

`graph.json` gains ATC from the validated engine — no re-run of the pipeline:

```sh
cargo run --example graph_atc -- web/scapis.graph.json phylum_group
```

This computes ATC on the graph's own edges, annotated by the chosen node field
(here `phylum_group`), and writes back per-node `atc` / `atc_last`, the `themes`
+ `theme_colors`, and an `atc_trace` (the same shape the propagation player
reads). Scores match the CLI / R reference.

```sh
cd web && python3 -m http.server      # http://localhost:8000/network.html
```

The bundled SCAPIS demo (`scapis.graph.json`) is **client data and is
gitignored** — generate it from the R analyses with
`.private/tools/scapis_to_graphjson.py`, or drop in any `graph.json` via
**Load graph.json**.

## graph previewer — `viewer.html`

A minimal previewer for a plain `--json` `graph.json` (nodes + edges + metrics +
theme network), kept for quick eyeballing of engine output. Superseded by
`index.html` for anything propagation-related.

## `atc_trace` schema

Produced by the demo generator; consumed by `index.html`.

```jsonc
{
  "meta":   { "title": "…", "source": "…", "n_steps": 14, "converged_at": 2 },
  "themes": ["brown", "black", "yellow"],
  "theme_colors": ["#9c6b3f", "#2b2b2b", "#e3b341"],
  "nodes":  [ { "id": "g1", "community": 0, "degree": 2, "betweenness": 0.0,
               "atc": 0.1127, "atc_last": 0.4028, "m0": [1,0,0] } ],
  "edges":  [ { "source": "g1", "target": "g2" } ],
  "atc_trace": {
    "steps": [
      // step 0 is T⁰ (initial, column-normalized); coinertia/delta are null.
      { "step": 0, "coinertia": null, "config": [[…themes…], …nodes…], "delta": null },
      // step p>0: config Tᵖ, RV co-inertia vs Tᵖ⁻¹, per-node L2 change.
      { "step": 1, "coinertia": 0.5136, "config": [[…], …], "delta": [0.25, …] }
    ]
  }
}
```

`config[i]` is node *i*'s theme-weight row at that step (the pie); `delta[i]` is
`‖rowᵢ(Tᵖ) − rowᵢ(Tᵖ⁻¹)‖₂`. The trace is captured by `AtcConfig::capture_trace`
in the Rust core.
