# The IntegraNet suite — organization & exchange format

xscalenet is one component of a small **network-biology platform**. This document
defines how the pieces fit together and the **`graph.json` interchange format**
that lets them interoperate — the lingua franca of the suite.

## Components

```
IntegraNet  (platform: analysis + annotation + propagation + UI)
├── inference   scalenet (R, reference)  ·  xscalenet (Rust, scalable engine)
├── analysis    ATC propagation · scale-free · modularity · enrichment · theme layer
└── ui          xscalenet-web (local server) + network browser + propagation player
```

| Piece | Language | Role |
|---|---|---|
| **scalenet** | R | Reference implementation of the reconstruction + SCS. |
| **xscalenet** | Rust | Scalable engine: `reconstruct`, `analyze`, ATC, scale-free, metrics. CLI + `pyo3` bindings + lib API. |
| **xscalenet-web** | Rust (axum) | Local server: `POST /run` reconstructs a raw matrix → rich `graph.json`; serves the frontend. |
| **network browser** | JS (static) | `web/network.html` — explore/annotate/analyse a `graph.json`. |
| **propagation player** | JS (static) | `web/index.html` — step-by-step ATC on small annotated graphs. |

Everything speaks **`graph.json`**. Producers: the CLI (`--json`), the `pyo3`
bindings, `examples/graph_atc` (add ATC to an existing graph), `xscalenet-web`
(`POST /run`), and adapters like `.private/tools/scapis_to_graphjson.py` (from R
output). Consumer: the web apps. This decoupling means any inference tool that
emits `graph.json` plugs into the whole UI + analysis stack.

## `graph.json` — the interchange format

A single self-contained object. All fields except `nodes`/`edges` are optional;
consumers degrade gracefully (e.g. no `color_map` → palette by field; no
`atc_trace` → no propagation UI). `null` encodes a missing value.

```jsonc
{
  "meta": {
    "title": "…", "source": "…",
    "directed": false,
    "color_by": "phylum_group",       // default node-colour field (a node key)
    "n_nodes": 426, "n_edges": 1245
  },

  // node-colour legend for `meta.color_by` (categorical)
  "color_groups": ["Bacillota", "Bacteroidota", …],       // order = legend order
  "color_map":   { "Bacillota": "#3fb950", … },           // group → hex

  // ATC theme axis (the non-missing subset of the colour-by groups); theme
  // colours MUST match color_map so legend / fill / pies agree.
  "themes":       ["Bacillota", "Bacteroidota", …],
  "theme_colors": ["#3fb950", "#2f6df6", …],

  "nodes": [
    { "id": "msp_0009",
      // any number of annotation columns (taxonomy ranks, functional flags, …)
      "phylum": "Bacillota", "genus": "Enterocloster", "family": "Lachnospiraceae",
      "butyrate": "1",
      // topology (from metrics::compute)
      "degree": 23, "betweenness": 0.033, "closeness": 0.54, "clustering": 0.33,
      "community": 1,
      // ATC scores (present when annotated)
      "atc": 0.0512, "atc_last": 0.069,
      // optional fixed layout (e.g. exported from an igraph layout matrix)
      "x": -0.24, "y": 10.85 }
  ],

  "edges": [
    { "source": "msp_0009", "target": "msp_0018",
      "ecorr": 0.33,          // signed Spearman — drives edge colour & filtering
      "score": 0.30,          // presence score
      "orient": 2,            // 1 forward · -1 backward · 2 undirected
      "orient_score": 0.5 }
  ],

  // step-by-step ATC propagation, for the player / in-context animation
  "atc_trace": {
    "annotation": "phylum_group",
    "n_iterations": 9,
    "steps": [
      // step 0 = initial config T⁰; coinertia/delta null
      { "step": 0, "coinertia": null, "config": [[…themes…], …nodes…], "delta": null },
      { "step": 1, "coinertia": 0.556, "config": [[…], …], "delta": [0.25, …] }
    ]
  }
}
```

Notes:
- `config[i]` is node *i*'s theme-weight row `Tᵖ` (the pie); `delta[i]` is
  `‖rowᵢ(Tᵖ) − rowᵢ(Tᵖ⁻¹)‖₂`. A node's running ATC is
  `Σₚ ‖Δrowᵢ‖·(1−coinertiaₚ)`, converging to `atc`.
- Values are rounded to 6 decimals by producers to keep the file small.
- Client-side features that don't need the engine (Louvain modules, Fisher
  enrichment, taxonomic-layer aggregation, layouts, filters) run purely on this
  object, so they work on **any** conforming graph.

## Conventions

- **Naming.** IntegraNet = platform / umbrella; `xscalenet` = the scalable engine
  sub-brand; `scalenet` stays the R reference. New tools are new `graph.json`
  producers/consumers, not forks.
- **Interop contract.** A tool joins the suite by reading or writing the format
  above — nothing else. Keep the R edge-list schema (`x y epresenceScore …
  ecorr`) mappable to it (a thin adapter, e.g. the SCAPIS converter).
- **Repo shape (proposed).** A meta-repo / GitHub org linking the components with
  a shared docs site and this spec as the versioned contract; each tool stays in
  its own language.

See also [integration.md](integration.md) (library API + `pyo3`), and the web
apps' own [../web/README.md](../web/README.md).
