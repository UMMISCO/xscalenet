# xscalenet web

Self-contained, dependency-free frontends for xscalenet output. No build step —
just static files.

## ATC propagation player — `index.html`

The real app. Animates the annotation-propagation dynamics (ATC, Prifti et al.
*Bioinformatics* 2010) **step by step**: each node is a pie chart of its theme
weights `Tᵖ`, the timeline scrubs through the propagation, the co-inertia curve
tracks convergence, and the bars rank nodes by ATC.

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

## network browser — `network.html`

Canvas view of a reconstructed network (`graph.json`) for exploring larger
graphs:

- **Colour by taxonomy / annotation** — a dropdown over any categorical node
  field (phylum, family, …) with a click-to-toggle legend.
- **Edge filters** — min `|ecorr|`, sign (both ± / positive / negative), min
  presence score, keep top-N by `|ecorr|`, and min node degree. All live.
- **Scale-free topology** — fits the degree distribution to a power law
  `P(k) ~ k^−γ` (log-log linear regression), reports **γ, R², k_max**, judges it
  against the scale-free criteria (Barabási's 2 ≤ γ ≤ 3; WGCNA's R² ≥ 0.80), and
  **recommends a min-|ecorr| threshold** that maximises the power-law fit while
  keeping a meaningful fraction of nodes. One click applies it.
- Hover for taxonomy + degree; drag nodes, scroll to zoom, pan the background.

Shareable views via URL query: `network.html?ec=0.4&sign=pos&colorby=phylum`.

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
