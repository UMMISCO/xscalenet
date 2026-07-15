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
