# Integration: library API + interchange format

xscalenet is usable as a **library** and produces a **self-contained JSON graph**
so a web frontend, notebook, or the IntegraNet app can consume its output without
re-parsing edge-list files.

## Python bindings

In-process bindings expose `reconstruct(input, ...)` returning the JSON below;
see [../bindings/python/README.md](../bindings/python/README.md). Built with
pyo3/maturin (abi3, one wheel for Python 3.8+).

```python
import json, xscalenet
g = json.loads(xscalenet.reconstruct("data.tsv", discretize=True, metrics=True,
              annotation="taxonomy.tsv", atc_col="genus", theme_network=True))
```

## Library API (Rust)

Two typed entry points, no I/O:

```rust
// 1. Reconstruct the network(s).
let result: PipelineResult = xscalenet::run(&data, &run_config, |_| {})?;

// 2. Analyse a network (metrics + ATC + enrichment + theme layer).
use xscalenet::analysis::{analyze, AnalysisConfig};
let net = &result.consensus[0].network;
let annot = Some((labels.as_slice(), theme_names.as_slice())); // node -> theme index
let a = analyze(net, result.nodes.len(), annot, &AnalysisConfig { weighted: false, theme_network: true });
// a.metrics (per-node degree/betweenness/community/modular_*), a.atc, a.enrichment,
// a.theme_network
```

`PipelineResult` and `NetworkAnalysis` are plain typed structs — the same values
the CLI formats and a `pyo3`/web layer would expose.

## JSON interchange (`--json` → `graph.json`)

A single object describing the network, its node-level analysis, and the
second-layer theme network:

```jsonc
{
  "nodes": [
    { "id": "GU:185", "community": 1, "annotation": "genus_1",
      "degree": 12, "strength": 12, "betweenness": 0.033, "closeness": 0.54,
      "clustering": 0.33, "modular_degree": 8, "modular_betweenness": 0.12,
      "atc": 0.0255, "atc_last": 0.41 }
  ],
  "edges": [
    { "source": "GU:185", "target": "GU:178",
      "orientation": "undirected", "score": 0.52, "ecorr": 0.91 }
  ],
  "theme_network": {                        // present only with --theme-network
    "themes": ["genus_0", "genus_1", "genus_2"],
    "edges": [ { "source": "genus_0", "target": "genus_1", "proximity": 0.91 } ]
  }
}
```

This is the shape a graph library (Cytoscape.js / Sigma.js) or IntegraNet needs:
node positions/colors from `community`/`annotation`, node size from a centrality
or `atc`, edge color from `orientation`/`ecorr`, and the theme layer as a second
graph. `null` encodes a missing value.

## Relation to IntegraNet

IntegraNet imports an edge list (`edge1`/`edge2`) plus an annotation table
(`id` + attributes). The `graph.json` above carries the same information in one
file (edges + per-node annotation + metrics), plus the functional layer
(enrichment, theme network) that IntegraNet's prototype computes separately — and
here computed by the validated Rust core. A thin adapter can feed `graph.json`
to the app, or emit the two TSVs from it.
