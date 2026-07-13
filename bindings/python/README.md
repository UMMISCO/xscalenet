# xscalenet — Python bindings

In-process Python access to the xscalenet engine (network reconstruction +
functional-layer analysis), built with pyo3/maturin.

## Build & install

```sh
pip install maturin
maturin develop --release -m bindings/python/Cargo.toml   # into the active venv
# or: maturin build --release -m bindings/python/Cargo.toml  (produces a wheel)
```

## Use

```python
import json, xscalenet

# Reconstruct a network (TSV or .gz; features in columns, obs in rows).
r = json.loads(xscalenet.reconstruct("data.tsv", methods="aracne,bayes_hc", discretize=True))
r["nodes"]                      # feature names
r["results"][0]["edges"]        # consensus edges at the first threshold

# Add the functional layer (metrics, ATC, per-module enrichment, theme network).
g = json.loads(xscalenet.reconstruct(
    "data.tsv", discretize=True, metrics=True,
    annotation="taxonomy.tsv", atc_col="genus", theme_network=True,
))["results"][0]
g["nodes"][0]                   # {id, community, annotation, betweenness, atc, ...}
g["theme_network"]["edges"]     # theme-proximity (second layer)
```

`reconstruct(input, methods, discretize, transpose, var_perc, eigen_perc, thresh,
epsilon, restart, seed, max_iter, iss, maxp, max_clusters, metrics, annotation,
atc_col, theme_network, weighted)` — keyword args mirror the CLI flags. Returns a
JSON string; parse with `json.loads`. The JSON schema is in
[../../docs/integration.md](../../docs/integration.md).
