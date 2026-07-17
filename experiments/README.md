# experiments — the concept on any annotated network

xscalenet's analysis layer (ATC propagation, the FunNet functional layer,
enrichment / guilds, the provenance cards) is **not microbiome-specific**. It
applies to *any* annotated network: nodes with a categorical annotation
(optionally a hierarchy) and edges with a weight/sign. These experiments show it
on classic SNAP graphs.

Each dataset has a `make_*_graph.py` that samples a visualizable subgraph from the
(large, local, gitignored) source data and writes a `graph.json` into `web/`
(also gitignored). The converters are committed; the data and generated graphs
are not.

| dataset | nodes = | annotation (coarse → fine) | edges = |
|---|---|---|---|
| **amazon** (SNAP co-purchase) | products | `group` (Book/Music/DVD/Video) → category hierarchy | co-purchased, weight = category overlap |
| **flickr** (SNAP image graph) | photos | (node features) | shared-context relationships |

`meta.ranks` in the graph.json declares the annotation hierarchy so the web app's
taxonomic-layer ladder works on *any* domain (not just species→phylum).

Run e.g.:
    python3 experiments/amazon/make_amazon_graph.py
    cargo run --example graph_atc -- web/amazon.graph.json group   # add ATC
    cd web && python3 -m http.server   # open network.html?graph=amazon.graph.json
