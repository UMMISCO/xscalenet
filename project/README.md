# project/ — the analysis workspace

One folder per **network** (dataset). xscalenet's analysis layer (ATC
propagation, the FunNet functional layer, enrichment / guilds, provenance cards)
is domain-agnostic — it runs on *any* annotated network — so each project here is
a different domain proving the same machinery.

```
project/
  projects.json        index (generated) — the app reads it for its network selector
  make_index.py        regenerates projects.json from every project.json
  <network>/
    project.json       manifest: title, domain, source, graph, color_by, ranks, note
    make_graph.py      converter: local source data → web/<network>.graph.json
    notes.md           what this network is / how it maps to the concept
```

The generated `graph.json` files live in `web/` (the served app root) and are
gitignored (large / client data), as is the local source data.

## Networks

| id | domain | nodes → | annotation |
|---|---|---|---|
| **scapis** | gut microbiome | MSP (species) | GTDB taxonomy + MSP traits · *client data* |
| **amazon** | e-commerce | products | group (Book/Music/DVD/Video) → category hierarchy |
| **flickr** | social / image | photos | dominant MIR tag |

## Add / rebuild a network

```sh
python3 project/<id>/make_graph.py            # → web/<id>.graph.json
cargo run --example graph_atc -- web/<id>.graph.json <color_by>   # add ATC (optional)
python3 project/make_index.py                 # refresh web/projects.json
cd web && python3 -m http.server              # network.html → pick it from "Networks ▾"
```
