#!/usr/bin/env python3
"""The R-validated 10-node ATC reference network, as a browsable project.

Enriches web/propagation.graph.json (the propagation-player demo) with a single
categorical annotation — each node's dominant theme (brown/black/yellow from its
m0 membership) — plus color_by/color_map, so it loads in the network browser and
colours/propagates like any other project. Small & committed (not client data).
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "web", "propagation.graph.json")
OUT = os.path.join(HERE, "..", "..", "web", "testnet.graph.json")

def main():
    g = json.load(open(SRC))
    themes = g.get("themes", []); tcol = g.get("theme_colors", [])
    for n in g["nodes"]:
        m0 = n.get("m0") or []
        if m0 and max(m0) > 0:
            n["theme"] = themes[m0.index(max(m0))]
        else:
            n["theme"] = "unannotated"
    groups = themes + (["unannotated"] if any(n["theme"] == "unannotated" for n in g["nodes"]) else [])
    cmap = {t: c for t, c in zip(themes, tcol)}
    cmap["unannotated"] = "#c9ced6"
    g["meta"] = {"title": "ATC 10-node reference network",
                 "source": "Prifti et al. 2010 dynamicPropagation (R-validated)",
                 "directed": False, "color_by": "theme", "annotation_source": "reference m0"}
    g["color_groups"] = groups
    g["color_map"] = cmap
    # the propagation player recolours nodes only when colour-by == the trace's
    # annotation; tie the trace to the `theme` field so the animation shows.
    if g.get("atc_trace") is not None:
        g["atc_trace"]["annotation"] = "theme"
    json.dump(g, open(OUT, "w"))
    print(f"wrote {os.path.abspath(OUT)} — {len(g['nodes'])} nodes, colour by theme")

if __name__ == "__main__":
    main()
