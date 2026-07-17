#!/usr/bin/env python3
"""Amazon SNAP co-purchase network -> a visualizable annotated graph.json.

Samples an induced subgraph over product ids [0, N) from amazon0601.txt, parses
group + category hierarchy from amazon-meta.txt, and writes web/amazon.graph.json:

  nodes: id, title, group (Book/Music/DVD/Video/…), cat1..cat3 (category levels),
         salesrank; color_map/themes over `group`; meta.ranks declares the
         group→cat1→cat2→cat3 hierarchy so the FunNet layer ladder works here too.
  edges: co-purchase links; weight = category overlap (Jaccard of the two
         products' category tokens), so topical similarity reads as a thicker edge.

The source data is large & local (gitignored); the output graph.json too.
"""
import json, os, re, sys

SNAP = "/Users/eprifti/Research/workspace_r/integranet/data/snap/amazon"
EDGES = os.path.join(SNAP, "amazon0601.txt")
META = os.path.join(SNAP, "amazon-meta.txt")
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "web", "amazon.graph.json")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000   # sample ids [0, N)

# group colours (coarse annotation, like phylum)
GROUP_COLORS = {
    "Book": "#2f6df6", "Music": "#3fb950", "DVD": "#db6d28", "Video": "#a371f7",
    "Toy": "#e3b341", "Software": "#e5534b", "Video Games": "#39c5cf",
    "Baby Product": "#db61a2", "CE": "#7ee787", "Sports": "#ffa657",
}

def parse_meta(ids):
    """Stream amazon-meta.txt; return {id: {group,title,cats,salesrank}} for ids."""
    want = set(ids); info = {}; cur = None; rec = None
    with open(META, encoding="latin-1") as f:
        for line in f:
            if line.startswith("Id:"):
                if rec is not None and cur in want:
                    info[cur] = rec
                cur = int(line.split(":")[1]); rec = {"cats": []}
                if cur >= max(want, default=0):   # stop once past the sampled range
                    if len(info) >= len(want): break
            elif rec is not None:
                s = line.strip()
                if s.startswith("title:"):
                    rec["title"] = s[6:].strip()
                elif s.startswith("group:"):
                    rec["group"] = s[6:].strip()
                elif s.startswith("salesrank:"):
                    try: rec["salesrank"] = int(s.split(":")[1])
                    except Exception: pass
                elif s.startswith("|"):
                    # |Books[283155]|Subjects[1000]|Religion...|  -> ["Books","Subjects",...]
                    toks = [re.sub(r"\[\d+\]$", "", t).strip() for t in s.split("|") if t.strip()]
                    if toks: rec["cats"].append(toks)
        if rec is not None and cur in want:
            info[cur] = rec
    return info

def main():
    # 1) induced subgraph on ids [0, N)
    adj = {}; edges = []
    with open(EDGES) as f:
        for line in f:
            if line.startswith("#"): continue
            a, b = line.split()
            a, b = int(a), int(b)
            if a < N and b < N and a != b:
                edges.append((a, b)); adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    nodes_with_deg = [i for i in range(N) if adj.get(i)]
    print(f"induced subgraph: {len(nodes_with_deg)} connected nodes, {len(edges)} edges", file=sys.stderr)

    # 2) annotation
    info = parse_meta(nodes_with_deg)
    print(f"parsed metadata for {len(info)} products", file=sys.stderr)

    def cats_of(i):
        rec = info.get(i, {})
        paths = rec.get("cats", [])
        # pick the deepest category path; expose the first 3 non-root levels
        best = max(paths, key=len) if paths else []
        # drop the leading "Books"/"Music"… root token (== group) if present
        levels = best[1:4] if len(best) > 1 else []
        return levels

    def catset(i):
        return {t for p in info.get(i, {}).get("cats", []) for t in p}

    keep = [i for i in nodes_with_deg if info.get(i, {}).get("group")]
    keepset = set(keep)
    nodes = []
    for i in keep:
        rec = info[i]; lv = cats_of(i)
        nd = {"id": f"p{i}", "label": (rec.get("title") or f"product {i}")[:60],
              "group": rec.get("group", "NA")}
        for k in range(3):
            nd[f"cat{k+1}"] = lv[k] if k < len(lv) else "NA"
        if "salesrank" in rec and rec["salesrank"] > 0:
            nd["salesrank"] = rec["salesrank"]
        nodes.append(nd)

    ej = []
    for a, b in edges:
        if a in keepset and b in keepset:
            ca, cb = catset(a), catset(b)
            jac = len(ca & cb) / len(ca | cb) if (ca or cb) else 0.0
            ej.append({"source": f"p{a}", "target": f"p{b}", "weight": round(jac, 4)})
    print(f"kept {len(nodes)} annotated nodes, {len(ej)} edges", file=sys.stderr)

    groups = sorted({n["group"] for n in nodes})
    palette = ["#2f6df6","#3fb950","#db6d28","#a371f7","#e3b341","#e5534b","#39c5cf","#db61a2","#7ee787","#ffa657"]
    color_map = {g: GROUP_COLORS.get(g, palette[i % len(palette)]) for i, g in enumerate(groups)}

    g = {
        "meta": {"title": f"Amazon co-purchase — {len(nodes)} products",
                 "source": "SNAP amazon0601 + amazon-meta", "directed": False,
                 "color_by": "group", "annotation_source": "amazon-meta",
                 # aggregation ladder, base (full network) first, then fine → coarse
                 "ranks": ["product", "cat3", "cat2", "cat1", "group"]},
        "color_groups": groups, "color_map": color_map,
        "nodes": nodes, "edges": ej,
    }
    with open(OUT, "w") as f:
        json.dump(g, f)
    print(f"wrote {os.path.abspath(OUT)}", file=sys.stderr)

if __name__ == "__main__":
    main()
