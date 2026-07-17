#!/usr/bin/env python3
"""Flickr SNAP image graph -> an annotated graph.json (tag as the annotation).

Nodes = photos, edges = shared-context relationships (flickrEdges.txt). The MIR
node-feature file gives each photo a bag of tags over a 20000-tag vocabulary; we
take each photo's **dominant tag** ("sky", "nikon", …) as its categorical
annotation. So the same analysis (ATC propagation, FunNet tag-proximity layer,
guilds, provenance) runs on a social/image graph — a very different domain.

  nodes: id, tag; color_map/color_by over `tag`; meta.ranks = ["photo","tag"].
  edges: co-context links (both endpoints tagged), ecorr = 1 (unweighted).

Source data is large & local (gitignored); output web/flickr.graph.json too.
"""
import json, os, sys, tarfile, io as _io

FDB = "/Users/eprifti/Research/workspace_r/integranet_serveur/datal/FlickrDB"
EDGES = os.path.join(FDB, "flickrEdges.txt")
FEATS_TAR = os.path.join(FDB, "nodeFeaturesFlickr.tar.gz")
FEATS_MEMBER = "nodeFeaturesFlickr/nodeFeaturesMIR.txt"
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "web", "flickr.graph.json")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500

def load_photo_tags():
    """photo_id -> dominant tag name, from the MIR feature file."""
    with tarfile.open(FEATS_TAR) as tf:
        f = tf.extractfile(FEATS_MEMBER)
        txt = _io.TextIOWrapper(f, encoding="latin-1")
        header = txt.readline().split()
        n_users, _, _ = int(header[0]), header[1], header[2]
        for _ in range(n_users):            # skip user id-map
            txt.readline()
        tag_by_pos = {}                     # position -> tag name
        line = txt.readline()
        pos0 = None
        while line:
            parts = line.split()
            # tag-vocabulary lines are `index tagname` (2 tokens, small integer)
            if len(parts) == 2 and parts[0].isdigit():
                idx = int(parts[0])
                if pos0 is None:
                    pos0 = idx
                tag_by_pos[idx - pos0] = parts[1]
                line = txt.readline()
            else:
                break                       # reached the matrix
        # matrix rows: `photo_id owner_id <countstring over the tag vocabulary>`
        photo_tag = {}
        while line:
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                pid = int(parts[0]); vec = parts[2]
                best_p, best_c = -1, 0
                for p, ch in enumerate(vec):
                    if ch != "." and ch.isdigit():
                        c = int(ch)
                        if c > best_c:
                            best_c, best_p = c, p
                if best_p >= 0 and best_p in tag_by_pos:
                    photo_tag[pid] = tag_by_pos[best_p]
            line = txt.readline()
        return photo_tag

def main():
    photo_tag = load_photo_tags()
    print(f"{len(photo_tag)} photos with a dominant tag "
          f"({len(set(photo_tag.values()))} distinct tags)", file=sys.stderr)
    tagged = set(photo_tag)

    # induced subgraph on tagged photos, grown by BFS from a well-connected seed
    adj = {}
    with open(EDGES) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a, b = line.split()
            a, b = int(a), int(b)
            if a in tagged and b in tagged and a != b:
                adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    if not adj:
        print("no tagged edges found", file=sys.stderr); sys.exit(1)
    seed = max(adj, key=lambda k: len(adj[k]))
    keep, frontier = {seed}, [seed]
    while frontier and len(keep) < N:
        nxt = []
        for u in frontier:
            for v in adj[u]:
                if v not in keep:
                    keep.add(v); nxt.append(v)
                    if len(keep) >= N:
                        break
            if len(keep) >= N:
                break
        frontier = nxt
    nodes = [{"id": f"ph{p}", "tag": photo_tag[p]} for p in keep]
    edges = []
    for u in keep:
        for v in adj[u]:
            if v in keep and u < v:
                edges.append({"source": f"ph{u}", "target": f"ph{v}", "ecorr": 1.0})
    print(f"subgraph: {len(nodes)} photos, {len(edges)} edges", file=sys.stderr)

    tags = sorted({n["tag"] for n in nodes})
    palette = ["#2f6df6","#3fb950","#db6d28","#a371f7","#e3b341","#e5534b","#39c5cf","#db61a2",
               "#7ee787","#ffa657","#79c0ff","#d2a8ff","#f78166","#56d364","#e3b341"]
    color_map = {t: palette[i % len(palette)] for i, t in enumerate(tags)}
    g = {
        "meta": {"title": f"Flickr image graph — {len(nodes)} photos",
                 "source": "SNAP flickrEdges + MIR tag features", "directed": False,
                 "color_by": "tag", "ranks": ["photo", "tag"]},
        "color_groups": tags, "color_map": color_map,
        "nodes": nodes, "edges": edges,
    }
    with open(OUT, "w") as f:
        json.dump(g, f)
    print(f"wrote {os.path.abspath(OUT)}", file=sys.stderr)

if __name__ == "__main__":
    main()
