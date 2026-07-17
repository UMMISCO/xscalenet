#!/usr/bin/env python3
"""Regenerate web/projects.json — the network index the app's selector reads.

Scans project/*/project.json and emits the app-facing subset to web/projects.json.
Run after adding a project or (re)building a graph.
"""
import json, os, glob

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "web", "projects.json")
WEB = os.path.join(HERE, "..", "web")
FIELDS = ["id", "title", "domain", "graph", "color_by", "note", "client_data"]

def main():
    projects = []
    for pj in sorted(glob.glob(os.path.join(HERE, "*", "project.json"))):
        m = json.load(open(pj))
        entry = {k: m[k] for k in FIELDS if k in m}
        # mark whether the generated graph is actually present (servable)
        entry["available"] = os.path.exists(os.path.join(WEB, m.get("graph", "")))
        projects.append(entry)
    # available first, then by title
    projects.sort(key=lambda p: (not p.get("available"), p.get("title", "")))
    json.dump({"projects": projects}, open(OUT, "w"), indent=1)
    print(f"wrote {os.path.abspath(OUT)} — {len(projects)} networks "
          f"({sum(p['available'] for p in projects)} available)")

if __name__ == "__main__":
    main()
