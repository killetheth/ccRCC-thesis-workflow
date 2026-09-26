#!/usr/bin/env python3
"""Rebuild the scaffold-comparison tables of the dissertation (Tables 4, B1, B2).

Every number in Results section "Reactome Was the Only Usable Scaffold at Module
Level" is derived here.  **One convention, applied throughout: what went into
Gabi is read from the scaffold build log, what came out is read from the signed
GraphML.**  The two differ because Gabi prunes indirect edges, and reporting
each side from its own source keeps that visible rather than hiding it.

Per dataset (Table B1):
  Modules      nodes in the directed network (isolated meta-nodes are absent)
  Edges in     candidate edges from logs/scaffold/<branch>/<CL>_scaffold*.log,
               before any pruning - the input side
  Pruned       Edges in minus the relationships Gabi returned; these are the
               edges its triangle-CMI test judged indirect
  Unresolved   edges returned without a direction.  An unoriented edge is stored
               as two reciprocal directed edges, so ordered pairs are collapsed
               onto {u,v} first or every one of them counts twice
  Oriented     returned relationships minus Unresolved - the output side
  Oriented %   100 * Oriented / Edges in, so the denominator is what Gabi was
               given rather than what survived its own pruning

Per branch (Table 4): the same quantities averaged over the eight datasets, plus
the highest- and lowest-orienting dataset of each branch.  Relay shares are also
computed here for reference, but note they use the strict definition (in-degree
> 0 AND out-degree > 0) rather than 5-sort_modules.R's, which counts an isolated
node as a Relay; Table B2 of the dissertation reports role composition from that
script instead, so do not cross-quote the two.

Run from gabi-maboss/scripts/query/ ; writes CSVs beside the other query output.
"""

import argparse                                  # command-line flags
import re                                        # parsing the scaffold logs
import csv                                       # table output
import os                                        # path handling
import statistics                                # branch means
import xml.etree.ElementTree as ET               # GraphML is plain XML

HERE = os.path.dirname(os.path.abspath(__file__))            # pin to this script
GABI = os.path.normpath(os.path.join(HERE, "..", "..", "output", "gabi"))
OUT  = os.path.normpath(os.path.join(HERE, "..", "..", "output", "module_sorting"))
LOGS = os.path.normpath(os.path.join(HERE, "..", "..", "logs", "scaffold"))
NS   = {"g": "http://graphml.graphdrawing.org/xmlns"}        # GraphML namespace
CELL_LINES = ["AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN"]
BRANCHES   = {"humannet": "HumanNet", "reactome": "Reactome", "huri": "HuRI"}


def read_scaffold_log(branch, cl):
    """Candidate edge count as built, before Gabi saw it. Returns (modules, edges)."""
    path = os.path.join(LOGS, branch, f"{cl}_scaffold_d0.log")      # d0 rerun, if present
    if not os.path.exists(path):
        path = os.path.join(LOGS, branch, f"{cl}_scaffold.log")     # single-threshold runs
    text = open(path).read()
    m = re.search(r"genes used: (\d+) \| meta-nodes: (\d+) \| candidate edges: (\d+)", text)
    return int(m.group(2)), int(m.group(3))                          # meta-nodes, candidate edges


def read_network(path):
    """Return (node ids, directed edges) from a signed GraphML file."""
    root = ET.parse(path).getroot()                          # parse the XML
    graph = root.find("g:graph", NS)                         # single <graph>
    nodes = [n.get("id") for n in graph.findall("g:node", NS)]
    edges = [(e.get("source"), e.get("target"))              # ordered pairs
             for e in graph.findall("g:edge", NS)]
    return nodes, edges


def summarise(nodes, edges):
    """Collapse directed edges to undirected relationships and score the roles."""
    directed = set(edges)                                    # de-duplicate first
    undirected = {frozenset(e) for e in directed if e[0] != e[1]}   # drop loops
    # An unresolved edge is present in both directions, so both ordered pairs
    # exist; a resolved one appears only once.
    unresolved = {p for p in undirected
                  if tuple(p) in directed and tuple(p)[::-1] in directed
                  or (len(p) == 2 and all((a, b) in directed
                                          for a, b in ((tuple(p)[0], tuple(p)[1]),
                                                       (tuple(p)[1], tuple(p)[0]))))}
    resolved_only = [e for e in directed                     # backbone edge list
                     if frozenset(e) not in unresolved and e[0] != e[1]]

    def relay_share(edge_list):
        """Share of nodes that both receive and send — i.e. relay nodes."""
        indeg  = {n: 0 for n in nodes}                       # start every node at 0
        outdeg = {n: 0 for n in nodes}
        for u, v in edge_list:                               # count each direction
            outdeg[u] = outdeg.get(u, 0) + 1
            indeg[v]  = indeg.get(v, 0) + 1
        relays = sum(1 for n in nodes if indeg[n] > 0 and outdeg[n] > 0)
        return 100.0 * relays / len(nodes) if nodes else float("nan")

    return {
        "modules":        len(nodes),
        "returned":       len(undirected),
        "unresolved":     len(unresolved),
        "oriented":       len(undirected) - len(unresolved),
        "relay_all":      relay_share(list(directed)),
        "relay_backbone": relay_share(resolved_only),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=OUT, help="where to write the CSVs")
    ap.add_argument("--quiet", action="store_true", help="suppress the printed tables")
    args = ap.parse_args()

    per_dataset, per_branch = [], []
    for branch, label in BRANCHES.items():                   # one block per scaffold
        rows = []
        for cl in CELL_LINES:
            path = os.path.join(GABI, branch, cl,
                                f"{cl}_{branch}_directed_signed_d0.graphml")
            if not os.path.exists(path):                     # Reactome omits the tag
                path = os.path.join(GABI, branch, cl,
                                    f"{cl}_directed_signed_d0.graphml")
            if not os.path.exists(path):                     # report, do not guess
                print(f"  ! missing: {path}")
                continue
            s = summarise(*read_network(path))               # the output side
            _, edges_in = read_scaffold_log(branch, cl)      # the input side
            s["edges_in"] = edges_in
            s["pruned"]   = edges_in - s["returned"]         # what the CMI test removed
            s["oriented_pct"] = 100.0 * s["oriented"] / edges_in
            s.update(scaffold=label, cell_line=cl)
            rows.append(s)
            per_dataset.append(s)
        if not rows:
            continue
        mean = lambda k: statistics.mean(r[k] for r in rows) # branch mean helper
        per_branch.append({
            "scaffold":     label,
            "modules":      round(mean("modules")),
            "edges_in":     round(mean("edges_in")),
            "pruned":       round(mean("pruned")),
            "returned":     round(mean("returned")),
            "unresolved":   round(mean("unresolved")),
            "oriented":     round(mean("oriented")),
            # Table 4's caption says "means across the eight datasets", so the
            # branch rate is the mean of the eight per-dataset rates. The ratio of
            # the two branch totals is also kept, since it differs whenever the
            # datasets are of unequal size, and quoting the wrong one is easy.
            "oriented_pct":       round(mean("oriented_pct"), 1),
            "oriented_pct_pooled": round(100.0 * mean("oriented") / mean("edges_in"), 1),
            "relay_all":      round(mean("relay_all"), 1),
            "relay_backbone": round(mean("relay_backbone"), 1),
            "highest": max(rows, key=lambda r: r["oriented_pct"])["cell_line"],
            "highest_pct": round(max(r["oriented_pct"] for r in rows), 1),
            "lowest":  min(rows, key=lambda r: r["oriented_pct"])["cell_line"],
            "lowest_pct": round(min(r["oriented_pct"] for r in rows), 1),
        })

    os.makedirs(args.outdir, exist_ok=True)                  # ensure the folder
    p1 = os.path.join(args.outdir, "scaffold_comparison_per_dataset.csv")
    p2 = os.path.join(args.outdir, "scaffold_comparison_per_branch.csv")
    with open(p1, "w", newline="") as f:                     # Table B1's source
        w = csv.DictWriter(f, fieldnames=["scaffold", "cell_line", "modules",
                                          "edges_in", "pruned", "returned",
                                          "unresolved", "oriented", "oriented_pct",
                                          "relay_all", "relay_backbone"])
        w.writeheader()
        for r in per_dataset:
            w.writerow({k: (round(v, 1) if isinstance(v, float) else v)
                        for k, v in r.items()})
    with open(p2, "w", newline="") as f:                     # Tables 4 and B2
        w = csv.DictWriter(f, fieldnames=list(per_branch[0].keys()))
        w.writeheader()
        w.writerows(per_branch)

    if not args.quiet:
        print(f"{'scaffold':<10}{'CL':<4}{'mods':>6}{'in':>7}{'pruned':>8}"
              f"{'unres':>7}{'oriented':>10}{'orient%':>9}")
        for r in per_dataset:
            print(f"{r['scaffold']:<10}{r['cell_line']:<4}{r['modules']:>6}"
                  f"{r['edges_in']:>7}{r['pruned']:>8}{r['unresolved']:>7}"
                  f"{r['oriented']:>10}{r['oriented_pct']:>9.1f}")
        print()
        for r in per_branch:
            print(f"{r['scaffold']:<10} modules={r['modules']:<5}in={r['edges_in']:<6}"
                  f"pruned={r['pruned']:<5}oriented={r['oriented']:<6}"
                  f"{r['oriented_pct']}% (pooled {r['oriented_pct_pooled']}%)"
                  f"   hi={r['highest']} ({r['highest_pct']})"
                  f"  lo={r['lowest']} ({r['lowest_pct']})")
    print(f"\nwrote {p1}\n      {p2}")


if __name__ == "__main__":
    main()
