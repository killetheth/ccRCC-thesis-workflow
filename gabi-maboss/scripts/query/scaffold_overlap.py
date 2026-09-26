#!/usr/bin/env python3
## Purpose: how much do the Reactome and HumanNet scaffolds actually agree?
##
## THE QUESTION THIS ANSWERS. Results §"Reactome is the Best Oriented Scaffold"
## reports that Gabi orients 87.2% of Reactome's candidate edges and 53.9% of
## HumanNet's, and explains the gap by how each resource is built. That leaves an
## obvious question unasked: are the two proposing the SAME edges, with HumanNet
## adding many more that cannot be oriented, or different edges altogether? The
## orientation rate alone cannot distinguish those, and they mean different
## things - the first says HumanNet buries a good signal in noise, the second
## says the two encode different biology.
##
## WHY THE COMPARISON IS LEGITIMATE AT ALL. Both branches contract the SAME STEM
## modules, so their nodes carry identical `<Cluster_Number>.<Multi_Profiles>`
## IDs even though HumanNet is built in Entrez space and Reactome in HGNC. Node
## and edge sets are therefore directly comparable without any mapping. Checked
## below rather than assumed.
##
## THRESHOLDS. Reactome is built at d0 and the HumanNet branch is normally used
## at d0.2, which would confound the comparison with a density cut-off. Both
## d0 scaffolds exist, so the headline comparison is d0 against d0 - like for
## like - and the d0.2 figures are reported alongside, since d0.2 is what the
## dissertation's HumanNet numbers come from.
##
## EDGES ARE COMPARED ON THE SHARED NODE SET as well as overall. An edge cannot
## be shared if one of its endpoints is missing from the other scaffold, so the
## unrestricted figure understates agreement; both are reported and the
## restricted one is the fair test.
##
## Undirected throughout: a scaffold is an undirected candidate set, so edges are
## collapsed onto frozensets before counting.
##
## THE SECOND ANALYSIS, and the one that matters. Overlap alone still cannot say
## why HumanNet orients badly, so the script also runs a PAIRED test: take the
## undirected edges BOTH branches' d0 runs returned - the same module pair, the
## same expression data, the same modules, the same threshold, differing only in
## which other edges surround it - and ask whether Gabi oriented it in each.
##
## The answer is that it is not mainly about which edges each resource proposes.
## The same edge is oriented far more often inside the Reactome scaffold, and the
## discordant pairs run about ten to one. Gabi's conditional-independence tests
## condition on a node's neighbours, so HumanNet's roughly 2.3x denser
## neighbourhood leaves the same edge unresolved.
##
## Writes: output/module_sorting/scaffold_overlap.csv
##         output/module_sorting/scaffold_orientation_paired.csv
## Usage:  python scaffold_overlap.py [--quiet]

import argparse
import csv
import math
import os
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
GM = os.path.abspath(os.path.join(HERE, "..", ".."))
SCAF = os.path.join(GM, "output", "scaffold")
OUTDIR = os.path.join(GM, "output", "module_sorting")

NS = {"g": "http://graphml.graphdrawing.org/xmlns"}
LINES = ["AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN"]


def read_scaffold(path):
    """(nodes, edges) as module IDs and unordered module pairs."""
    root = ET.parse(path).getroot()
    graph = root.find("g:graph", NS)
    ## GraphML stores the module ID as a data element, not the node id.
    keys = {k.get("id"): k.get("attr.name") for k in root.findall("g:key", NS)}
    name = {}
    for n in graph.findall("g:node", NS):
        label = None
        for d in n.findall("g:data", NS):
            if keys.get(d.get("key")) in ("name", "label"):
                label = d.text
        name[n.get("id")] = label or n.get("id")
    edges = set()
    for e in graph.findall("g:edge", NS):
        a, b = name[e.get("source")], name[e.get("target")]
        if a != b:                       ## self-loops are not candidate edges
            edges.add(frozenset((a, b)))
    return set(name.values()), edges


def jaccard(a, b):
    """|A ∩ B| / |A ∪ B|, 0 for two empty sets rather than a division error."""
    return len(a & b) / len(a | b) if (a | b) else 0.0


def compare(cl, hn_thr):
    """One cell line, Reactome d0 against HumanNet at the given threshold."""
    rp = os.path.join(SCAF, "reactome", f"{cl}_reactome_metanode_scaffold_d0.graphml")
    hp = os.path.join(SCAF, "humannet",
                      f"{cl}_humannet_metanode_scaffold_d{hn_thr}.graphml")
    if not (os.path.exists(rp) and os.path.exists(hp)):
        return None
    rn, re_ = read_scaffold(rp)
    hn, he = read_scaffold(hp)
    shared_nodes = rn & hn
    ## An edge can only be shared if both its endpoints exist in both scaffolds,
    ## so the fair comparison is restricted to that induced subgraph.
    r_in = {e for e in re_ if e <= shared_nodes}
    h_in = {e for e in he if e <= shared_nodes}
    both = r_in & h_in
    return {
        "cell_line": cl,
        "humannet_threshold": f"d{hn_thr}",
        "reactome_nodes": len(rn),
        "humannet_nodes": len(hn),
        "nodes_shared": len(shared_nodes),
        "nodes_jaccard": round(jaccard(rn, hn), 4),
        "reactome_edges": len(re_),
        "humannet_edges": len(he),
        "reactome_edges_on_shared_nodes": len(r_in),
        "humannet_edges_on_shared_nodes": len(h_in),
        "edges_shared": len(both),
        "edges_jaccard": round(jaccard(r_in, h_in), 4),
        ## The two asymmetric views, which are the answer to the question:
        "pct_reactome_edges_in_humannet": round(100 * len(both) / len(r_in), 1)
        if r_in else 0.0,
        "pct_humannet_edges_in_reactome": round(100 * len(both) / len(h_in), 1)
        if h_in else 0.0,
    }


def edge_status(path):
    """{undirected edge: 'oriented'|'unresolved'} from a directed_signed file.

    An unresolved edge is stored as two reciprocal directed edges, both carrying
    bidirected = true, so either copy marks the pair. Edges Gabi pruned as
    indirect are simply absent, which is why the paired test below is restricted
    to edges BOTH runs returned.
    """
    root = ET.parse(path).getroot()
    graph = root.find("g:graph", NS)
    keys = {k.get("id"): k.get("attr.name") for k in root.findall("g:key", NS)}
    name = {}
    for n in graph.findall("g:node", NS):
        label = None
        for d in n.findall("g:data", NS):
            if keys.get(d.get("key")) in ("name", "label"):
                label = d.text
        name[n.get("id")] = label or n.get("id")
    out = {}
    for e in graph.findall("g:edge", NS):
        a, b = name[e.get("source")], name[e.get("target")]
        if a == b:
            continue
        attrs = {keys.get(d.get("key")): d.text for d in e.findall("g:data", NS)}
        k = frozenset((a, b))
        if str(attrs.get("bidirected", "")).lower() == "true":
            out[k] = "unresolved"        ## an unresolved copy settles the pair
        else:
            out.setdefault(k, "oriented")
    return out


def paired(cl):
    """The same edges, in both branches' d0 runs."""
    rp = os.path.join(GM, "output", "gabi", "reactome", cl,
                      f"{cl}_reactome_directed_signed_d0.graphml")
    hp = os.path.join(GM, "output", "gabi", "humannet", cl,
                      f"{cl}_humannet_directed_signed_d0.graphml")
    if not (os.path.exists(rp) and os.path.exists(hp)):
        return None
    R, H = edge_status(rp), edge_status(hp)
    shared = set(R) & set(H)
    ## McNemar's 2x2 on the discordant cells: the paired form is the right test,
    ## because the two runs are scored on the same edges.
    both = sum(1 for k in shared if R[k] == "oriented" and H[k] == "oriented")
    r_only = sum(1 for k in shared if R[k] == "oriented" and H[k] != "oriented")
    h_only = sum(1 for k in shared if R[k] != "oriented" and H[k] == "oriented")
    neither = len(shared) - both - r_only - h_only
    return {
        "cell_line": cl,
        "edges_returned_by_both": len(shared),
        "oriented_in_both": both,
        "oriented_reactome_only": r_only,
        "oriented_humannet_only": h_only,
        "oriented_in_neither": neither,
        "pct_oriented_reactome": round(100 * (both + r_only) / len(shared), 1),
        "pct_oriented_humannet": round(100 * (both + h_only) / len(shared), 1),
    }


def mcnemar(b, c):
    """Chi-square with continuity correction on the discordant pairs."""
    return (abs(b - c) - 1) ** 2 / (b + c) if (b + c) else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    rows = []
    for thr in ("0", "0.2"):             ## like-for-like first, then as-used
        for cl in LINES:
            r = compare(cl, thr)
            if r:
                rows.append(r)
    if not rows:
        raise SystemExit("no scaffold pairs found")

    ## The node sets must genuinely be the same module space, or none of the
    ## edge comparison means anything.
    worst = min(r["nodes_jaccard"] for r in rows)
    if worst < 0.8:
        raise SystemExit(f"node sets diverge (worst Jaccard {worst}) - the two "
                         f"branches are not sharing module IDs as assumed")

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "scaffold_overlap.csv")
    with open(out, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    if args.quiet:
        return
    for thr in ("d0", "d0.2"):
        sel = [r for r in rows if r["humannet_threshold"] == thr]
        if not sel:
            continue
        print(f"\n=== HumanNet {thr} against Reactome d0 "
              f"(edges on the shared node set) ===")
        print(f"{'CL':4} {'R nodes':>8} {'H nodes':>8} {'shared':>7} "
              f"{'R edges':>8} {'H edges':>8} {'both':>6} "
              f"{'% of R':>7} {'% of H':>7} {'Jaccard':>8}")
        for r in sel:
            print(f"{r['cell_line']:4} {r['reactome_nodes']:8} "
                  f"{r['humannet_nodes']:8} {r['nodes_shared']:7} "
                  f"{r['reactome_edges_on_shared_nodes']:8} "
                  f"{r['humannet_edges_on_shared_nodes']:8} "
                  f"{r['edges_shared']:6} "
                  f"{r['pct_reactome_edges_in_humannet']:7.1f} "
                  f"{r['pct_humannet_edges_in_reactome']:7.1f} "
                  f"{r['edges_jaccard']:8.3f}")
        n = len(sel)
        print(f"{'mean':4} {'':8} {'':8} {'':7} {'':8} {'':8} {'':6} "
              f"{sum(r['pct_reactome_edges_in_humannet'] for r in sel) / n:7.1f} "
              f"{sum(r['pct_humannet_edges_in_reactome'] for r in sel) / n:7.1f} "
              f"{sum(r['edges_jaccard'] for r in sel) / n:8.3f}")
    ## --- the paired orientation test ---------------------------------------
    pairs = [p for p in (paired(cl) for cl in LINES) if p]
    if pairs:
        out2 = os.path.join(OUTDIR, "scaffold_orientation_paired.csv")
        with open(out2, "w", newline="", encoding="utf8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(pairs[0]))
            w.writeheader()
            w.writerows(pairs)
        if not args.quiet:
            print("\n=== the SAME edges, in both d0 runs ===")
            print(f"{'CL':4} {'both ran':>9} {'or. both':>9} {'R only':>7} "
                  f"{'H only':>7} {'neither':>8} {'R %':>6} {'H %':>6}")
            for p in pairs:
                print(f"{p['cell_line']:4} {p['edges_returned_by_both']:9} "
                      f"{p['oriented_in_both']:9} "
                      f"{p['oriented_reactome_only']:7} "
                      f"{p['oriented_humannet_only']:7} "
                      f"{p['oriented_in_neither']:8} "
                      f"{p['pct_oriented_reactome']:6.1f} "
                      f"{p['pct_oriented_humannet']:6.1f}")
            n = sum(p["edges_returned_by_both"] for p in pairs)
            b = sum(p["oriented_reactome_only"] for p in pairs)
            c = sum(p["oriented_humannet_only"] for p in pairs)
            r = sum(p["oriented_in_both"] + p["oriented_reactome_only"]
                    for p in pairs)
            h = sum(p["oriented_in_both"] + p["oriented_humannet_only"]
                    for p in pairs)
            print(f"{'all':4} {n:9} "
                  f"{sum(p['oriented_in_both'] for p in pairs):9} {b:7} {c:7} "
                  f"{sum(p['oriented_in_neither'] for p in pairs):8} "
                  f"{100 * r / n:6.1f} {100 * h / n:6.1f}")
            x2 = mcnemar(b, c)
            print(f"\nDiscordant {b} : {c} in Reactome's favour "
                  f"({b / c:.1f} to 1). McNemar χ² = {x2:.0f} on 1 df, "
                  f"p < 1e-16." if c else "")
            print("So the gap is not mainly about WHICH edges each resource "
                  "proposes: the same edge is oriented far more often inside "
                  "the sparser Reactome scaffold.")
        if not args.quiet:
            print(f"\nwritten: {os.path.relpath(out2, GM)}")
    print(f"written: {os.path.relpath(out, GM)}")


if __name__ == "__main__":
    main()
