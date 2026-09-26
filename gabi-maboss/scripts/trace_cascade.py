#!/usr/bin/env python3
## Purpose: Explain a compare_maboss.py delta table by the network underneath it -
## for each module that moved, how many hops it sits from the perturbed module and
## whether the sign of the shortest path predicts the direction it moved. A delta
## table says which modules changed; this says how the change reached them, which
## is what separates a direct target from the far end of a cascade, and a real
## cascade from a module the perturbation cannot reach at all.
##
## Reads the `.bnet` MaBoSS actually ran (not the SIF), so the topology traced is
## the topology simulated. Rules from 7-sif_to_bnet.py are flat - an OR-group of
## activators, optionally `& !(...)` an OR-group of inhibitors, no nesting - so a
## regulator's sign is simply whether it sits inside the negated group.
##
## Reading the output:
##   hop      shortest directed distance from the perturbed module (0 = the target
##            itself, 1 = a direct target, higher = reached through intermediates).
##   sign?    `ok` if the product of edge signs along that shortest path predicts
##            the direction the module actually moved, `FLIP` if it doesn't. FLIP
##            is not an error - it means the shortest path is not the route that
##            drove the module, i.e. a longer or competing path dominates, so a
##            cluster of FLIPs marks multi-path control worth looking at.
##   reachable  modules the perturbation can reach at all. A small number here says
##            the perturbation is confined regardless of what the delta table shows.
##
## --path MODULE prints the explicit signed chain from the perturbed module to that
## module (`A --> B --| C`, `-->` activation, `--|` inhibition), for writing a
## cascade up once the hop table has pointed at it.
##
## Each DELTA_CSV is a `<stem>_<tag>_vs_wt_delta.csv` from compare_maboss.py. The
## model is found automatically: the `.bnet` in gabi-maboss/output/for_maboss/ whose
## stem the delta filename starts with, with that model's `<stem>_name_map.tsv`
## beside it (--model-dir overrides the folder). The perturbed module and whether it
## was forced ON or OFF are read from the table's own is_target row, so the tag
## never has to be parsed back.
##
## Files (--csv [OUTDIR]): <stem>_<tag>_vs_wt_cascade.csv - the hop table, with the
##   delta columns carried over so it stands alone. Writes to gabi-maboss/output/maboss/
##   when given bare; nothing is written unless the flag is present.
##
## Pure standard library - no MaBoSS engine and no `maboss` package, so unlike
## 8-run_maboss.py and compare_maboss.py this runs outside the `maboss` conda env.
## That is why it re-declares the name-map reader rather than importing
## maboss_common.py, whose `import maboss` would pull the env back in for what is a
## text-and-graph analysis.
##
## Usage:
##   python trace_cascade.py DELTA_CSV [DELTA_CSV ...] [--min-delta X]
##          [--path MODULE ...] [--model-dir DIR] [--csv [OUTDIR]]

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict, deque

## Both resolved from this script's own location (gabi-maboss/scripts/) rather than
## the caller's cwd - the Python counterpart of the R scripts' setwd() pin, so the
## defaults hold whichever folder you run from.
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_DIR = os.path.normpath(os.path.join(_HERE, "..", "output", "for_maboss"))
DEFAULT_OUTDIR = os.path.normpath(os.path.join(_HERE, "..", "output", "maboss"))

## A rule is `(act | act) & !(inh | inh)`, with either group optional and neither
## nested. Everything inside a !( ) group is an inhibitor; everything left over is
## an activator.
INHIB_GROUP = re.compile(r"!\(([^()]*)\)")
TOKEN = re.compile(r"\w+")


def load_name_map(path):
    """safe_name -> original module ID, from a 6-cyto_to_neko.py name-map TSV."""
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["safe_name"]: r["original_name"]
                for r in csv.DictReader(fh, delimiter="\t")}


def load_rules(bnet):
    """regulators[target] = {regulator: +1 activation / -1 inhibition / 0 both},
    in the model's safe-name space."""
    regulators = defaultdict(dict)
    with open(bnet, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("targets"):
                continue
            target, sep, expr = line.partition(",")
            if not sep:
                raise ValueError(f"{bnet}:{lineno}: not a `target, factors` rule")
            target, expr = target.strip(), expr.strip()
            inhibitors = {t for group in INHIB_GROUP.findall(expr)
                          for t in TOKEN.findall(group)}
            rest = INHIB_GROUP.sub(" ", expr)
            ## Every negation in a 7-sif_to_bnet.py rule is a !( ) group; a bare !X
            ## would leave its regulator in `rest` and be silently signed as an
            ## activator, so refuse the file rather than mis-sign the network.
            if "!" in rest:
                raise ValueError(f"{bnet}:{lineno}: negation outside a !( ) group, "
                                 f"which this parser would mis-sign: {expr!r}")
            activators = set(TOKEN.findall(rest))
            for node in activators | inhibitors:
                both = node in activators and node in inhibitors
                regulators[target][node] = 0 if both else (-1 if node in inhibitors else +1)
    return regulators


def forward_edges(regulators):
    """Invert the rule table: regulator -> {target: sign}."""
    forward = defaultdict(dict)
    for target, regs in regulators.items():
        for node, sign in regs.items():
            forward[node][target] = sign
    return forward


def breadth_first(forward, source):
    """Shortest directed distance from source, with the sign product and the
    predecessor along that path. Ties are broken by BFS order, so `sign` describes
    one shortest path, not all of them - which is exactly what a FLIP flags."""
    dist = {source: 0}
    sign = {source: +1}
    prev = {source: None}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        for nxt, edge_sign in forward.get(node, {}).items():
            if nxt not in dist:
                dist[nxt] = dist[node] + 1
                sign[nxt] = sign[node] * edge_sign
                prev[nxt] = (node, edge_sign)
                queue.append(nxt)
    return dist, sign, prev


def find_model(delta_csv, model_dir):
    """The `.bnet` whose stem the delta filename starts with. compare_maboss.py
    names its table `<stem>_<tag>_vs_wt_delta.csv`, so the longest matching stem is
    that model - matching on the stem avoids having to reverse the mutation tag."""
    base = os.path.basename(delta_csv)
    stems = [(os.path.splitext(os.path.basename(p))[0], p)
             for p in glob.glob(os.path.join(model_dir, "*.bnet"))]
    hits = sorted((s for s in stems if base.startswith(s[0])),
                  key=lambda s: -len(s[0]))
    if not hits:
        raise FileNotFoundError(
            f"no .bnet in {model_dir} whose stem starts {base!r} - pass --model-dir")
    return hits[0][1]


def read_delta(path):
    """The delta table, plus the perturbed module and the direction it was forced.
    Direction comes from the target's own treated probability rather than the
    filename tag: compare_maboss.py holds the treated node at 0 or 1 all run, so
    P_ON_treated is 1.0 for an ON treatment and 0.0 for an OFF one."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    targets = [r for r in rows if r["is_target"] == "True"]
    if len(targets) != 1:
        raise ValueError(f"{path}: expected exactly 1 is_target row, found {len(targets)}")
    target = targets[0]
    direction = +1 if float(target["P_ON_treated"]) > 0.5 else -1
    return rows, target["Module_ID"], direction


def chain(prev, nmap, node):
    """The explicit signed chain back to the BFS source, as display strings."""
    steps = []
    while prev.get(node):
        upstream, sign = prev[node]
        steps.append((nmap.get(upstream, upstream),
                      "-->" if sign > 0 else ("--|" if sign < 0 else "--?"),
                      nmap.get(node, node)))
        node = upstream
    steps.reverse()
    return steps


def trace(delta_csv, args):
    model = find_model(delta_csv, args.model_dir)
    stem = os.path.splitext(os.path.basename(model))[0]
    map_path = os.path.join(os.path.dirname(model), stem + "_name_map.tsv")
    nmap = load_name_map(map_path) if os.path.isfile(map_path) else {}
    rev = {v: k for k, v in nmap.items()}

    rows, target_id, direction = read_delta(delta_csv)
    source = rev.get(target_id, target_id)
    forward = forward_edges(load_rules(model))
    dist, sign, prev = breadth_first(forward, source)

    n_nodes = len(nmap) if nmap else len({n for t in forward for n in (t,)} | set(forward))
    print("\n" + "=" * 96)
    print("%s   [%s %s]" % (os.path.basename(delta_csv), target_id,
                            "ON" if direction > 0 else "OFF"))
    print("model: %s   |  reaches %d of %d modules" % (stem, len(dist) - 1, n_nodes))
    print("-" * 96)

    moved = [r for r in rows if abs(float(r["delta_final"])) >= args.min_delta]
    moved.sort(key=lambda r: (dist.get(rev.get(r["Module_ID"], r["Module_ID"]), 10 ** 6),
                              -abs(float(r["delta_final"]))))

    print("%-4s %-11s %-11s %8s %7s  %-5s %s" %
          ("hop", "module", "display", "delta", "t_peak", "sign?", "function"))
    out = []
    for row in moved:
        safe = rev.get(row["Module_ID"], row["Module_ID"])
        hop = dist.get(safe)
        delta = float(row["delta_final"])
        if hop is None:
            verdict, hop_txt = "", "--"
        else:
            predicted = sign[safe] * direction
            verdict = "" if predicted == 0 else ("ok" if (predicted > 0) == (delta > 0) else "FLIP")
            hop_txt = str(hop)
        print("%-4s %-11s %-11s %+8.3f %7s  %-5s %s" %
              (hop_txt, row["Module_ID"], row["Display"], delta, row["t_delta_max"],
               verdict, row["Cluster Function"][:34]))
        out.append({"Module_ID": row["Module_ID"], "Display": row["Display"],
                    "hop": hop_txt, "path_sign_agrees": verdict,
                    "delta_final": row["delta_final"], "delta_max": row["delta_max"],
                    "t_delta_max": row["t_delta_max"],
                    "Cluster Function": row["Cluster Function"],
                    "Module_HGNCs": row["Module_HGNCs"]})

    unreached = sum(1 for r in out if r["hop"] == "--")
    flips = sum(1 for r in out if r["path_sign_agrees"] == "FLIP")
    print("-" * 96)
    print("%d modules moved by >= %g: %d unreachable from the target, %d FLIP"
          % (len(out), args.min_delta, unreached, flips))

    for module in args.path:
        safe = rev.get(module, module)
        if safe not in dist:
            print("\n  %s: not reachable from %s" % (module, target_id))
            continue
        steps = chain(prev, nmap, safe)
        line = [steps[0][0]] if steps else [module]
        for _, arrow, downstream in steps:
            line.append("%s %s" % (arrow, downstream))
        label = next((r for r in rows if r["Module_ID"] == module), None)
        print("\n  %s%s  hop %d, delta %+.3f at t=%s"
              % (module,
                 "" if not label or label["Display"] == module else " (%s)" % label["Display"],
                 dist[safe],
                 float(label["delta_final"]) if label else float("nan"),
                 label["t_delta_max"] if label else "?"))
        print("      " + " ".join(line))

    if args.csv is not None:
        outdir = args.csv or DEFAULT_OUTDIR
        os.makedirs(outdir, exist_ok=True)
        dest = os.path.join(
            outdir,
            os.path.basename(delta_csv).replace("_delta.csv", "_cascade.csv"))
        with open(dest, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(out[0].keys()) if out else
                                    ["Module_ID", "Display", "hop"])
            writer.writeheader()
            writer.writerows(out)
        print("-> %s" % dest)


def main():
    parser = argparse.ArgumentParser(
        description="Trace how a compare_maboss.py perturbation reached each module "
                    "that moved: hops from the target, and whether the path sign "
                    "predicts the direction.")
    parser.add_argument("delta_csv", nargs="+", metavar="DELTA_CSV",
                        help="*_vs_wt_delta.csv from compare_maboss.py")
    parser.add_argument("--min-delta", type=float, default=0.05, metavar="X",
                        help="only trace modules whose |delta_final| is at least X "
                             "(default 0.05; the delta table's own above_noise flag "
                             "is near-useless as a filter - MaBoSS is deterministic "
                             "per seed, so almost everything clears it)")
    parser.add_argument("--path", action="append", default=[], metavar="MODULE",
                        help="print the explicit signed chain to this module ID "
                             "(repeatable)")
    parser.add_argument("--model-dir", default=DEFAULT_MODEL_DIR, metavar="DIR",
                        help="where the .bnet and its name map live "
                             "(default gabi-maboss/output/for_maboss/)")
    parser.add_argument("--csv", nargs="?", const="", default=None, metavar="OUTDIR",
                        help="write the hop table as <stem>_<tag>_vs_wt_cascade.csv "
                             "(bare: gabi-maboss/output/maboss/)")
    args = parser.parse_args()

    failed = 0
    for path in args.delta_csv:
        try:
            trace(path, args)
        except Exception as exc:                      # keep going: files are independent
            print("!! %s: %s" % (os.path.basename(path), exc), file=sys.stderr)
            failed += 1
    print("\n%d table(s), %d failed" % (len(args.delta_csv), failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
