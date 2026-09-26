#!/usr/bin/env python3
## Purpose: Run one or more MaBoSS `.bnet` models and report the result - the
## final-state probability table, the fixed points (optionally written to CSV),
## and an optional node-trajectory graph, under optional module knock-outs/gains.
## Generalises the old UH-only smoke test (`6-run_uh.py`) into a small CLI, so you
## can analyse just the networks you picked out in Cytoscape and converted through
## 6-cyto_to_neko.py -> 7-sif_to_bnet.py.
##
## To simulate a treatment (wild-type vs treated, per-module change in P(ON)),
## use compare_maboss.py instead - it shares this script's model-building code
## (maboss_common.py), so its wild-type arm is built exactly as a plain run here.
##
## Each MODEL is a `.bnet` from 7-sif_to_bnet.py. Its node names are the MaBoSS-safe
## ones 6-cyto_to_neko.py produced (`m6_34`, …), and that network's
## `<stem>_name_map.tsv` translates them back to the original module IDs
## (`6.34`, …) - relabelling every table and graph, and letting --mutate take
## module IDs. The map is found automatically: for each model, the script looks
## for `<stem>_name_map.tsv` beside it, which is where 6-cyto_to_neko.py puts it
## (model and map both default to gabi-maboss/output/for_maboss/, the folder for
## what MaBoSS is fed). So several models can be run in one go and each is
## labelled with its OWN map. --name-map FILE overrides that for every model, for
## the odd case where the map sits elsewhere.
##
## Reading the output:
##   * A state is written as its ON nodes, `--`-separated (`<nil>` = all OFF).
##   * The terminal shows the fixed points (one row each) and only a COUNT of the
##     final states - get_last_states_probtraj returns them as a single row with a
##     column per attractor, which prints as an unreadable wall once there are
##     more than a handful. Use --csv to get that table as <stem>_last_states.csv.
##   * Self-input modules (logic = themselves) hold their initial value forever;
##     loadBNet starts every node OFF, so by default you see the single all-inputs
##     -OFF attractor. --free-inputs starts each such input at 50/50 instead, so
##     the run samples every input combination and reveals the full landscape.
##
## Perturbations:
##   --mutate 18.37=OFF   knock a module out (locked OFF); =ON locks it ON.
##   Repeatable. The node accepts a module ID (with --name-map) or a safe name.
##   Mutated nodes are never freed by --free-inputs. Output files for a mutated run
##   are tagged with the mutation (e.g. <stem>_18_37OFF_fixpoints.csv) so they
##   don't overwrite the wild-type run.
##
## Files (--csv [OUTDIR]):  <stem>_last_states.csv, <stem>_fixpoints.csv
## Graphs (--plot [OUTDIR]): <stem>_node_traj.png - per-module P(ON) over time,
##   the only graph of a plain run. State-level plots (top states over time,
##   final-state pie) were dropped as not relevant here; the state *tables* are
##   unaffected. --fig-width/--fig-height set its size in inches; a long --max-time
##   or a crowded legend usually wants a wider figure than the 9x5 default.
## Both write to gabi-maboss/output/maboss/ when given bare; name a folder to
## override. Nothing is written unless the flag is present - a plain run just prints.
##
## Usage:
##   python 8-run_maboss.py MODEL.bnet [MODEL.bnet ...]
##          [--name-map FILE] [--mutate NODE=ON|OFF ...] [--free-inputs] [--seed N]
##          [--plot [OUTDIR]] [--fig-width IN] [--fig-height IN] [--csv [OUTDIR]]
##          [--max-time T] [--sample-count N] [--threads N]
## Run in the activated `maboss` conda env so the `maboss` package AND the MaBoSS
## engine binary (`MaBoSS`, on the env's PATH) are both found.

import argparse
import os
import sys

from maboss_common import (DEFAULT_OUTDIR, build_sim, find_name_map, load_name_map,
                            mutation_tag, node_traj)


def relabel_state(state, nmap):
    """Translate a MaBoSS state label (`m21_14 -- m6_34`) to module IDs."""
    if not nmap or state == "<nil>":
        return state
    return " -- ".join(nmap.get(tok, tok) for tok in state.split(" -- "))


def relabelled_fptable(result, nmap):
    """The fixed-point table with node columns and the State column as module IDs,
    or None if there are no fixed points."""
    fp = result.get_fptable()
    if fp is None or len(fp) == 0:
        return None
    fp = fp.rename(columns=lambda c: nmap.get(c, c)).copy()
    if "State" in fp:
        fp["State"] = fp["State"].map(lambda s: relabel_state(s, nmap))
    return fp


def save_plots(result, stem, outdir, nmap, figsize):
    """Write the node-trajectory graph for one result.

    Per-module activation is the readout this work is after, so the graph side is
    deliberately just that one plot - the state-level views (top states over time,
    final-state pie) were dropped. The state-level *tables* are untouched: the
    attractor probabilities still print, and --csv still writes them.

    `figsize` is the (width, height) in inches from --fig-width/--fig-height: a
    long run, or one with enough modules to fill the legend, needs a wider figure
    than the default before the lines stop overlapping into a solid band.
    """
    import matplotlib
    matplotlib.use("Agg")                          # headless: write files, no display
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    written = []

    nodes = node_traj(result, nmap)
    ax = nodes.plot(figsize=figsize)
    ax.set_xlabel("time"); ax.set_ylabel("P(node = ON)")
    ax.set_title(f"{stem}: module activation over time")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, title="module")
    p = os.path.join(outdir, f"{stem}_node_traj.png")
    ax.figure.savefig(p, bbox_inches="tight", dpi=150); plt.close(ax.figure); written.append(p)

    return written


def run_one(bnet, args):
    """Load, run, and report one `.bnet` model. Returns True on success."""
    ## This model's own name map, so a multi-model run relabels each network
    ## correctly rather than applying one map to all of them.
    map_path = find_name_map(bnet, args.name_map)
    nmap = load_name_map(map_path) if map_path else {}
    rev  = {orig: safe for safe, orig in nmap.items()}   # module ID -> safe name, for --mutate

    stem = os.path.splitext(os.path.basename(bnet))[0]
    print(f"\n=== {bnet} ===")

    ## Say which map was used - a silently missing one would otherwise show up
    ## only as tables full of safe names, and would make --mutate reject module IDs.
    if map_path:
        how = "given" if args.name_map else "found beside the model"
        print(f"name map ({how}): {map_path}")
    else:
        print("name map: none found - labels stay as the MaBoSS-safe node names")

    sim, mutations, freed = build_sim(bnet, args, args.mutate, rev, nmap, args.seed)
    for disp, state in mutations:
        print(f"mutation: {disp} {state}")
    if args.free_inputs:
        shown = ", ".join(nmap.get(n, n) for n in freed) or "none found"
        print(f"free inputs (istate 50/50): {shown}")

    result = sim.run()

    ## Final-state / attractor table. Still computed (--csv writes it), but not
    ## dumped to the terminal: it is one row with a column per attractor, so a
    ## network with any number of them prints as an unreadable wall. Only the
    ## count goes to the screen; use --csv to actually read them.
    probtraj = result.get_last_states_probtraj().rename(columns=lambda c: relabel_state(c, nmap))
    print(f"final states: {probtraj.shape[1]} (--csv writes the table)")

    ## Fixed points - the stable states and their basin sizes. Kept on screen:
    ## one row each, so this stays legible where the table above does not.
    fp = relabelled_fptable(result, nmap)
    if fp is not None:
        print(f"\nfixed points ({len(fp)}):")
        print(fp.to_string(index=False))

    out_stem = stem + mutation_tag(mutations)       # tag mutated runs so they don't clobber WT

    if args.csv:
        os.makedirs(args.csv, exist_ok=True)
        states_out = os.path.join(args.csv, f"{out_stem}_last_states.csv")
        probtraj.to_csv(states_out)
        print(f"-> {states_out}")
        if fp is not None:
            fp_out = os.path.join(args.csv, f"{out_stem}_fixpoints.csv")
            fp.to_csv(fp_out, index=False)
            print(f"-> {fp_out}")

    if args.plot:
        for p in save_plots(result, out_stem, args.plot, nmap,
                            (args.fig_width, args.fig_height)):
            print(f"-> {p}")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run MaBoSS on one or more .bnet models and report the "
                    "attractors, fixed points, and an optional node-trajectory "
                    "graph, under optional module knock-outs/gains.")
    parser.add_argument("bnet", nargs="+", help="MaBoSS .bnet model file(s)")
    parser.add_argument("--name-map", default=None,
                        help="6-cyto_to_neko.py <stem>_name_map.tsv, applied to every model; "
                             "omit to use each model's own <stem>_name_map.tsv if one sits "
                             "beside it")
    parser.add_argument("--mutate", action="append", default=[], metavar="NODE=ON|OFF",
                        help="lock a module ON (constitutive) or OFF (knock-out); repeatable")
    parser.add_argument("--seed", type=int, default=None, metavar="N",
                        help="MaBoSS seed_pseudorandom (default: MaBoSS's own)")
    parser.add_argument("--free-inputs", action="store_true",
                        help="start each self-input module at 50/50 to sample the whole landscape")
    parser.add_argument("--plot", metavar="OUTDIR", nargs="?", const=DEFAULT_OUTDIR, default=None,
                        help=f"save the node-trajectory graph, <stem>_node_traj.png "
                             f"(default: {DEFAULT_OUTDIR})")
    parser.add_argument("--fig-width", type=float, default=9.0, metavar="IN",
                        help="width of the --plot figure in inches (default 9)")
    parser.add_argument("--fig-height", type=float, default=5.0, metavar="IN",
                        help="height of the --plot figure in inches (default 5)")
    parser.add_argument("--csv", metavar="OUTDIR", nargs="?", const=DEFAULT_OUTDIR, default=None,
                        help=f"save <stem>_last_states.csv and <stem>_fixpoints.csv "
                             f"(default: {DEFAULT_OUTDIR})")
    parser.add_argument("--max-time", type=float, default=None,
                        help="MaBoSS max_time (default: leave the model's own value)")
    parser.add_argument("--sample-count", type=int, default=None,
                        help="MaBoSS sample_count / number of trajectories")
    parser.add_argument("--threads", type=int, default=None, metavar="N",
                        help="MaBoSS thread_count - trajectories are independent, so this "
                             "scales near-linearly (default: MaBoSS's own, which is 1)")
    args = parser.parse_args(argv)

    ran, failed = 0, []
    for bnet in args.bnet:
        if not os.path.isfile(bnet):
            print(f"{bnet}: FAILED - not a file")
            failed.append(bnet)
            continue
        try:
            run_one(bnet, args)                     # finds this model's own name map
            ran += 1
        except Exception as e:                      # keep going if one model fails
            print(f"{bnet}: FAILED - {e}")
            failed.append(bnet)

    print(f"\n{ran} run, {len(failed)} failed")
    if failed:
        print("Failed models:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
