#!/usr/bin/env python3
## Purpose: Simulate a treatment on a MaBoSS `.bnet` model - run a wild-type arm
## and a treated arm of the same model and report the per-module change in
## P(ON), which is the actual question when a drug targets a gene, and is far
## more readable than two sets of 180 trajectories. Shares its model-building
## code with 8-run_maboss.py (both import maboss_common.py), so the wild-type arm
## here is built exactly as a plain 8-run_maboss.py run.
##
## Each MODEL is a `.bnet` from 7-sif_to_bnet.py, found and labelled exactly as in
## 8-run_maboss.py - see that script's header for the name-map lookup. --mutate
## and --free-inputs work the same way here too, applied as the BACKGROUND (both
## arms) rather than the treatment.
##
## Treatment (--compare 18.37=OFF):
##   Both arms are built by the same code path so they can only differ by the
##   treatment. Note the treated arm is locked OFF for the whole run, i.e. a drug
##   held at dose; that is NOT the same experiment as starting the module at 0 and
##   letting its regulators switch it back on, which models a single transient dose.
##
##   The delta table carries each module's curated labels - Cluster Function and
##   Module_HGNCs, joined from data/cyto_exp/<CL>_display.csv by the cell-line code
##   the model filename starts with - so a module that moved arrives with its
##   function and gene list attached. --display-csv overrides the lookup.
##
##   Both delta_final (change at the end, the attractor readout) and delta_max
##   (largest change at any time, with t_delta_max) are reported: a treatment can
##   shift a module hard mid-run and still wash out by the end, which the final
##   value alone would report as "nothing happened".
##
##   --noise-floor re-runs wild-type at a second seed and treats the largest
##   wild-type-vs-wild-type difference as the floor below which a delta is
##   sampling noise. Needed because MaBoSS is deterministic per seed - repeating
##   the same command gives a bit-identical answer, so the floor has to come from
##   changing the seed, not from running twice.
##
## Files (--csv [OUTDIR]): <stem>_<tag>_vs_wt_delta.csv (the ranked table). <tag>
##   carries the --mutate background as well as the --compare treatment, so a
##   treatment run and the same treatment on a mutant background land in different
##   files rather than the second overwriting the first.
## Graphs (--plot [OUTDIR]): <stem>_<tag>_vs_wt_compare_traj.png (the --top
##   most-changed modules, wild-type solid vs treated dashed) and
##   <stem>_<tag>_vs_wt_compare_heatmap.png (the most-changed modules' delta over
##   time - at least 40 rows, more with a larger --top; the full set is in the
##   delta CSV, and on a 160-module network the rest is a block of zeroes).
##   --fig-width/--fig-height size both figures; leaving the height off keeps the
##   heatmap's own height, which scales with how many rows it is drawing.
## Both write to gabi-maboss/output/maboss/ when given bare; name a folder to
## override. Nothing is written unless the flag is present - a plain run just prints.
##
## Usage:
##   python compare_maboss.py MODEL.bnet [MODEL.bnet ...] --compare NODE=ON|OFF ...
##          [--name-map FILE] [--mutate NODE=ON|OFF ...] [--top N] [--noise-floor]
##          [--display-csv FILE] [--seed N] [--free-inputs]
##          [--plot [OUTDIR]] [--fig-width IN] [--fig-height IN] [--csv [OUTDIR]]
##          [--max-time T] [--sample-count N] [--threads N]
## Run in the activated `maboss` conda env so the `maboss` package AND the MaBoSS
## engine binary (`MaBoSS`, on the env's PATH) are both found.

import argparse
import os
import sys

import pandas as pd

from maboss_common import (DEFAULT_OUTDIR, build_sim, find_name_map, load_name_map,
                            mutation_tag, node_traj)


## Curated module labels, for answering "what does the module that changed
## actually do". build_display_csv.R writes one <CL>_display.csv per cell line;
## a model's cell line is the code its filename starts with (UH_reactome_... -> UH).
CELL_LINES = ("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")
DISPLAY_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "cyto_exp"))
ANNOT_COLS = ("Display", "Cluster Function", "Module_HGNCs")


def find_display_csv(stem, explicit):
    """Locate the display CSV for one model, or None.

    The cell line is taken from the leading code of the model filename, which is
    how every network in this pipeline is named. Returns None (rather than
    raising) when the code isn't one of the eight or the file isn't there - the
    comparison still works, it just goes out unannotated.
    """
    if explicit:
        return explicit
    code = stem.split("_")[0].upper()
    if code not in CELL_LINES:
        return None
    path = os.path.join(DISPLAY_DIR, f"{code}_display.csv")
    return path if os.path.isfile(path) else None


def load_annotation(path):
    """Module_ID -> curated labels, one row per module.

    The display CSV has a row per *gene*, so Module_ID repeats; the label columns
    are constant within a module, so the first row of each is the module's.
    """
    d = pd.read_csv(path).drop_duplicates(subset="Module_ID").set_index("Module_ID")
    return d[[c for c in ANNOT_COLS if c in d.columns]]


def compute_delta(wt, tx):
    """Per-module change in P(ON), treated minus wild-type.

    Returns (delta_over_time, summary). The summary carries both the change at
    the end of the run (the attractor readout) and the largest change at any
    time point - a module can be pushed hard mid-run and still recover, which the
    final value alone would hide, and which matters when the question is what a
    treatment does *over time*.
    """
    mods = wt.columns.intersection(tx.columns)
    times = wt.index.intersection(tx.index)
    wt, tx = wt.loc[times, mods], tx.loc[times, mods]

    d = tx - wt
    t_at_max = d.abs().idxmax()                       # time of the biggest swing
    summary = pd.DataFrame({
        "P_ON_wt": wt.iloc[-1],
        "P_ON_treated": tx.iloc[-1],
        "delta_final": d.iloc[-1],
        "delta_max": pd.Series([d.at[t_at_max[m], m] for m in mods], index=mods),
        "t_delta_max": t_at_max,
    })
    return d, summary.reindex(summary["delta_max"].abs().sort_values(ascending=False).index)


def annotate_delta(summary, annot, targets, floor):
    """Attach curated labels, the treatment target flag, and the noise verdict."""
    out = summary.copy()
    out.insert(0, "is_target", [m in targets for m in out.index])
    if annot is not None:
        joined = out.join(annot, how="left")
        for c in reversed([c for c in ANNOT_COLS if c in joined.columns]):
            out.insert(1, c, joined[c])
    if floor is not None:
        out["above_noise"] = out["delta_max"].abs() > floor
    out.index.name = "Module_ID"
    return out


def _label(mod, tbl):
    """Plot label for one module: its curated Display name if the join found one."""
    if "Display" in tbl.columns:
        disp = tbl.at[mod, "Display"]
        if isinstance(disp, str) and disp and disp != mod:
            return f"{mod} ({disp})"
    return str(mod)


def save_compare_plots(d, wt, tx, tbl, stem, outdir, top, width, height):
    """Write the two comparison graphs: paired trajectories, and a delta heatmap.

    Both are limited to the most-changed modules - the paired plot to `top`, the
    heatmap to at least 40 since it carries more rows legibly. Drawing all of them
    is what makes a 180-module figure unreadable, and the modules that didn't move
    are exactly the ones with nothing to show.

    `width`/`height` are --fig-width/--fig-height, in inches, applied to both
    figures. `height` is None unless asked for, because the two have different
    sensible defaults: the paired plot is a fixed 6, while the heatmap's grows
    with its row count so the labels stay apart. An explicit --fig-height wins
    over both.
    """
    import matplotlib
    matplotlib.use("Agg")                          # headless: write files, no display
    import matplotlib.pyplot as plt
    import numpy as np

    os.makedirs(outdir, exist_ok=True)
    written = []
    sel = list(tbl.head(top).index)

    ## 1. Wild-type vs treated for the most-changed modules. One colour per
    ## module, solid = wild-type and dashed = treated, so each pair reads as one
    ## module rather than as two unrelated lines.
    fig, ax = plt.subplots(figsize=(width, height or 6.0))
    colours = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(sel), 1)))
    for colour, mod in zip(colours, sel):
        ax.plot(wt.index, wt[mod], color=colour, lw=1.6, label=_label(mod, tbl))
        ax.plot(tx.index, tx[mod], color=colour, lw=1.6, ls="--")
    ax.set_xlabel("time"); ax.set_ylabel("P(module = ON)"); ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{stem}: {len(sel)} most-changed modules, wild-type vs treated")
    ## Two legends stacked down the right margin. Both anchor to the TOP: the
    ## module legend used to be centred and the arms legend bottom-aligned, which
    ## collided as soon as --top made the module list tall enough to reach the
    ## bottom of the axes. Anchored this way the module list grows downwards into
    ## empty margin instead, whatever --top is.
    arms = [plt.Line2D([], [], color="0.3", lw=1.6, label="wild-type"),
            plt.Line2D([], [], color="0.3", lw=1.6, ls="--", label="treated")]
    arms_legend = ax.legend(handles=arms, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    ax.add_artist(arms_legend)
    mods_legend = ax.legend(loc="upper left", bbox_to_anchor=(1.0, 0.88),
                            fontsize=8, title="module")
    p = os.path.join(outdir, f"{stem}_compare_traj.png")
    ## Both legends sit outside the axes, and a legend re-added with add_artist is
    ## not picked up by the tight bbox on its own - name them explicitly or long
    ## module labels get sliced off at the figure edge.
    fig.savefig(p, bbox_inches="tight", dpi=150,
                bbox_extra_artists=(mods_legend, arms_legend))
    plt.close(fig); written.append(p)

    ## 2. The most-changed modules' delta over time as a heatmap, ordered by size
    ## of change. This is the view that scales - many more rows are legible here
    ## than as lines - so it shows more than the paired plot above, but not every
    ## module: on a 160-module network only a handful move, and drawing the rest
    ## is a block of zeroes that also pushes the row count past the point where
    ## labels fit, leaving unreadable bands. Floor it at 40 so a small --top still
    ## gets the wider view the heatmap is for.
    order = list(tbl.head(max(top, 40)).index)
    M = d[order].T
    labelled = len(order) <= 40                    # labels stop being legible past this
    ## Only reserve per-row height when the rows are actually labelled; otherwise
    ## a 160-module figure is mostly empty space.
    per_row = 0.16 if labelled else 0.06
    auto_height = max(3.0, min(per_row * len(order) + 1.5, 18.0))
    fig, ax = plt.subplots(figsize=(width, height or auto_height))

    ## Scale the colour to the data, symmetric about zero. A fixed +-1 scale is
    ## the honest one but renders a uniformly white plot whenever the treatment
    ## shifts things by a tenth or so - which is the normal case - so the range
    ## is stated on the colourbar instead of being assumed.
    lim = max(float(M.abs().max().max()), 0.02)
    ## shading="nearest" takes cell CENTRES, so the coordinate arrays match the
    ## data shape - with "flat" they would each have to be one longer.
    im = ax.pcolormesh(M.columns.values, np.arange(len(order)), M.values,
                       cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    ax.set_xlabel("time"); ax.set_ylabel("module (most changed at top)")
    ax.invert_yaxis()
    if labelled:
        ax.set_yticks(np.arange(len(order)))
        ax.set_yticklabels([_label(m, tbl) for m in order], fontsize=7)
    else:
        ax.set_yticks([])
    ax.set_title(f"{stem}: change in P(module = ON), treated - wild-type")
    fig.colorbar(im, ax=ax, label=f"Δ P(ON), scale ±{lim:.3f}", fraction=0.03, pad=0.02)
    p = os.path.join(outdir, f"{stem}_compare_heatmap.png")
    fig.savefig(p, bbox_inches="tight", dpi=150); plt.close(fig); written.append(p)

    return written


def run_one(bnet, args):
    """Load, run, and report one model's wild-type-vs-treated comparison.

    --mutate (if given) is background applied to BOTH arms; --compare is the
    treatment, applied only to the treated arm. So the delta isolates the
    treatment even when the model is being run on a mutant background.
    """
    map_path = find_name_map(bnet, args.name_map)
    nmap = load_name_map(map_path) if map_path else {}
    rev  = {orig: safe for safe, orig in nmap.items()}   # module ID -> safe name, for --mutate

    stem = os.path.splitext(os.path.basename(bnet))[0]
    print(f"\n=== {bnet} ===")

    if map_path:
        how = "given" if args.name_map else "found beside the model"
        print(f"name map ({how}): {map_path}")
    else:
        print("name map: none found - labels stay as the MaBoSS-safe node names")

    wt_sim, bg, freed = build_sim(bnet, args, args.mutate, rev, nmap, args.seed)
    tx_sim, tx_muts, _ = build_sim(bnet, args, args.mutate + args.compare, rev, nmap, args.seed)
    treatment = tx_muts[len(bg):]

    for disp, state in bg:
        print(f"background mutation (both arms): {disp} {state}")
    print("treatment (treated arm only): "
          + ", ".join(f"{d} {s}" for d, s in treatment))
    if args.free_inputs:
        shown = ", ".join(nmap.get(n, n) for n in freed) or "none found"
        print(f"free inputs (istate 50/50): {shown}")

    wt = node_traj(wt_sim.run(), nmap)
    tx = node_traj(tx_sim.run(), nmap)
    d, summary = compute_delta(wt, tx)

    ## Noise floor: the same wild-type arm at a different seed. Two runs at the
    ## SAME seed are bit-identical (MaBoSS is deterministic per seed), so a plain
    ## repeat would report a floor of zero and flatter every small delta.
    floor = None
    if args.noise_floor:
        alt_sim, _, _ = build_sim(bnet, args, args.mutate, rev, nmap, (args.seed or 0) + 1)
        alt = node_traj(alt_sim.run(), nmap)
        floor = float((alt - wt).abs().max().max())
        print(f"noise floor (wild-type at two seeds): max |Δ| = {floor:.4f}")

    ## Which module names the treatment actually hit, for the is_target flag.
    targets = {disp for disp, _ in treatment}
    annot_path = find_display_csv(stem, args.display_csv)
    annot = load_annotation(annot_path) if annot_path else None
    print(f"module labels: {annot_path or 'none found - delta table goes out unannotated'}")

    tbl = annotate_delta(summary, annot, targets, floor)

    moved = tbl[tbl["delta_max"].abs() > (floor or 0)]
    print(f"\nmodules changed ({len(moved)} of {len(tbl)} above "
          f"{'the noise floor' if floor else 'zero'}), top {min(args.top, len(tbl))}:")
    cols = [c for c in ("Display", "Cluster Function", "P_ON_wt", "P_ON_treated",
                        "delta_final", "delta_max") if c in tbl.columns]
    shown = tbl.head(args.top)[cols].copy()
    if "Cluster Function" in shown:                # keep the terminal narrow
        shown["Cluster Function"] = shown["Cluster Function"].fillna("").str.slice(0, 32)
    print(shown.to_string(float_format=lambda v: f"{v:6.3f}"))

    ## The background goes in the filename as well as the treatment. Tagging by
    ## treatment alone means `--compare X=ON` and `--compare X=ON --mutate Y=OFF`
    ## write to the same stem, so a batch running both silently overwrites the
    ## first with the second - and the survivor carries no sign of which it is.
    ## mutation_tag([]) is "", so runs without a background keep their old names.
    out_stem = stem + mutation_tag(bg) + mutation_tag(treatment) + "_vs_wt"

    if args.csv:
        os.makedirs(args.csv, exist_ok=True)
        p = os.path.join(args.csv, f"{out_stem}_delta.csv")
        tbl.to_csv(p)
        print(f"-> {p}")
    if args.plot:
        for p in save_compare_plots(d, wt, tx, tbl, out_stem, args.plot, args.top,
                                    args.fig_width, args.fig_height):
            print(f"-> {p}")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Simulate a treatment on one or more MaBoSS .bnet models: run "
                    "wild-type and treated arms and report the per-module change "
                    "in P(ON).")
    parser.add_argument("bnet", nargs="+", help="MaBoSS .bnet model file(s)")
    parser.add_argument("--name-map", default=None,
                        help="6-cyto_to_neko.py <stem>_name_map.tsv, applied to every model; "
                             "omit to use each model's own <stem>_name_map.tsv if one sits "
                             "beside it")
    parser.add_argument("--mutate", action="append", default=[], metavar="NODE=ON|OFF",
                        help="background applied to BOTH arms (not the treatment); repeatable")
    parser.add_argument("--compare", action="append", default=[], required=True,
                        metavar="NODE=ON|OFF",
                        help="the treatment: locks a module ON or OFF in the treated arm "
                             "only; repeatable for a combination")
    parser.add_argument("--top", type=int, default=12, metavar="N",
                        help="how many of the most-changed modules to print and plot (default 12)")
    parser.add_argument("--noise-floor", action="store_true",
                        help="also run wild-type at a second seed and treat the largest "
                             "wild-type-vs-wild-type change as the floor below which a delta "
                             "is sampling noise")
    parser.add_argument("--display-csv", default=None, metavar="FILE",
                        help="curated module labels; omit to use data/cyto_exp/<CL>_display.csv "
                             "for the cell line the model filename starts with")
    parser.add_argument("--seed", type=int, default=None, metavar="N",
                        help="MaBoSS seed_pseudorandom, applied to every arm (default: MaBoSS's)")
    parser.add_argument("--free-inputs", action="store_true",
                        help="start each self-input module at 50/50 to sample the whole landscape")
    parser.add_argument("--plot", metavar="OUTDIR", nargs="?", const=DEFAULT_OUTDIR, default=None,
                        help=f"save the comparison graphs (default: {DEFAULT_OUTDIR})")
    parser.add_argument("--fig-width", type=float, default=10.0, metavar="IN",
                        help="width of both --plot figures in inches (default 10)")
    parser.add_argument("--fig-height", type=float, default=None, metavar="IN",
                        help="height of both --plot figures in inches; omit to keep the "
                             "defaults (6 for the paired plot, row-count-scaled for the heatmap)")
    parser.add_argument("--csv", metavar="OUTDIR", nargs="?", const=DEFAULT_OUTDIR, default=None,
                        help=f"save the ranked delta table (default: {DEFAULT_OUTDIR})")
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
