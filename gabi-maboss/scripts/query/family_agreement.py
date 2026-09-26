#!/usr/bin/env python3
## Purpose: measure how far the backbone and bidirectional model families agree,
## and write the answer out as tables the dissertation can cite.
##
## THE QUESTION. Every network here exists as two models. The BACKBONE drops the
## edges Gabi could not orient; the BIDIRECTIONAL keeps them as reciprocal pairs.
## The bi models were carried forward, and that choice needs evidence: where both
## families contain a module, do they say the same thing about it?
##
## WHAT IT COMPARES. For one perturbation run there are two delta tables, one per
## family, each with a row per module and a `delta_final` column - the change in
## P(ON) from wild type at the end of the run. Modules are matched by ID, so the
## comparison is over the INTERSECTION of the two models' nodes, which is the
## backbone's node set (it is the smaller model; the modules only the bi model
## has are exactly what the backbone cannot report on, and are counted here but
## compared nowhere, since there is nothing to compare them with).
##
## THE THREE POPULATIONS, and why the distinction decides the answer:
##
##   shared     every module both families contain
##   unmoved    |delta| < THRESHOLD - counted per family as well as for both at
##              once, since a module the backbone leaves alone while the bi model
##              moves it is a disagreement of a kind the sign columns cannot show
##   moved      |delta| >= THRESHOLD in BOTH - the only modules on which the two
##              families can be said to agree or disagree about direction
##
## Correlation is reported over both `shared` and `moved`, because the first is
## misleading on its own and the difference is the reason:
##
##   In a typical run most modules do not move at all. For the A498 18.48 ON run,
##   157 of 180 shared modules sit still in both families, so a correlation over
##   all shared modules is largely zero agreeing with zero. The clearest case is
##   the CAKI-1 21.36 knockout, which moves NOTHING in either family and still
##   scores r = 0.911 over all shared modules - a correlation between one family's
##   sampling noise and the other's. Over the moved set it is undefined, which is
##   the honest answer. Quote `r_moved`; `r_shared` is kept only so the gap can
##   be seen.
##
## Sign agreement is over `moved` alone, for the same reason: a module that did
## not move has no direction to agree about.
##
## WHAT THIS DOES NOT ESTABLISH, and it matters. The comparison covers only the
## eleven TARGETED runs, because those are the only perturbations simulated on
## both families. They were chosen from topology - mostly high-betweenness hubs -
## so they are neither a random nor a complete sample of the ~200 single-module
## perturbations available per model, and they are biased towards perturbations
## with large effects. A complete answer needs `sweep_maboss.py` run over the
## backbone models as well as the bi models, and this script pointed at the two
## sweep tables instead; that has not been done, and the dissertation reports the
## gap as a limitation rather than papering over it. The `--sweep` path below is
## written for that comparison so it is a re-run and not a rewrite when the
## backbone sweep exists.
##
## Pearson, not Spearman: the quantity of interest is whether the two families
## report the same SIZE of change, not merely the same ranking, and both vectors
## are on the same bounded scale already.
##
## Outputs (to gabi-maboss/output/maboss/family_agreement/ unless --outdir):
##   family_agreement_per_run.csv    one row per run - the numbers behind the
##                                   dissertation's model-family table
##   family_agreement_disagreements.csv  every module whose direction differs,
##                                   with both deltas, so each can be looked at
##
## Runs outside the `maboss` conda env: it reads finished delta tables and
## simulates nothing.
##
## Usage:
##   python family_agreement.py                      # the eleven targeted runs
##   python family_agreement.py --threshold 0.1
##   python family_agreement.py --sweep              # when both families are swept

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GM = os.path.abspath(os.path.join(HERE, "..", ".."))
TARGETED = os.path.join(GM, "output", "maboss", "top_6-time_20")
SWEEP = os.path.join(GM, "output", "maboss", "sweep")
OUTDIR = os.path.join(GM, "output", "maboss", "family_agreement")

## The paired model stems, per cell line: (backbone, bidirectional). A498's
## carry the `drop_unsigned_edge` tag because that network had eight unsignable
## edges to drop; the other two did not.
FAMILIES = {
    "A498": ("AH_reactome_backbone_subnet_drop_unsigned_edge",
             "AH_reactome_bi_subnet_drop_unsigned_edge"),
    "CAKI-1": ("CH_reactome_backbone_subnet_edge",
               "CH_reactome_bi_subnet_edge"),
    "UMRC2": ("UH_reactome_backbone_subnet_edge",
              "UH_reactome_bi_subnet_edge"),
}

## Module IDs are <Cluster_Number>.<Multi_Profiles>; an all-numeric one becomes a
## float unless forced to str, turning 18.30 into 18.3 and matching nothing.
IDCOL = {"Module_ID": str, "perturbed": str, "readout": str}


## One run applies two perturbations, and its filename cannot say which is the
## treatment and which the background - both are just <module><STATE> segments.
## It is named explicitly rather than parsed, so the label matches the
## dissertation's (treatment first, background second) instead of file order.
COMBINED = {"18_36OFF_26_8ON": "26.8 ON / 18.36 OFF"}


def pretty_run(tag):
    """'18_3_12OFF' -> '18.3_12 OFF'.

    Only the FIRST underscore is the decimal point: a sanitised module is
    <cluster>_<profiles>, so 18_3_12 is 18.3_12 and not 18.3.12. A tag with two
    perturbations in it is looked up rather than parsed - see COMBINED.
    """
    if tag in COMBINED:
        return COMBINED[tag]
    m = re.match(r"^(.*?)(ON|OFF)$", tag)
    if not m:
        return tag
    if re.search(r"(ON|OFF)_", tag):
        raise SystemExit(f"{tag!r} holds more than one perturbation and is not "
                         f"in COMBINED - add it rather than letting it mangle")
    return m.group(1).replace("_", ".", 1).rstrip(".") + " " + m.group(2)


def corr(x, y):
    """Pearson r, or NaN where it is not defined rather than a spurious value."""
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def compare(x, y, threshold):
    """The three populations and the numbers taken over them."""
    moved = (np.abs(x) >= threshold) & (np.abs(y) >= threshold)
    quiet_bb = np.abs(x) < threshold
    quiet_bi = np.abs(y) < threshold
    n = int(moved.sum())
    same = int((np.sign(x[moved]) == np.sign(y[moved])).sum()) if n else 0
    return {
        "n_shared": len(x),
        ## Per family as well as jointly: the gap between the two is the count of
        ## modules one family moves and the other does not, which is a real
        ## disagreement that never reaches the sign comparison.
        "n_unmoved_backbone": int(quiet_bb.sum()),
        "n_unmoved_bi": int(quiet_bi.sum()),
        "n_unmoved_both": int((quiet_bb & quiet_bi).sum()),
        "n_moved_both": n,
        "r_shared": corr(x, y),
        "r_moved": corr(x[moved], y[moved]),
        "n_same_direction": same,
        "pct_same_direction": (100.0 * same / n) if n else np.nan,
        "n_disagree": n - same,
    }, moved


def targeted_pairs():
    """Every run simulated on both families, as (line, tag, backbone, bi) paths."""
    out = []
    for line, (bb, bi) in FAMILIES.items():
        for f in sorted(glob.glob(os.path.join(TARGETED,
                                               f"{bb}_*_vs_wt_delta.csv"))):
            tag = os.path.basename(f)[len(bb) + 1:-len("_vs_wt_delta.csv")]
            mate = os.path.join(TARGETED, f"{bi}_{tag}_vs_wt_delta.csv")
            if os.path.isfile(mate):
                out.append((line, tag, f, mate))
            else:
                print(f"  note: {line} {tag} has no bidirectional counterpart",
                      file=sys.stderr)
    return out


def run_targeted(threshold):
    rows, bad = [], []
    for line, tag, fbb, fbi in targeted_pairs():
        a = pd.read_csv(fbb, dtype=IDCOL).set_index("Module_ID")
        b = pd.read_csv(fbi, dtype=IDCOL).set_index("Module_ID")
        shared = a.index.intersection(b.index)
        x = a.loc[shared, "delta_final"].to_numpy()
        y = b.loc[shared, "delta_final"].to_numpy()
        stats, moved = compare(x, y, threshold)
        stats.update(cell_line=line, run=pretty_run(tag), run_tag=tag,
                     n_bi_only=len(b.index.difference(a.index)))
        rows.append(stats)
        for mod, xx, yy in zip(shared[moved], x[moved], y[moved]):
            if np.sign(xx) != np.sign(yy):
                bad.append({"cell_line": line, "run": pretty_run(tag),
                            "module": mod,
                            "genes": b.loc[mod, "Module_HGNCs"]
                            if "Module_HGNCs" in b.columns else "",
                            "delta_backbone": round(float(xx), 4),
                            "delta_bidirectional": round(float(yy), 4)})
    cols = ["cell_line", "run", "run_tag", "n_shared", "n_bi_only",
            "n_unmoved_backbone", "n_unmoved_bi", "n_unmoved_both",
            "n_moved_both", "r_shared", "r_moved",
            "n_same_direction", "pct_same_direction", "n_disagree"]
    return pd.DataFrame(rows)[cols], pd.DataFrame(bad)


def run_sweep(threshold):
    """The complete comparison, once both families have been swept.

    Not reachable today: `sweep_maboss.py` has only been run on the bidirectional
    models, so there is no backbone sweep to read. Written now so that producing
    one is a re-run of this script rather than a new piece of analysis.
    """
    rows = []
    for line, (bb, bi) in FAMILIES.items():
        fbb = os.path.join(SWEEP, f"{bb}_sweep_OFF.csv")
        fbi = os.path.join(SWEEP, f"{bi}_sweep_OFF.csv")
        if not (os.path.isfile(fbb) and os.path.isfile(fbi)):
            print(f"  {line}: no backbone sweep at {os.path.basename(fbb)} "
                  f"- run sweep_maboss.py on the backbone model first",
                  file=sys.stderr)
            continue
        a = pd.read_csv(fbb, dtype=IDCOL)
        b = pd.read_csv(fbi, dtype=IDCOL)
        key = ["perturbed", "readout"]
        m = a[key + ["delta_final"]].merge(b[key + ["delta_final"]], on=key,
                                           suffixes=("_bb", "_bi"))
        stats, _ = compare(m.delta_final_bb.to_numpy(),
                           m.delta_final_bi.to_numpy(), threshold)
        stats.update(cell_line=line, run="whole sweep", run_tag="sweep_OFF",
                     n_bi_only=np.nan)
        rows.append(stats)
    return pd.DataFrame(rows), pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=0.05,
                    help="a module counts as moved at |delta_final| >= this "
                         "(default 0.05, the value used throughout)")
    ap.add_argument("--sweep", action="store_true",
                    help="compare whole sweeps instead of the targeted runs "
                         "(needs sweep_maboss.py to have been run on both)")
    ap.add_argument("--outdir", default=OUTDIR)
    a = ap.parse_args()

    per_run, disagree = (run_sweep(a.threshold) if a.sweep
                         else run_targeted(a.threshold))
    if per_run.empty:
        raise SystemExit("nothing to compare")

    os.makedirs(a.outdir, exist_ok=True)
    p1 = os.path.join(a.outdir, "family_agreement_per_run.csv")
    per_run.to_csv(p1, index=False)
    print(f"wrote {p1}  ({len(per_run)} runs, threshold {a.threshold})")
    if not disagree.empty:
        p2 = os.path.join(a.outdir, "family_agreement_disagreements.csv")
        disagree.to_csv(p2, index=False)
        print(f"wrote {p2}  ({len(disagree)} module responses)")

    def fmt(v):
        return "  -  " if pd.isna(v) else f"{v:5.3f}"

    print()
    print(f"{'Line':8s} {'Run':22s} {'shared':>6s} {'un.bb':>6s} {'un.bi':>6s} "
          f"{'moved':>6s} {'r shared':>8s} {'r moved':>8s} {'same dir':>9s}")
    for _, r in per_run.iterrows():
        pct = ("   -  " if pd.isna(r.pct_same_direction)
               else f"{r.pct_same_direction:5.1f}%")
        print(f"{r.cell_line:8s} {r.run:22s} {r.n_shared:6d} "
              f"{r.n_unmoved_backbone:6d} {r.n_unmoved_bi:6d} "
              f"{r.n_moved_both:6d} "
              f"{fmt(r.r_shared):>8s} {fmt(r.r_moved):>8s} {pct:>9s}")

    tot = int(per_run.n_moved_both.sum())
    dis = int(per_run.n_disagree.sum())
    print()
    print(f"across {len(per_run)} runs: {tot} module responses move in both "
          f"families, {dis} differ in direction "
          f"({100 * (tot - dis) / tot:.1f}% agree)")
    if not disagree.empty:
        rep = disagree.module.value_counts()
        multi = rep[rep > 1]
        if len(multi):
            print("  modules disagreeing in more than one run: "
                  + ", ".join(f"{m} ({n})" for m, n in multi.items()))
    print()
    print("NB the targeted runs are a topology-chosen sample of the possible "
          "perturbations,\n    not a random or complete one - see the header, "
          "and --sweep for the full test.")


if __name__ == "__main__":
    main()
