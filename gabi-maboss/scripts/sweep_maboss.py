#!/usr/bin/env python3
## Purpose: Perturb EVERY module in a model, one at a time, and write the whole
## result out as a long table - the systematic version of compare_maboss.py, which
## answers one perturbation at a time. Where that script asks "what does knocking
## out CREB1 do", this builds the table that answers "what happens to the apoptosis
## modules when any module is knocked out, in each network", which is a query rather
## than a run.
##
## The numbers are compare_maboss.py's numbers: the wild-type arm is built by the
## same maboss_common.build_sim(), and each perturbed arm's delta comes from the
## same compare_maboss.compute_delta(). So a row of this table and a single
## compare_maboss.py run of the same perturbation agree by construction, and the
## targeted runs already written up stay comparable with the sweep.
##
## What makes the sweep affordable is that the wild-type arm does not depend on the
## perturbation: compare_maboss.py runs it again for every treatment, this runs it
## once per model and reuses it for all ~200 perturbed arms. The perturbed arms are
## independent, so they are spread over --workers processes (MaBoSS's own threading
## scales poorly on these models - one thread per worker and many workers in
## parallel is several times faster than the reverse).
##
## Output (three CSVs per model, into output/maboss/sweep/ unless --outdir):
##   <stem>_sweep_<STATE>.csv          the database: one row per
##       (perturbed module, readout module) pair - P(ON) wild-type and perturbed,
##       delta_final and delta_max, and the topology columns below. ~200 x 200 rows
##       per model, and the three models' files concatenate (the `model` and
##       `cell_line` columns keep them apart).
##   <stem>_sweep_<STATE>_summary.csv  one row per perturbed module: how many
##       modules it moved, how far, and how many it can reach at all - the
##       "which knockouts matter" ranking.
##   <stem>_sweep_modules.csv          one row per module: curated labels, wild-type
##       P(ON), in/out degree. The dimension table the other two join back to, and
##       the only place the full Module_HGNCs gene list is kept.
##   <stem>_sweep_<STATE>_traj.npz     every arm's P(ON) curves - `wt` (time x
##       module), `traj` (perturbation x time x module), and the `perturbed`,
##       `modules` and `times` labels. The three CSVs summarise each pair to a few
##       numbers, which cannot be drawn as a curve; this is what lets any
##       perturbation be plotted over time afterwards, or a wild-type-vs-perturbed
##       pair redrawn, without paying for the simulation again. A few MB per model
##       (float32, compressed). --no-traj skips it.
##
## Topology columns come from the `.bnet` itself, via trace_cascade.py's parser, so
## every pair arrives with `hop` (shortest directed distance from the perturbed
## module, 1 = direct target, blank = unreachable) and `path_sign_agrees` (`ok` /
## `FLIP` - whether the sign product along that path predicts the direction the
## module actually moved). That is what separates a direct effect from the far end
## of a cascade when querying, and it costs nothing to compute.
##
## The noise floor is always measured (wild-type at two seeds, one extra run out of
## ~200) and applied as the `above_noise` flag, since a sweep is exactly where an
## unfiltered table of near-zero deltas would be read as signal.
##
## Usage:
##   python sweep_maboss.py MODEL.bnet [MODEL.bnet ...] [--state OFF|ON]
##          [--workers N] [--mutate NODE=ON|OFF ...] [--free-inputs] [--seed N]
##          [--only MODULE ...] [--resume] [--no-traj] [--checkpoint N]
##          [--outdir DIR] [--display-csv FILE] [--name-map FILE]
##          [--max-time T] [--sample-count N] [--threads N]
## Run in the activated `maboss` conda env so the `maboss` package AND the MaBoSS
## engine binary (`MaBoSS`, on the env's PATH) are both found.

import argparse
import multiprocessing as mp
import os
import sys
import time

import pandas as pd

from maboss_common import build_sim, find_name_map, load_name_map, node_traj
from compare_maboss import (CELL_LINES, compute_delta, find_display_csv,
                            load_annotation)
from trace_cascade import breadth_first, forward_edges, load_rules


## A sweep writes three files per model and is a different kind of artefact from a
## targeted compare_maboss.py run, so it gets its own subfolder of what MaBoSS
## produces rather than being scattered among the per-treatment deltas.
DEFAULT_SWEEP_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "output", "maboss", "sweep"))

## Module IDs are `<Cluster_Number>.<Multi_Profiles>`, so plenty of them (21.14,
## 18.37) look like numbers. Every read-back of these columns forces text - left to
## infer, a table whose IDs happened to be all-numeric would come back as floats,
## which silently loses `18.30` to `18.3` and matches nothing on join. Anything
## else querying these files should do the same: `dtype=str` in pandas,
## `colClasses = "character"` in R.
ID_COLS = {"model": str, "cell_line": str, "state": str,
           "perturbed": str, "readout": str}

## Rough peak memory of one worker: MaBoSS's own run is small, but reading its
## probtraj back into a per-module table is not, and that peak is what decides how
## many arms can run at once. Cores are rarely the limit - a box with more cores
## than (RAM / this) will swap and lock up long before it runs out of them. It is a
## rule of thumb from the models here (~160-200 modules); check actual worker RSS
## on the first few runs and raise --workers if there is headroom.
WORKER_GB = 2.5

## Column order of the long table, fixed so appended chunks line up and so a
## resumed run keeps writing the same shape.
COLUMNS = ["model", "cell_line", "state",
           "perturbed", "perturbed_display", "perturbed_function",
           "readout", "readout_display", "readout_function",
           "P_ON_wt", "P_ON_perturbed", "delta_final", "delta_max", "t_delta_max",
           "hop", "path_sign_agrees", "above_noise", "is_target"]


def default_workers():
    """Cores, but capped by memory - the constraint that actually bites.

    Reading MemAvailable rather than total: a sweep is normally started on a box
    already doing something else, and it is that free headroom the workers have to
    fit inside. Falls back to the core count where /proc isn't there.
    """
    cores = os.cpu_count() or 1
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    available_gb = int(line.split()[1]) / 1024 / 1024
                    return max(1, min(cores, int(available_gb / WORKER_GB)))
    except OSError:
        pass
    return cores


## --------------------------------------------------------------------------
## Worker side. One process per core, each holding the model path and the shared
## wild-type trajectory; the perturbed arms are independent runs, so nothing
## crosses between them.
## --------------------------------------------------------------------------

_W = {}


def _init_worker(bnet, args, rev, nmap, wt):
    _W.update(bnet=bnet, args=args, rev=rev, nmap=nmap, wt=wt)


def _run_perturbation(safe_name):
    """Run one perturbed arm; return its delta table and its P(ON) trajectory.

    Returns (safe_name, (summary, trajectory)) or (safe_name, error string) - a
    model node whose perturbation MaBoSS refuses shouldn't cost the other 200.

    The trajectory is reindexed onto the wild-type arm's own times and modules
    before it leaves the worker, so every perturbation's array lands on the same
    grid and the stack can be a plain 3-D array rather than 200 ragged frames.
    """
    try:
        spec = [f"{safe_name}={_W['args'].state}"]
        sim, _, _ = build_sim(_W["bnet"], _W["args"], _W["args"].mutate + spec,
                              _W["rev"], _W["nmap"], _W["args"].seed)
        result = sim.run()
        tx = node_traj(result, _W["nmap"])
        del result                                    # drops its temp dir now
        _, summary = compute_delta(_W["wt"], tx)
        wt = _W["wt"]
        grid = (tx.reindex(index=wt.index, columns=wt.columns)
                  .to_numpy(dtype="float32") if _W["args"].trajectories else None)
        return safe_name, (summary, grid)
    except Exception as e:                            # pragma: no cover - engine-side
        return safe_name, f"FAILED - {e}"


## --------------------------------------------------------------------------
## Parent side.
## --------------------------------------------------------------------------

def save_trajectories(path, wt, traj):
    """Write the wild-type arm and every perturbed arm's P(ON) curves to one .npz.

    The delta table keeps only summary statistics per pair, which cannot be drawn
    as a curve; storing the trajectories is what lets any perturbation be plotted
    over time afterwards without paying for the simulation again. float32 is well
    inside MaBoSS's own sampling noise (~1e-4 here) and halves the file.

    Written as one array rather than a file per perturbation: ~200 x 40 x 200
    float32 is only a few MB compressed, and it loads in a single call.
    """
    import numpy as np

    perturbed = sorted(traj)
    np.savez_compressed(
        path,
        wt=wt.to_numpy(dtype="float32"),
        traj=np.stack([traj[m] for m in perturbed]) if perturbed
             else np.empty((0, len(wt.index), len(wt.columns)), dtype="float32"),
        perturbed=np.array(perturbed, dtype=str),
        modules=np.array(list(wt.columns), dtype=str),
        times=wt.index.to_numpy(dtype="float64"))


def load_trajectories(path):
    """Perturbed module -> (time x module) array, from a previous run's .npz.

    Only used to carry a --resume forward: the trajectories already computed are
    read back so the rewritten file keeps them.
    """
    import numpy as np

    with np.load(path, allow_pickle=False) as z:
        return {m: z["traj"][i] for i, m in enumerate(z["perturbed"])}


def degrees(bnet):
    """(in_degree, out_degree) per safe node name, from the model's own rules."""
    regulators = load_rules(bnet)
    forward = forward_edges(regulators)
    indeg = {t: len(r) for t, r in regulators.items()}
    outdeg = {n: len(t) for n, t in forward.items()}
    return indeg, outdeg


def module_table(wt, nmap, annot, indeg, outdeg, freed):
    """The dimension table: one row per module, with everything constant about it.

    Keyed on the module ID the other two tables use, but keeping the safe name too,
    since that is what the `.bnet` and any MaBoSS output outside this script speak.
    """
    rows = []
    for safe, module in sorted(nmap.items(), key=lambda kv: kv[1]):
        rows.append({
            "Module_ID": module,
            "safe_name": safe,
            "P_ON_wt": float(wt[module].iloc[-1]) if module in wt.columns else None,
            "in_degree": indeg.get(safe, 0),
            "out_degree": outdeg.get(safe, 0),
            "is_free_input": safe in freed,
        })
    tbl = pd.DataFrame(rows).set_index("Module_ID")
    if annot is not None:
        tbl = tbl.join(annot, how="left")
    return tbl.reset_index()


def annotate(summary, source_module, dist, sign, ctx):
    """Turn one perturbation's delta summary into rows of the long table."""
    out = summary.copy()
    out.index.name = "readout"
    out = out.reset_index().rename(columns={"P_ON_treated": "P_ON_perturbed"})

    out.insert(0, "model", ctx["stem"])
    out.insert(1, "cell_line", ctx["cell_line"])
    out.insert(2, "state", ctx["state"])
    out.insert(3, "perturbed", source_module)

    rev, direction = ctx["rev"], ctx["direction"]

    ## hop / path_sign_agrees use trace_cascade.py's convention: the perturbation's
    ## own direction times the sign product along the shortest path is the predicted
    ## direction, and a mixed-sign edge anywhere on it (product 0) predicts nothing.
    hops, verdicts = [], []
    for module, delta in zip(out["readout"], out["delta_final"]):
        safe = rev.get(module, module)
        hop = dist.get(safe)
        hops.append(hop)
        if hop is None:
            verdicts.append("")
            continue
        predicted = sign[safe] * direction
        verdicts.append("" if predicted == 0
                        else ("ok" if (predicted > 0) == (delta > 0) else "FLIP"))
    out["hop"] = hops
    out["path_sign_agrees"] = verdicts
    out["above_noise"] = out["delta_max"].abs() > ctx["floor"]
    out["is_target"] = out["readout"] == source_module

    annot = ctx["annot"]
    if annot is not None:
        labels = annot.reindex(out["readout"])
        for col, name in (("Display", "readout_display"),
                          ("Cluster Function", "readout_function")):
            out[name] = labels[col].values if col in labels.columns else None
        src = annot.reindex([source_module])
        for col, name in (("Display", "perturbed_display"),
                          ("Cluster Function", "perturbed_function")):
            out[name] = src[col].iloc[0] if col in src.columns else None
    return out.reindex(columns=COLUMNS)


def summarise(delta_csv, modules, ctx):
    """Per-perturbation summary, read back from the long table on disk.

    Reading the file rather than accumulating in memory is what makes --resume
    correct: the summary then covers everything in the table, not just the
    perturbations this invocation happened to run.
    """
    long = pd.read_csv(delta_csv, dtype=ID_COLS)
    off_target = long[~long["is_target"]]
    grouped = off_target.groupby("perturbed")

    summary = pd.DataFrame({
        "n_moved": off_target[off_target["above_noise"]].groupby("perturbed").size(),
        "n_readouts": grouped.size(),
        "n_reachable": grouped["hop"].apply(lambda s: s.notna().sum()),
        "max_abs_delta_final": grouped["delta_final"].apply(lambda s: s.abs().max()),
        "max_abs_delta_max": grouped["delta_max"].apply(lambda s: s.abs().max()),
        "sum_abs_delta_final": grouped["delta_final"].apply(lambda s: s.abs().sum()),
        "n_flip": off_target[off_target["path_sign_agrees"] == "FLIP"]
                  .groupby("perturbed").size(),
        "n_up": off_target[off_target["above_noise"] & (off_target["delta_final"] > 0)]
                .groupby("perturbed").size(),
        "n_down": off_target[off_target["above_noise"] & (off_target["delta_final"] < 0)]
                  .groupby("perturbed").size(),
    })
    ## Perturbations that moved nothing have no rows in the filtered groupings, so
    ## they come back as NaN rather than 0 - which would read as "not measured".
    for c in ("n_moved", "n_flip", "n_up", "n_down"):
        summary[c] = summary[c].fillna(0).astype(int)

    summary.insert(0, "state", ctx["state"])
    summary.insert(0, "cell_line", ctx["cell_line"])
    summary.insert(0, "model", ctx["stem"])
    summary = summary.join(modules.set_index("Module_ID"), how="left")
    summary.index.name = "perturbed"
    return summary.sort_values("sum_abs_delta_final", ascending=False)


def stamp():
    """`[2026-08-07 14:32:01]`, matching the run_all_gabi_<branch>.sh log format.

    A sweep runs for hours under nohup, so the log is the only record of how long
    it took - and "how long does this take" is a question the write-up has to
    answer. Averages and an ETA don't survive as evidence; wall-clock lines do.
    """
    return time.strftime("[%Y-%m-%d %H:%M:%S]")


def elapsed_str(seconds):
    """`1h 42m` / `9m 30s` - a duration to quote, not a float of seconds."""
    if seconds >= 3600:
        return f"{int(seconds // 3600)}h {int(seconds % 3600 // 60)}m"
    if seconds >= 60:
        return f"{int(seconds // 60)}m {int(seconds % 60)}s"
    return f"{seconds:.0f}s"


def sweep(bnet, args):
    """Sweep every module of one model."""
    stem = os.path.splitext(os.path.basename(bnet))[0]
    cell_line = stem.split("_")[0].upper()
    cell_line = cell_line if cell_line in CELL_LINES else ""

    map_path = find_name_map(bnet, args.name_map)
    nmap = load_name_map(map_path) if map_path else {}
    rev = {orig: safe for safe, orig in nmap.items()}

    model_started = time.time()
    print(f"\n{stamp()} === {bnet} ===")
    print(f"name map: {map_path or 'none found - modules stay as MaBoSS-safe names'}")

    ## The wild-type arm, once. Every perturbed arm is compared against this one
    ## run, which is the whole reason a sweep is affordable.
    t0 = time.time()
    wt_sim, background, freed = build_sim(bnet, args, args.mutate, rev, nmap, args.seed)
    wt = node_traj(wt_sim.run(), nmap)
    for disp, state in background:
        print(f"background mutation (all arms): {disp} {state}")
    if args.free_inputs:
        print(f"free inputs (istate 50/50): "
              f"{', '.join(nmap.get(n, n) for n in freed) or 'none found'}")

    ## Noise floor from a second seed, since MaBoSS is deterministic per seed - the
    ## same reasoning as compare_maboss.py --noise-floor, but always on here: one
    ## extra run buys the `above_noise` flag on every row of a 40,000-row table.
    alt_sim, _, _ = build_sim(bnet, args, args.mutate, rev, nmap, (args.seed or 0) + 1)
    floor = float((node_traj(alt_sim.run(), nmap) - wt).abs().max().max())
    print(f"wild-type + noise floor: {time.time() - t0:.0f}s, "
          f"floor (two seeds) max |Δ| = {floor:.4f}")

    annot_path = find_display_csv(stem, args.display_csv)
    annot = load_annotation(annot_path) if annot_path else None
    if annot is not None:
        ## Same reasoning as ID_COLS: the display CSV is read by compare_maboss.py's
        ## loader, which infers, so pin the join key to text here rather than find
        ## out through a table of silently blank labels.
        annot.index = annot.index.astype(str)
    print(f"module labels: {annot_path or 'none found - table goes out unannotated'}")

    indeg, outdeg = degrees(bnet)
    forward = forward_edges(load_rules(bnet))
    modules = module_table(wt, nmap, annot, indeg, outdeg, freed)

    os.makedirs(args.outdir, exist_ok=True)
    mod_path = os.path.join(args.outdir, f"{stem}_sweep_modules.csv")
    modules.to_csv(mod_path, index=False)
    print(f"-> {mod_path}  ({len(modules)} modules)")

    ## What to perturb: every node of the model except those the background already
    ## holds fixed, since perturbing those is not an experiment.
    held = {n for n, _ in background}
    targets = [(safe, nmap.get(safe, safe)) for safe in wt_sim.network
               if nmap.get(safe, safe) not in held and safe not in held]
    if args.only:
        wanted = set(args.only)
        targets = [(s, m) for s, m in targets if m in wanted or s in wanted]
        missing = wanted - {m for _, m in targets} - {s for s, _ in targets}
        if missing:
            print(f"--only: not in this model, ignored: {', '.join(sorted(missing))}")
    targets.sort(key=lambda sm: sm[1])

    delta_path = os.path.join(args.outdir, f"{stem}_sweep_{args.state}.csv")
    done = set()
    if args.resume and os.path.isfile(delta_path):
        done = set(pd.read_csv(delta_path, usecols=["perturbed"],
                               dtype=str)["perturbed"])
        print(f"--resume: {len(done)} perturbations already in {delta_path}")
    todo = [(s, m) for s, m in targets if m not in done]

    if not todo:
        print("nothing left to run")
    else:
        print(f"{stamp()} sweeping {len(todo)} modules {args.state} "
              f"on {args.workers} worker(s)")

    ctx = {"stem": stem, "cell_line": cell_line, "state": args.state, "rev": rev,
           "annot": annot, "floor": floor,
           "direction": +1 if args.state == "ON" else -1}

    ## Header only on a genuinely new file, so --resume appends to what is there.
    header = not (args.resume and os.path.isfile(delta_path))
    if header and os.path.isfile(delta_path):
        os.remove(delta_path)

    traj_path = (os.path.join(args.outdir, f"{stem}_sweep_{args.state}_traj.npz")
                 if args.trajectories else None)
    traj = {}
    if traj_path and args.resume and os.path.isfile(traj_path):
        traj = load_trajectories(traj_path)
        print(f"--resume: {len(traj)} trajectories already in {traj_path}")

    ## Results are buffered and the two files written together, so the CSV never
    ## gets ahead of the .npz. If it did, a crash would leave perturbations that
    ## --resume skips (they are in the table) but whose curves were never saved,
    ## and the only way back would be re-running them. Within a checkpoint the
    ## .npz goes first: extra trajectories are harmless, missing ones are not.
    pending = []

    def checkpoint():
        nonlocal header
        if traj_path:
            save_trajectories(traj_path, wt, traj)
        if pending:
            pd.concat(pending).to_csv(delta_path, mode="a", header=header, index=False)
            header = False
            pending.clear()

    failures, started = [], time.time()
    ## chunksize 1: the runs are ~1 minute each and vary, so handing them out one at
    ## a time keeps the workers evenly loaded and the progress line truthful.
    with mp.Pool(args.workers, _init_worker, (bnet, args, rev, nmap, wt)) as pool:
        for i, (safe, result) in enumerate(
                pool.imap_unordered(_run_perturbation, [s for s, _ in todo], 1), 1):
            module = nmap.get(safe, safe)
            if isinstance(result, str):
                print(f"  [{i}/{len(todo)}] {module}: {result}")
                failures.append(module)
                continue
            summary, grid = result

            dist, sign, _ = breadth_first(forward, safe)
            rows = annotate(summary, module, dist, sign, ctx)
            pending.append(rows)
            if grid is not None:
                traj[module] = grid
            if i % args.checkpoint == 0:
                checkpoint()

            moved = int((rows["above_noise"] & ~rows["is_target"]).sum())
            elapsed = time.time() - started
            eta = elapsed / i * (len(todo) - i)
            print(f"{stamp()}   [{i}/{len(todo)}] {module}: {moved} modules moved, "
                  f"{len(dist) - 1} reachable  ({elapsed / i:.0f}s/run, "
                  f"ETA {elapsed_str(eta)})", flush=True)
    checkpoint()
    if traj_path:
        print(f"-> {traj_path}  ({len(traj)} trajectories, "
              f"{os.path.getsize(traj_path) / 1e6:.1f} MB)")

    if os.path.isfile(delta_path):
        print(f"-> {delta_path}")
        sum_path = os.path.join(args.outdir, f"{stem}_sweep_{args.state}_summary.csv")
        summarise(delta_path, modules, ctx).to_csv(sum_path)
        print(f"-> {sum_path}")

    if failures:
        print(f"{len(failures)} perturbation(s) failed: {', '.join(failures)}")
    print(f"{stamp()} {stem} done in {elapsed_str(time.time() - model_started)}"
          f" ({len(todo)} perturbations)")
    return not failures


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Perturb every module of a MaBoSS .bnet model in turn and write "
                    "the per-module effects out as one long, queryable table.")
    parser.add_argument("bnet", nargs="+", help="MaBoSS .bnet model file(s)")
    parser.add_argument("--state", choices=("OFF", "ON"), default="OFF",
                        help="how each module is perturbed in its own arm: OFF = "
                             "knock-out, ON = constitutive (default OFF)")
    parser.add_argument("--workers", type=int, default=None, metavar="N",
                        help=f"perturbed arms to run in parallel (default: cores, "
                             f"capped at one worker per {WORKER_GB} GB of free "
                             f"memory - reading a run's output back is memory-hungry "
                             f"enough to lock up a box that only counts cores). "
                             f"Independent arms parallelise better than MaBoSS's own "
                             f"--threads, so prefer workers over threads")
    parser.add_argument("--only", action="append", default=[], metavar="MODULE",
                        help="perturb only this module (repeatable) instead of all "
                             "of them - for a smoke test before committing the hours")
    parser.add_argument("--resume", action="store_true",
                        help="keep the existing table and run only the modules "
                             "missing from it")
    parser.add_argument("--no-traj", dest="trajectories", action="store_false",
                        help="skip the <stem>_sweep_<STATE>_traj.npz of P(ON) curves. "
                             "Only worth it if disk is short: the curves are a few MB "
                             "and are the only way to plot a perturbation over time "
                             "afterwards without re-running the simulation")
    parser.add_argument("--checkpoint", type=int, default=25, metavar="N",
                        help="flush results to disk every N perturbations (default 25) "
                             "- a crash costs at most this many runs")
    parser.add_argument("--outdir", default=DEFAULT_SWEEP_DIR, metavar="DIR",
                        help=f"where the three CSVs go (default: {DEFAULT_SWEEP_DIR})")
    parser.add_argument("--name-map", default=None,
                        help="6-cyto_to_neko.py <stem>_name_map.tsv, applied to every "
                             "model; omit to use each model's own map beside it")
    parser.add_argument("--display-csv", default=None, metavar="FILE",
                        help="curated module labels; omit to use "
                             "data/cyto_exp/<CL>_display.csv for the model's cell line")
    parser.add_argument("--mutate", action="append", default=[], metavar="NODE=ON|OFF",
                        help="background applied to EVERY arm including wild-type, so "
                             "the whole sweep runs on a mutant background; repeatable")
    parser.add_argument("--free-inputs", action="store_true",
                        help="start each self-input module at 50/50 to sample the "
                             "whole landscape")
    parser.add_argument("--seed", type=int, default=None, metavar="N",
                        help="MaBoSS seed_pseudorandom, applied to every arm")
    parser.add_argument("--max-time", type=float, default=None,
                        help="MaBoSS max_time (default: leave the model's own value)")
    parser.add_argument("--sample-count", type=int, default=None,
                        help="MaBoSS sample_count / number of trajectories")
    parser.add_argument("--threads", type=int, default=1, metavar="N",
                        help="MaBoSS thread_count WITHIN each arm (default 1 - the "
                             "cores are better spent on --workers)")
    args = parser.parse_args(argv)

    ## A sweep runs for hours and is normally watched by tailing a redirected log,
    ## where Python's default block buffering would leave that log empty for the
    ## first hour. Line buffering makes progress visible as it happens.
    sys.stdout.reconfigure(line_buffering=True)

    if args.workers is None:
        args.workers = default_workers()
        print(f"--workers not given: using {args.workers} of {os.cpu_count()} cores "
              f"(memory-capped at ~{WORKER_GB} GB per worker)")

    batch_started = time.time()
    print(f"{stamp()} BATCH START - {len(args.bnet)} model(s), state {args.state}")

    ok, failed = 0, []
    for bnet in args.bnet:
        if not os.path.isfile(bnet):
            print(f"{bnet}: FAILED - not a file")
            failed.append(bnet)
            continue
        try:
            if sweep(bnet, args):
                ok += 1
            else:
                failed.append(bnet)
        except Exception as e:                        # one bad model, not the batch
            print(f"{bnet}: FAILED - {e}")
            failed.append(bnet)

    print(f"\n{stamp()} BATCH END - {ok} model(s) swept, {len(failed)} with "
          f"failures, total {elapsed_str(time.time() - batch_started)}")
    if failed:
        print("Failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
