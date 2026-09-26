#!/usr/bin/env python3
## Purpose: Turn a Cytoscape edge-table CSV export into a NeKo-importable signed
## SIF, with MaBoSS-safe node names - fixing the two things a plain Cytoscape
## `.sif` export can't carry through to a Boolean model:
##
##   1. SIGN. Cytoscape's "Export Network as SIF" writes only the generic edge
##      type ("interacts with") in the middle column, so the activation/inhibition
##      `sign` that 3-sign_export_<branch>.R put on every edge is lost. NeKo reads
##      the middle SIF column through its effect map, and "interacts with" is not
##      in that map, so every edge lands as Effect "undefined" and the Boolean
##      rules come out unsigned. The edge *table* CSV, unlike the SIF, still holds
##      the `sign` column - so we read that and write the middle column as
##      `stimulation`/`inhibition` (both are keys NeKo recognises; note the CSV's
##      own word "activation" is NOT one, hence the translation below).
##
##   2. NODE NAMES. Modules are named `<Cluster_Number>.<Multi_Profiles>` (e.g.
##      `6.3_12`). A MaBoSS identifier must start with a LETTER and then hold only
##      letters/digits/underscores (its grammar is Word(alphas, alphanums+'_')),
##      so `6.3_12` is rejected downstream on two counts: the dot, and the leading
##      digit. We replace every non-word character with `_` and prefix `m` (for
##      module) whenever the result wouldn't start with a letter: `6.3_12` ->
##      `m6_3_12`. NB this is a LETTER prefix, not the underscore the ginsim tests
##      used - MaBoSS's CFG parser rejects a leading `_` too (that scheme only
##      suited GinSim). A `<stem>_name_map.tsv` is written alongside so MaBoSS
##      results (which speak safe names) can be translated back to module IDs.
##
## Input: one or more Cytoscape edge-table CSVs (File > Export > Table to File,
## on the edge table of the selection you want to model). Needs `source`,
## `target`, and a sign column (`sign` of activation/inhibition preferred; falls
## back to a numeric `sign_num` of +1/-1).
##
## Output - the two files serve different halves of the bridge, so they are
## written to the folder each belongs to (both paths resolved from this script's
## own location rather than the caller's cwd, so they hold wherever you run from):
##
##   <stem>_neko.sif       source <TAB> stimulation|inhibition <TAB> target
##                         -> gabi-maboss/data/for_neko/   (-o/--outdir)
##                         NeKo's input, consumed by 7-sif_to_bnet.py.
##
##   <stem>_name_map.tsv   safe_name <TAB> original_name
##                         -> gabi-maboss/output/for_maboss/  (-m/--map-outdir)
##                         NeKo never reads this; only 8-run_maboss.py does, to
##                         turn safe names back into module IDs. It goes where
##                         7-sif_to_bnet.py writes `<stem>.bnet` so that step 8
##                         finds `<stem>_name_map.tsv` beside the model on its own.
##
## The `_neko.sif` is what you feed NeKo, e.g.
##   from neko.core.network import Network
##   from neko._outputs.exports import Exports
##   net = Network(sif_file="<stem>_neko.sif", resources="omnipath")
##   Exports(net).export_bnet("<stem>")          # -> <stem>.bnet for MaBoSS
##
## Usage:
##   python 6-cyto_to_neko.py EDGE_CSV [EDGE_CSV ...] [-o OUTDIR] [-m MAP_OUTDIR]
## Run from anywhere - paths are taken as given. Pure standard library.

import argparse
import csv
import os
import re
import sys


## Both defaults are relative to this script (gabi-maboss/scripts/) rather than
## the caller's cwd - so the flags are only needed when you deliberately want the
## output somewhere else. The SIF is a prepared NeKo input, hence data/for_neko;
## the name map is read at the MaBoSS end, hence output/for_maboss, beside the
## `.bnet` 7-sif_to_bnet.py builds from that same SIF.
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SIF_DIR = os.path.normpath(os.path.join(_HERE, "..", "data", "for_neko"))
DEFAULT_MAP_DIR = os.path.normpath(os.path.join(_HERE, "..", "output", "for_maboss"))


def make_safe(name):
    """Sanitise one module name to a MaBoSS-legal identifier.

    Replace every non-word character with `_` (kills the `.` in
    `<cluster>.<profile>`), then prefix `m` (module) whenever the result wouldn't
    start with a letter - because MaBoSS's grammar (Word(alphas, alphanums+'_'))
    requires a leading letter, rejecting both a leading digit and a leading `_`.
    `6.3_12` -> `m6_3_12`; an existing symbol like `VHL` is left as `VHL`.
    """
    s = re.sub(r"[^0-9A-Za-z_]", "_", str(name).strip())
    if not s or not s[0].isalpha():
        s = "m" + s
    return s


## The CSV's `sign` column says activation/inhibition; NeKo's effect map keys are
## stimulation/inhibition (it does NOT recognise the word "activation"). Translate.
SIGN_TO_EFFECT = {"activation": "stimulation", "inhibition": "inhibition"}


def resolve_effect(row):
    """Return the NeKo effect string for one edge row, or None if unresolved.

    Prefer the textual `sign` column; fall back to numeric `sign_num` (+1/-1).
    Blank/NA/anything else -> None (reported as unsigned, written as `undefined`).
    """
    sign = (row.get("sign") or "").strip().lower()
    if sign in SIGN_TO_EFFECT:
        return SIGN_TO_EFFECT[sign]

    ## The numeric fallback matches an explicit whitelist and must keep doing so.
    ## An edge that could not be signed carries sign = NA in R, and igraph writes
    ## the accompanying NA_integer_ out as -2147483648, the minimum 32-bit int.
    ## Testing the number's sign instead - float(num) > 0, or similar - would
    ## therefore classify every unsignable edge as INHIBITION, which under the
    ## 'any activator unless any inhibitor' rule is the most damaging error
    ## available: one spurious inhibitor holds a node off permanently. There are
    ## 362 such edges across the 24 signed networks (103 Reactome, 251 HumanNet,
    ## 8 HuRI), so this is not hypothetical. Anything unrecognised must fall
    ## through to None and be written as `undefined`.
    num = (row.get("sign_num") or "").strip()
    if num in ("1", "1.0", "+1"):
        return "stimulation"
    if num in ("-1", "-1.0"):
        return "inhibition"
    return None


def convert_one(csv_path, sif_dir, map_dir):
    """Convert a single Cytoscape edge CSV to a signed SIF + name map.

    The SIF lands in sif_dir (NeKo's input) and the name map in map_dir (read at
    the MaBoSS end); they are separate folders by design, see the header.
    Returns True on success, False if the file was unusable (missing columns).
    """
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        for needed in ("source", "target"):
            if needed not in cols:
                print(f"{csv_path}: FAILED - no '{needed}' column (found: {', '.join(cols)})")
                return False
        if "sign" not in cols and "sign_num" not in cols:
            print(f"{csv_path}: FAILED - no 'sign' or 'sign_num' column to read the edge sign from")
            return False
        rows = list(reader)

    safe_of = {}       # original name -> safe name
    edges = []         # (safe_source, effect, safe_target)
    n_stim = n_inhib = n_undef = 0

    for row in rows:
        src_raw, tgt_raw = row["source"].strip(), row["target"].strip()
        if not src_raw or not tgt_raw:
            continue                                   # skip isolated-node / blank rows
        for raw in (src_raw, tgt_raw):
            safe_of.setdefault(raw, make_safe(raw))

        effect = resolve_effect(row)
        if effect == "stimulation":
            n_stim += 1
        elif effect == "inhibition":
            n_inhib += 1
        else:
            effect = "undefined"                       # keep the edge, flag it below
            n_undef += 1
        edges.append((safe_of[src_raw], effect, safe_of[tgt_raw]))

    ## A safe name that stands for two different originals would silently merge
    ## two modules in the model - refuse rather than corrupt the topology.
    collisions = {}
    for original, safe in safe_of.items():
        collisions.setdefault(safe, []).append(original)
    clashing = {s: o for s, o in collisions.items() if len(o) > 1}
    if clashing:
        print(f"{csv_path}: FAILED - name collisions after sanitising:")
        for safe, originals in clashing.items():
            print(f"    {safe}  <-  {', '.join(sorted(originals))}")
        return False

    stem = os.path.splitext(os.path.basename(csv_path))[0]
    os.makedirs(sif_dir, exist_ok=True)
    os.makedirs(map_dir, exist_ok=True)
    sif_path = os.path.join(sif_dir, f"{stem}_neko.sif")
    map_path = os.path.join(map_dir, f"{stem}_name_map.tsv")

    with open(sif_path, "w", encoding="utf-8") as fh:
        for s, effect, t in edges:
            fh.write(f"{s}\t{effect}\t{t}\n")

    with open(map_path, "w", encoding="utf-8") as fh:
        fh.write("safe_name\toriginal_name\n")
        for original, safe in sorted(safe_of.items(), key=lambda kv: kv[1]):
            fh.write(f"{safe}\t{original}\n")

    print(f"{csv_path}  ->  {sif_path}")
    print(f"    name map  ->  {map_path}")
    print(f"    {len(edges)} edges, {len(safe_of)} nodes | "
          f"{n_stim} stimulation, {n_inhib} inhibition, {n_undef} undefined")
    if n_undef:
        print(f"    WARNING: {n_undef} edge(s) had no usable sign - written as "
              f"'undefined'; resolve these by hand before running MaBoSS.")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert Cytoscape edge-table CSV exports into NeKo-importable "
                    "signed SIFs with MaBoSS-safe node names.")
    parser.add_argument("csv", nargs="+", help="Cytoscape edge-table CSV export(s)")
    parser.add_argument("-o", "--outdir", default=DEFAULT_SIF_DIR,
                        help=f"write the _neko.sif here (default: {DEFAULT_SIF_DIR})")
    parser.add_argument("-m", "--map-outdir", default=DEFAULT_MAP_DIR,
                        help="write the _name_map.tsv here - keep it wherever "
                             "7-sif_to_bnet.py puts the .bnet, so 8-run_maboss.py "
                             f"finds it on its own (default: {DEFAULT_MAP_DIR})")
    args = parser.parse_args(argv)

    converted, failed = 0, []
    for path in args.csv:
        if not os.path.isfile(path):
            print(f"{path}: FAILED - not a file")
            failed.append(path)
            continue
        if convert_one(path, args.outdir, args.map_outdir):
            converted += 1
        else:
            failed.append(path)

    print(f"\n{converted} converted, {len(failed)} failed")
    if failed:
        print("Failed files:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
