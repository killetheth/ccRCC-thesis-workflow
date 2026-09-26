#!/usr/bin/env python3
## Purpose: Build a MaBoSS `.bnet` model from a NeKo-importable signed SIF, using
## NeKo to turn the topology + edge signs into Boolean rules. The middle step of
## the Cytoscape -> NeKo -> MaBoSS route:
##
##   6-cyto_to_neko.py  ->  <stem>_neko.sif   (signed, MaBoSS-safe names)
##   7-sif_to_bnet.py   ->  <stem>.bnet       (this script)
##   8-run_maboss.py  ->  final-state table
##
## For each SIF: load it into a NeKo Network (OmniPath is named only to satisfy
## the engine - the topology comes entirely from the SIF), strip any bimodal
## edges (each one would otherwise double the number of exported files), and write
## one clean `.bnet`. The `_neko` suffix 6-cyto_to_neko.py adds is dropped from the
## output stem, so `<stem>.bnet` pairs with that network's `<stem>_name_map.tsv`.
##
## Give the SIF file(s) as arguments to convert only the networks you chose; with
## no arguments it falls back to every `*.sif` in the current folder (alphabetical).
##
## Outputs go to gabi-maboss/output/for_maboss/ - the folder holding what MaBoSS
## is fed, as against gabi-maboss/output/maboss/ where 8-run_maboss.py writes what
## MaBoSS produces. It is resolved from this script's own location, so it holds
## wherever you run from; -o OUTDIR overrides. Keeping the `.bnet` beside its
## `<stem>_name_map.tsv` there (6-cyto_to_neko.py writes the map to that same
## folder) is what lets 8-run_maboss.py find the map on its own.
##
## Usage:
##   python 7-sif_to_bnet.py [SIF ...] [-o OUTDIR]
## Run from wherever the SIFs are (gabi-maboss/data/for_neko/ by default), or give
## full paths. Needs the `neko` Python package.

import argparse
import glob
import os
import sys

from neko.core.network import Network      # loads each SIF into a Network object
from neko._outputs.exports import Exports  # the BNET/SIF writer


## Canonical home for what MaBoSS is fed (the `.bnet` and its name map), relative
## to this script (gabi-maboss/scripts/) rather than the caller's cwd - so `-o` is
## only needed when you deliberately want the output somewhere else.
DEFAULT_OUTDIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "for_maboss"))


def bnet_stem(sif_path, outdir):
    """Output stem for one SIF: basename minus `.sif` and a trailing `_neko`,
    placed in outdir."""
    stem = os.path.basename(sif_path).removesuffix(".sif").removesuffix("_neko")
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, stem)


def convert_one(sif, outdir):
    """Convert a single SIF to a `.bnet`. Returns True on success."""
    net = Network(sif_file=sif, resources="omnipath")       # topology + signs from the SIF

    ## Bimodal edges make NeKo emit two BNETs per model; strip them for one clean file.
    bimodal = net.edges[net.edges["Effect"] == "bimodal"]
    if not bimodal.empty:
        print(f"{sif}: {len(bimodal)} bimodal edge(s) - removing")
        net.remove_bimodal_interactions()

    stem = bnet_stem(sif, outdir)
    ## glob.escape() so a stem containing `[` or `]` (legal in a Cytoscape network
    ## name, and a character class to glob) still matches only its own files.
    neko_files = glob.escape(stem) + "_[0-9]*.bnet"

    ## Clear this stem's NeKo intermediates from an earlier run before exporting.
    ## Without this, a leftover `<stem>_2.bnet` makes the post-export glob find two
    ## files, the rename below is skipped, and no clean <stem>.bnet is produced.
    ## Only `<stem>_<n>.bnet` is touched - never another network's files, and never
    ## <stem>.bnet itself, which the rename overwrites at the end anyway.
    stale = sorted(glob.glob(neko_files))
    for path in stale:
        os.remove(path)
    if stale:
        print(f"{sif}: cleared {len(stale)} NeKo intermediate(s) from a previous run")

    Exports(net).export_bnet(stem)                          # NeKo writes <stem>_<n>.bnet

    ## NeKo always suffixes the file index (`_1`, `_2`, …). With bimodal edges
    ## stripped there is exactly one file - rename it to the clean <stem>.bnet so
    ## it pairs with that network's <stem>_name_map.tsv.
    produced = sorted(glob.glob(neko_files))
    final = f"{stem}.bnet"
    if len(produced) == 1:
        os.replace(produced[0], final)
        print(f"{sif}  ->  {final}")
        return True
    if not produced:
        ## NeKo returned without writing anything (e.g. every edge dropped as
        ## bimodal, leaving no network to export). Report it rather than counting
        ## a run that produced no model as a success.
        print(f"{sif}: FAILED - NeKo exported no .bnet file")
        return False
    ## More than one file despite the bimodal strip: keep them and say so, since
    ## there is no single <stem>.bnet for 8-run_maboss.py to pick up.
    print(f"{sif}  ->  {', '.join(produced)}")
    print(f"    WARNING: NeKo split this model across {len(produced)} files, so no "
          f"clean {os.path.basename(final)} was written - resolve by hand before running MaBoSS.")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build MaBoSS .bnet model(s) from NeKo-importable signed SIF(s).")
    parser.add_argument("sif", nargs="*",
                        help="SIF file(s); default: every *.sif in the current folder")
    parser.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                        help=f"write the .bnet files here (default: {DEFAULT_OUTDIR})")
    args = parser.parse_args(argv)

    sif_files = args.sif if args.sif else sorted(glob.glob("*.sif"))
    if not sif_files:
        print("No .sif files given or found in this folder.")
        return 1

    converted, failed = 0, []
    for sif in sif_files:
        if not os.path.isfile(sif):
            print(f"{sif}  FAILED: not a file")
            failed.append(sif)
            continue
        try:                                                # keep going if one file fails
            if convert_one(sif, args.outdir):
                converted += 1
            else:
                failed.append(sif)
        except Exception as e:
            print(f"{sif}  FAILED: {e}")
            failed.append(sif)

    print(f"\n{converted} converted, {len(failed)} failed")
    if failed:
        print("Failed files:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
