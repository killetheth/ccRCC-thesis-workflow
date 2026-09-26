## Purpose: Shared helpers for 8-run_maboss.py (plain runs) and compare_maboss.py
## (treatment-vs-wild-type runs) - name-map lookup, mutation handling, and
## simulation setup, so both scripts build a model identically. Not runnable on
## its own; imported by both.

import csv
import os
import re

import maboss


## Canonical home for what MaBoSS produces, relative to this script's own
## location (gabi-maboss/scripts/) rather than the caller's cwd - so a bare
## --csv/--plot lands there wherever you run from.
DEFAULT_OUTDIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "maboss"))


def load_name_map(path):
    """safe_name -> original module ID, from a 6-cyto_to_neko.py name-map TSV."""
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["safe_name"]: r["original_name"]
                for r in csv.DictReader(fh, delimiter="\t")}


def find_name_map(bnet, explicit):
    """Locate the name map for one model: an explicit --name-map if given, else
    `<stem>_name_map.tsv` beside the model. Returns the path, or None if there
    isn't one (the model then reports its safe node names as-is).

    Auto-discovery is per model, so a multi-model run labels each network with its
    own map instead of forcing one map onto all of them. The stems line up because
    7-sif_to_bnet.py drops 6-cyto_to_neko.py's `_neko` suffix, leaving `<stem>.bnet`
    next to that same network's `<stem>_name_map.tsv`.
    """
    if explicit:
        return explicit
    beside = os.path.splitext(os.path.abspath(bnet))[0] + "_name_map.tsv"
    return beside if os.path.isfile(beside) else None


def resolve_node(token, sim, rev):
    """Map a user-given node (safe name or module ID) to its network name."""
    if token in sim.network:
        return token
    if token in rev:                                # rev: module ID -> safe name
        return rev[token]
    raise ValueError(f"unknown node '{token}' - not a model node, and not a module "
                     f"ID in the name map (pass --name-map, or use the safe name)")


def apply_mutations(sim, specs, rev, nmap):
    """Apply each `NODE=ON|OFF` mutation; return the list of (display, state)."""
    applied = []
    for spec in specs:
        node_tok, sep, state = spec.partition("=")
        state = state.strip().upper()
        if not sep or state not in ("ON", "OFF"):
            raise ValueError(f"--mutate {spec!r}: expected NODE=ON or NODE=OFF")
        node = resolve_node(node_tok.strip(), sim, rev)
        sim.mutate(node, state)                     # OFF = knock-out, ON = constitutive
        applied.append((nmap.get(node, node), state))
    return applied


def free_self_inputs(sim):
    """Set every non-mutated self-input node (logic == its own name) to istate
    50/50, so the run samples both ON and OFF starts. Returns the names freed."""
    freed = []
    for name in sim.network:
        node = sim.network[name]
        if node.logExp.strip() == name and not node.is_mutant:
            sim.network.set_istate(name, [0.5, 0.5])
            freed.append(name)
    return freed


def mutation_tag(mutations):
    """Filename-safe suffix for a mutated run, e.g. [('18.37','OFF')] -> _18_37OFF."""
    if not mutations:
        return ""
    return "_" + "_".join(re.sub(r"[^0-9A-Za-z]", "_", f"{d}{s}") for d, s in mutations)


def build_sim(bnet, args, mutate_specs, rev, nmap, seed=None):
    """Load one model and set it up: parameters, mutations, then free inputs.

    Factored out so the two arms of a compare_maboss.py run are built by the same
    code - a wild-type and a treated arm that differed in max_time, sample_count or
    seed would produce a delta that is partly just a parameter change. Order
    matters: mutations first, because free_self_inputs() skips nodes already
    mutated. Returns (sim, mutations, freed).
    """
    sim = maboss.loadBNet(bnet)

    ## Only touch a parameter the caller actually set, so the default run stays
    ## byte-for-byte the MaBoSS default (what the original smoke test did).
    params = {}
    if args.max_time is not None:
        params["max_time"] = args.max_time
    if args.sample_count is not None:
        params["sample_count"] = args.sample_count
    if seed is not None:
        params["seed_pseudorandom"] = seed
    if args.threads is not None:
        ## MaBoSS samples independent stochastic trajectories, so splitting
        ## sample_count over threads is near-linear and changes nothing about
        ## the result beyond how it is divided up.
        params["thread_count"] = args.threads
    if params:
        sim.update_parameters(**params)

    mutations = apply_mutations(sim, mutate_specs, rev, nmap)
    freed = free_self_inputs(sim) if args.free_inputs else []
    return sim, mutations, freed


def node_traj(result, nmap):
    """Per-module P(ON) over time, columns relabelled to module IDs."""
    return result.get_nodes_probtraj().rename(columns=lambda c: nmap.get(c, c))
