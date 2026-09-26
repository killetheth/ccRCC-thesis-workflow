#!/usr/bin/env python3
## Purpose: Turn sweep_maboss.py's whole-model knockout database into a process-level
## report - specifically "which module knockouts move APOPTOSIS, and which move
## PROLIFERATION, in each network". The sweep tables are one row per
## (knocked-out module -> readout module) pair; this script adds the process
## annotation those rows lack, fixes the significance flag they carry, aggregates
## the pairs up to process level, and draws the figures for the strongest hits.
##
## Reads (from gabi-maboss/output/maboss/sweep/):
##   <stem>_sweep_modules.csv    - one row per module: labels, P_ON_wt, degree
##   <stem>_sweep_OFF.csv        - the pair database
##   <stem>_sweep_OFF_traj.npz   - every arm's P(ON) curves (for the figures)
##
## Writes:
##   <query outdir>/  the derived tables (see TABLES below)
##   <figure outdir>/ paired trajectory plots + delta heatmaps, styled to match
##                    compare_maboss.py's save_compare_plots() so a figure from
##                    here and one from a targeted run are directly comparable
##
## Three things this script does that the raw sweep does not:
##
## 1. PROCESS ANNOTATION. The sweep carries `Cluster Function`, but that
##    vocabulary has no apoptosis category at all - apoptotic modules sit inside
##    the catch-all "Regulation (metabolism, development, apoptosis)" - so
##    apoptosis has to be recovered from the gene lists in `Module_HGNCs`. Both a
##    gene-level and a Cluster-Function-level call are computed, and both are
##    reported, because they disagree and the disagreement is informative.
##
## 2. SIGNIFICANCE RECALIBRATION. The sweep's own `above_noise` flag is far too
##    permissive: its floor comes from re-running wild-type at a second seed, and
##    with --free-inputs only a handful of nodes are stochastic, so two wild-type
##    runs agree to ~1e-4 while *mutating* any node shifts every module by ~1e-3.
##    The result is knockouts credited with moving 184 of 186 modules while
##    reaching none of them. The fix used here needs no re-running: for each
##    knockout, the modules with no directed path from it (blank `hop`) provably
##    cannot have been affected, so THEIR delta spread is that knockout's own
##    noise - an empirical null measured on the real thing, per perturbation.
##
## 3. AGGREGATION. Pair rows are rolled up to (knockout x process) so a knockout
##    can be ranked by what it does to apoptosis as a whole, and to proliferation
##    as a whole, rather than to one module at a time.
##
## Usage:
##   python sweep_process_report.py [--sweep-dir DIR] [--fig-dir DIR]
##                                  [--query-dir DIR] [--top N] [--no-figures]
## Pure pandas/numpy/matplotlib - no MaBoSS engine and no `maboss` package, so
## unlike compare_maboss.py this runs outside the `maboss` conda env.

import argparse
import os
import sys

import numpy as np
import pandas as pd


## ---------------------------------------------------------------------------
## Paths. Resolved from this script's own location rather than the caller's cwd,
## matching the convention every other script in the repo uses (setwd() in R,
## __file__ in the Python pipeline steps), so the defaults hold whichever folder
## the script is run from.
## ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))              # gabi-maboss/scripts/query/
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))  # repository root
SWEEP_DIR = os.path.join(REPO, "gabi-maboss", "output", "maboss", "sweep")
QUERY_DIR = os.path.join(SWEEP_DIR, "process_query")           # derived tables (gitignored)
FIG_DIR = os.path.join(REPO, "dissertation", "figures", "sweep_process")
## Where gabi-maboss/scripts/annotate_go.R writes its output. Optional: if the
## folder is absent the hand-curated gene sets below are used instead, so this
## script runs either before or after the GO re-annotation.
GO_DIR = os.path.join(REPO, "gabi-maboss", "data", "go")

## The three models the sweep covers. Each is (cell-line code, file stem). The
## stems differ in more than the code - AH's network had unsigned edges dropped -
## so they are listed rather than templated.
MODELS = [
    ("AH", "AH_reactome_bi_subnet_drop_unsigned_edge"),
    ("CH", "CH_reactome_bi_subnet_edge"),
    ("UH", "UH_reactome_bi_subnet_edge"),
]

## Cell-line codes -> the line and condition they stand for, for figure titles.
CELL_LINE_NAMES = {
    "AH": "A498, hypoxia",
    "CH": "CAKI-1, hypoxia",
    "UH": "UMRC2, hypoxia",
}


## ---------------------------------------------------------------------------
## Process gene sets.
##
## Curated by hand against the 650-gene universe these three networks actually
## contain, so every symbol below is present in at least one module - the sets
## are deliberately not general-purpose gene sets copied wholesale, because a
## set whose members are absent inflates nothing but the methods section.
##
## Membership follows the canonical apoptotic machinery (intrinsic and extrinsic
## arms plus their direct BCL-2-family and IAP regulators) and the canonical
## cell-cycle/mitotic machinery. Two tiers each: CORE is the machinery itself and
## carries the headline analysis; EXTENDED adds upstream regulators and
## stress-response effectors whose apoptotic or proliferative role is real but
## context-dependent, and is reported separately so the headline never rests on
## a debatable call.
##
## CAVEAT worth carrying into the write-up: these are hand-curated, not drawn
## from a versioned resource (the sandbox has no network access to GO/MSigDB),
## so they are auditable but not automatically reproducible. Both tiers are
## written out to the query folder so the exact membership used is on record.
## ---------------------------------------------------------------------------

## Apoptosis, tier 1: the apoptotic machinery itself.
##   APAF1        apoptosome scaffold, intrinsic pathway
##   BAX          pro-apoptotic BCL-2 effector, MOMP
##   BCL2         anti-apoptotic BCL-2 family
##   BCL2L1       BCL-xL, anti-apoptotic BCL-2 family
##   BIRC5        survivin, IAP family (also a mitotic CPC component - see note)
##   CASP4        inflammatory/ER-stress caspase
##   CASP6        effector caspase
##   FAS          death receptor, extrinsic pathway
##   PMAIP1       NOXA, BH3-only pro-apoptotic
##   TNFSF10      TRAIL, death-receptor ligand
APOPTOSIS_CORE = {
    "APAF1", "BAX", "BCL2", "BCL2L1", "BIRC5", "CASP4", "CASP6",
    "FAS", "PMAIP1", "TNFSF10",
}

## Apoptosis, tier 2: regulators and effectors that gate apoptosis without being
## the machinery. DDIT3 (CHOP) and TRIB3 are the ER-stress arm; NFKBIA/RELA the
## NF-kB survival axis; JUN/ATF3 stress-activated transcription; IGFBP3/BTG2
## p53-responsive pro-apoptotic/anti-proliferative targets; PRKCD and RIPK2
## apoptotic and necroptotic signalling; HMOX1 and NQO1 the oxidative-stress
## arm; MYC sensitises to apoptosis as well as driving cycle entry.
APOPTOSIS_EXTENDED = APOPTOSIS_CORE | {
    "DDIT3", "TRIB3", "NFKBIA", "RELA", "JUN", "ATF3", "IGFBP3",
    "BTG2", "PRKCD", "RIPK2", "HMOX1", "NQO1", "MYC", "HMGB1",
}

## Proliferation, tier 1: cell-cycle drivers, the mitotic apparatus and the
## kinetochore - the modules whose activity IS proliferative capacity.
PROLIFERATION_CORE = {
    # cyclins, CDKs and their direct regulators
    "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCND1", "CCNE2", "CDK6",
    "CDKN1B", "CDKN2C", "CDKN3", "CDC25A", "CDC25B", "CKS1B", "CKS2",
    "SKP2", "WEE1", "CHEK1", "E2F7", "TFDP1", "CDC7", "DBF4",
    # anaphase-promoting complex and mitotic exit
    "CDC20", "CDC27", "PTTG1", "FBXO5", "GMNN",
    # mitotic kinases, kinesins and spindle
    "AURKB", "MELK", "NEK2", "TTK", "PBK", "TPX2", "TACC3", "STMN1",
    "KIF11", "KIF14", "KIF15", "KIF23", "KIF2C", "KIFC1", "PRR11",
    "ANLN", "ECT2", "ASPM", "HMMR", "TROAP", "PSRC1", "CKAP2L",
    "GAS2L3", "STIL", "NCAPD2", "NCAPG2", "CDCA3", "CDCA7", "CDCA8",
    # kinetochore and centromere
    "BUB1", "MAD2L1", "NDC80", "KNTC1", "MIS12", "MIS18BP1", "ZWINT",
    "SPC25", "CENPA", "CENPB", "CENPF", "CENPK", "CENPL", "CENPN",
    "CENPW", "CENPX", "HJURP", "TOP2A",
}

## Proliferation, tier 2: DNA replication licensing and synthesis. Real
## proliferation readouts, but they are also the substance of the networks'
## "DNA metabolic process" cluster, so folding them into the headline set would
## make roughly a third of every network "proliferation" and flatten the
## contrast the report is built on.
PROLIFERATION_EXTENDED = PROLIFERATION_CORE | {
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MCM10",
    "GINS2", "GINS3", "ORC5", "ORC6", "CDC45", "POLA1", "POLA2",
    "POLD2", "POLD4", "POLE2", "POLE3", "PRIM1", "RFC2", "RFC4",
    "RFC5", "RPA2", "FEN1", "DTL", "DSCC1", "TYMS", "DHFR", "DUT",
    "ASF1B", "CHAF1B", "UHRF1", "NASP", "TIMELESS", "TIPIN", "EXO1",
    "MYC", "PCNA", "RRM2",
}

## The Cluster Function label that comes closest to each process in the
## networks' own vocabulary, used as the independent second call. Note the
## asymmetry: proliferation has its own category, apoptosis does not - it is
## bundled with metabolism and development in one label, which is exactly why
## the gene-level call carries the headline.
CF_PROLIFERATION = "Proliferation"
CF_APOPTOSIS_MIXED = "Regulation (metabolism, development, apoptosis)"

## Percentile of the unreachable-module null taken as the significance floor.
## 99 rather than the maximum: the maximum of a ~120-value null is itself noisy
## and a single outlier would set the threshold for the whole perturbation.
NULL_PERCENTILE = 99.0


## ---------------------------------------------------------------------------
## Loading
## ---------------------------------------------------------------------------

def load_modules(stem, sweep_dir):
    """One row per module: labels, wild-type P(ON), in/out degree.

    Module IDs are `<Cluster_Number>.<Multi_Profiles>`, so an all-numeric ID
    like `18.30` reads as a float and becomes `18.3`, matching nothing on the
    join back to the pair table. dtype=str on the ID columns is not optional.
    """
    path = os.path.join(sweep_dir, f"{stem}_sweep_modules.csv")
    df = pd.read_csv(path, dtype={"Module_ID": str, "safe_name": str, "Display": str})
    df["Module_HGNCs"] = df["Module_HGNCs"].fillna("")          # modules with no gene list
    df["Cluster Function"] = df["Cluster Function"].fillna("")  # unannotated clusters
    return df


def load_pairs(stem, sweep_dir):
    """The pair database: one row per (knocked-out module, readout module).

    Same string-ID discipline as above, on both ends of the pair. `hop` is left
    as a float because it is genuinely missing for unreachable pairs, and NaN is
    the honest representation of "no directed path" - that missingness is the
    signal the noise recalibration below is built on.
    """
    path = os.path.join(sweep_dir, f"{stem}_sweep_OFF.csv")
    id_cols = {"perturbed": str, "readout": str,
               "perturbed_display": str, "readout_display": str}
    df = pd.read_csv(path, dtype=id_cols, low_memory=False)
    return df


def load_traj(stem, sweep_dir):
    """The P(ON) curves: wild-type, every perturbed arm, and the label arrays.

    Returns (wt, traj, perturbed, modules, times) with wt as time x module and
    traj as perturbation x time x module. This is what makes the figures
    possible without re-simulating anything.
    """
    path = os.path.join(sweep_dir, f"{stem}_sweep_OFF_traj.npz")
    if not os.path.isfile(path):
        return None
    z = np.load(path, allow_pickle=True)
    return {
        "wt": z["wt"],                                  # (n_times, n_modules)
        "traj": z["traj"],                              # (n_perturbed, n_times, n_modules)
        "perturbed": [str(x) for x in z["perturbed"]],  # row order of traj
        "modules": [str(x) for x in z["modules"]],      # column order of wt/traj
        "times": z["times"],
    }


## ---------------------------------------------------------------------------
## Process annotation
## ---------------------------------------------------------------------------

def split_genes(cell):
    """`Module_HGNCs` is a comma-joined gene list; return it as a set."""
    return {g.strip() for g in str(cell).split(",") if g.strip()}


def load_go_annotation(go_dir, cl):
    """Load annotate_go.R's output for one cell line, or (None, None).

    Returns (gene_flags, module_labels):
      gene_flags    HGNC symbol -> set of process names it belongs to, taken
                    from GO membership with ancestors included
      module_labels Module_ID -> (label, label_source, n_annotated) - the
                    GO-derived replacement for the legacy Cluster Function

    Returning None rather than raising when the folder is absent is deliberate:
    the report has to be runnable before the R re-annotation has been done, and
    the fallback (hand-curated sets) is stated in the output either way.
    """
    flags_path = os.path.join(go_dir, "gene_process_flags.tsv")
    labels_path = os.path.join(go_dir, "all_module_go_labels.tsv")

    gene_flags = None
    if os.path.isfile(flags_path):
        f = pd.read_csv(flags_path, sep="\t", dtype={"HGNC": str, "ENTREZID": str})
        ## Every non-bookkeeping boolean column is a process. Reading them out of
        ## the file rather than hardcoding the two names means adding a process
        ## to PROCESS_ROOTS in the R script needs no change here.
        procs = [c for c in f.columns
                 if c not in ("HGNC", "ENTREZID") and not c.endswith("_terms")
                 and f[c].dtype == bool]
        gene_flags = {}
        for _, r in f.iterrows():
            hit = {p for p in procs if bool(r[p])}
            if hit:
                gene_flags[str(r["HGNC"]).strip()] = hit

    module_labels = None
    if os.path.isfile(labels_path):
        L = pd.read_csv(labels_path, sep="\t",
                        dtype={"Module_ID": str, "Cell_Line": str})
        L = L[L["Cell_Line"] == cl]
        module_labels = {
            str(r["Module_ID"]): (r.get("label"), r.get("label_source"),
                                  r.get("n_annotated"))
            for _, r in L.iterrows()
        }
    return gene_flags, module_labels


def annotate_modules(mods, gene_flags=None, module_labels=None):
    """Tag every module with its apoptosis / proliferation membership.

    Two independent annotation systems are computed and both are kept:

      *_curated   from the hand-curated gene sets in this file
      *_go        from GO membership, if annotate_go.R has been run
      *_cf        the networks' own legacy Cluster Function agrees

    `*_core` - the column every downstream table and figure uses - is the GO
    call when GO is available and the curated call otherwise. Keeping the
    superseded call in the output rather than overwriting it is what lets the
    two be diffed, which matters here: the legacy labels came from a ~2013 GO
    snapshot that could not annotate most of the genes in question, so the
    disagreement is itself a result.

    `*_extended` is the curated sensitivity tier, unchanged by GO.
    """
    out = mods.copy()
    genes = out["Module_HGNCs"].map(split_genes)

    ## Which curated genes each module actually contributed, kept as text so the
    ## published tables are self-explaining rather than a column of True/False.
    out["apoptosis_genes"] = genes.map(
        lambda g: ", ".join(sorted(g & APOPTOSIS_CORE)))
    out["apoptosis_genes_ext"] = genes.map(
        lambda g: ", ".join(sorted(g & (APOPTOSIS_EXTENDED - APOPTOSIS_CORE))))
    out["proliferation_genes"] = genes.map(
        lambda g: ", ".join(sorted(g & PROLIFERATION_CORE)))
    out["proliferation_genes_ext"] = genes.map(
        lambda g: ", ".join(sorted(g & (PROLIFERATION_EXTENDED - PROLIFERATION_CORE))))

    out["apoptosis_curated"] = out["apoptosis_genes"] != ""
    out["apoptosis_extended"] = out["apoptosis_curated"] | (out["apoptosis_genes_ext"] != "")
    out["apoptosis_cf"] = out["Cluster Function"] == CF_APOPTOSIS_MIXED

    out["proliferation_curated"] = out["proliferation_genes"] != ""
    out["proliferation_extended"] = out["proliferation_curated"] | (out["proliferation_genes_ext"] != "")
    out["proliferation_cf"] = out["Cluster Function"] == CF_PROLIFERATION

    ## --- GO membership, when annotate_go.R has been run ------------------
    if gene_flags:
        for proc in ("apoptosis", "proliferation"):
            ## Which of the module's genes GO puts in this process, kept as text
            ## for the same reason as the curated columns: a flag nobody can
            ## trace back to a gene is a flag nobody can check.
            matched = genes.map(lambda g: ", ".join(
                sorted(x for x in g if proc in gene_flags.get(x, ()))))
            out[f"{proc}_go_genes"] = matched
            out[f"{proc}_go"] = matched != ""
            out[f"{proc}_core"] = out[f"{proc}_go"]
        out["annotation_source"] = "GO"
    else:
        for proc in ("apoptosis", "proliferation"):
            out[f"{proc}_core"] = out[f"{proc}_curated"]
        out["annotation_source"] = "curated"

    ## --- GO-derived functional label, replacing Cluster Function ---------
    if module_labels:
        out["GO_Function"] = out["Module_ID"].map(
            lambda m: (module_labels.get(m) or (None,))[0])
        out["GO_label_source"] = out["Module_ID"].map(
            lambda m: (module_labels.get(m) or (None, None))[1])
        ## The legacy label is kept under its own name so the two can be
        ## compared row by row; `Cluster Function` itself is left untouched
        ## because other scripts in the pipeline still read that header.
        out["legacy_cluster_function"] = out["Cluster Function"]

    ## BIRC5 (survivin) is both an IAP and a chromosomal-passenger-complex
    ## subunit, so a module carrying it is legitimately both processes at once.
    ## Flagged rather than arbitrated - a module that is both is a real finding
    ## about the clustering, not a bug to be resolved by dropping it from one set.
    out["dual_process"] = out["apoptosis_core"] & out["proliferation_core"]

    ## --- how confidently can this module be READ as its process? ---------
    ## A module is contracted to one Boolean node, so a knockout acts on all of
    ## its genes at once. That makes the interpretation of a result depend on
    ## what else is in the node, not only on whether a known process gene is
    ## present:
    ##
    ##   single_gene  the node IS one known process gene - the result means
    ##                exactly what it appears to mean
    ##   clean        several genes, at least one known for this process and
    ##                none known for the other - reads as the process, with the
    ##                caveat that the other genes come along
    ##   mixed        carries known genes of BOTH processes, so a change in the
    ##                node cannot be attributed to either
    ##   none         no known gene of this process
    ##
    ## `process_purity` is the share of the node's genes that are known for the
    ## process, which is the continuous version of the same question: a lone
    ## BAX in a two-gene node is a stronger readout than a lone BIRC5 in a
    ## 62-gene one.
    n_genes = genes.map(len).clip(lower=1)
    for proc in ("apoptosis", "proliferation"):
        other = "proliferation" if proc == "apoptosis" else "apoptosis"
        col = f"{proc}_go_genes" if f"{proc}_go_genes" in out else f"{proc}_genes"
        ocol = f"{other}_go_genes" if f"{other}_go_genes" in out else f"{other}_genes"
        n_hit = out[col].map(lambda s: len(split_genes(s)) if s else 0)
        out[f"{proc}_n_genes"] = n_hit
        out[f"{proc}_purity"] = (n_hit / n_genes).round(3)
        out[f"{proc}_tier"] = [
            "none" if not h else
            "mixed" if o else
            "single_gene" if tot == 1 else
            "clean"
            for h, o, tot in zip(out[col], out[ocol], n_genes)
        ]

    ## The STEM temporal profile, recovered from the module ID (which is
    ## `<Cluster_Number>.<Multi_Profiles>`). Genes merged into one node share
    ## this profile as well as the inherited functional cluster, so the merge is
    ## constrained by measured co-expression and not by the functional call
    ## alone - worth carrying, because it is the main reason a merged node is
    ## defensible at all.
    out["stem_profile"] = out["Module_ID"].map(
        lambda m: m.split(".", 1)[1] if "." in str(m) else "")
    out["cluster_number"] = out["Module_ID"].map(
        lambda m: str(m).split(".", 1)[0])
    return out


## ---------------------------------------------------------------------------
## Significance recalibration
## ---------------------------------------------------------------------------

def empirical_floors(pairs, percentile=NULL_PERCENTILE):
    """Per-knockout significance floor, measured from the unreachable modules.

    For one knocked-out module, every readout with no directed path from it
    (`hop` is NaN) cannot have been affected through the network, so whatever
    delta those rows carry is that run's sampling noise. Taking a high
    percentile of their |delta| gives a floor derived from the run itself,
    rather than extrapolated from a seed change on the wild-type arm - which is
    what the sweep's own `above_noise` does, and why it flags almost everything.

    Returns one row per knocked-out module: the floor on both delta measures,
    plus how many null rows it was computed from (a knockout that reaches most
    of the network has a small null and a correspondingly shaky floor).
    """
    p = pairs.copy()
    p["abs_final"] = p["delta_final"].abs()
    p["abs_max"] = p["delta_max"].abs()
    unreachable = p[p["hop"].isna()]

    floors = unreachable.groupby("perturbed").agg(
        floor_final=("abs_final", lambda s: float(np.percentile(s, percentile))),
        floor_max=("abs_max", lambda s: float(np.percentile(s, percentile))),
        null_n=("abs_final", "size"),
        null_median=("abs_final", "median"),
    )
    return floors


def apply_significance(pairs, floors):
    """Add reachability-aware significance columns to the pair table.

    A pair is significant only if it is both reachable (there is a path the
    effect could have travelled) and larger than that knockout's own empirical
    floor. Knockouts with too small a null to calibrate (fewer than MIN_NULL
    unreachable readouts) fall back to the pooled median floor rather than
    being dropped, so a highly connected hub is not silently excluded from the
    ranking for the very reason that makes it interesting.
    """
    MIN_NULL = 10                       # below this the percentile is meaningless
    out = pairs.merge(floors, left_on="perturbed", right_index=True, how="left")

    ## Pooled fallback, from the knockouts that DID have enough null rows.
    ok = floors[floors["null_n"] >= MIN_NULL]
    pooled_final = float(ok["floor_final"].median()) if len(ok) else 0.0
    pooled_max = float(ok["floor_max"].median()) if len(ok) else 0.0
    thin = out["null_n"].isna() | (out["null_n"] < MIN_NULL)
    out.loc[thin, "floor_final"] = pooled_final
    out.loc[thin, "floor_max"] = pooled_max
    out["floor_is_pooled"] = thin

    out["reachable"] = out["hop"].notna()
    out["sig_final"] = out["reachable"] & (out["delta_final"].abs() > out["floor_final"])
    out["sig_max"] = out["reachable"] & (out["delta_max"].abs() > out["floor_max"])
    ## The self-row (a knockout's effect on itself) is trivially the largest
    ## thing in the table and would top every ranking; excluded from aggregates.
    out["is_self"] = out["perturbed"] == out["readout"]
    return out, pooled_final, pooled_max


## ---------------------------------------------------------------------------
## Aggregation to process level
## ---------------------------------------------------------------------------

def aggregate_by_process(pairs, mods, process, membership_col):
    """Roll pair rows up to one row per knockout, for one process.

    Only readouts belonging to the process are counted, and the knockout's own
    self-row is dropped. The columns answer the three different senses of
    "this knockout affects apoptosis":
      n_sig / n_readouts   how many of the process's modules it moved at all
      sum_abs_delta_final  its total effect, which favours broad shallow effects
      max_abs_delta_final  its single strongest effect, which favours deep narrow ones
      net_delta_final      the DIRECTION, summed - the column that separates a
                           knockout that switches apoptosis on from one that
                           switches it off, and the only one that answers the
                           therapeutic question
      n_up / n_down        how that net figure was arrived at, since a net near
                           zero can mean "no effect" or "half up, half down"
    """
    readout_is_process = set(mods.loc[mods[membership_col], "Module_ID"])
    if not readout_is_process:
        return pd.DataFrame()

    sub = pairs[pairs["readout"].isin(readout_is_process) & ~pairs["is_self"]].copy()
    if sub.empty:
        return pd.DataFrame()

    g = sub.groupby("perturbed")
    agg = pd.DataFrame({
        "n_readouts": g.size(),
        "n_reachable": g["reachable"].sum(),
        "n_sig": g["sig_final"].sum(),
        "n_up": g.apply(lambda d: int(((d["delta_final"] > 0) & d["sig_final"]).sum()),
                        include_groups=False),
        "n_down": g.apply(lambda d: int(((d["delta_final"] < 0) & d["sig_final"]).sum()),
                          include_groups=False),
        ## The magnitude columns count only significant pairs, so a knockout
        ## cannot accumulate a large total out of 40 sub-threshold readouts.
        "sum_abs_delta_final": g.apply(
            lambda d: float(d.loc[d["sig_final"], "delta_final"].abs().sum()),
            include_groups=False),
        "net_delta_final": g.apply(
            lambda d: float(d.loc[d["sig_final"], "delta_final"].sum()),
            include_groups=False),
        "max_abs_delta_final": g.apply(
            lambda d: float(d.loc[d["sig_final"], "delta_final"].abs().max())
            if d["sig_final"].any() else 0.0,
            include_groups=False),
        ## The transient view. A knockout can push a process hard mid-run and
        ## have it recover by the attractor, which delta_final alone reports as
        ## "nothing happened" - and on these networks that is the common case,
        ## not the exception, so it gets its own columns rather than a footnote.
        "n_sig_max": g["sig_max"].sum(),
        "sum_abs_delta_max": g.apply(
            lambda d: float(d.loc[d["sig_max"], "delta_max"].abs().sum()),
            include_groups=False),
        "net_delta_max": g.apply(
            lambda d: float(d.loc[d["sig_max"], "delta_max"].sum()),
            include_groups=False),
        "max_abs_delta_max": g.apply(
            lambda d: float(d.loc[d["sig_max"], "delta_max"].abs().max())
            if d["sig_max"].any() else 0.0,
            include_groups=False),
        ## When the biggest swing happened, for the pair that swung most. Near
        ## the end of the run means the effect is still developing at max_time;
        ## early means it was a spike the network then absorbed.
        "t_at_max_swing": g.apply(
            lambda d: float(d.loc[d["delta_max"].abs().idxmax(), "t_delta_max"])
            if len(d) else np.nan,
            include_groups=False),
        ## Nearest readout of this process, in hops. 1 means the knockout sits
        ## directly on a module of the process; a large value means whatever it
        ## does, it does at the far end of a cascade.
        "min_hop_sig": g.apply(
            lambda d: float(d.loc[d["sig_final"], "hop"].min())
            if d["sig_final"].any() else np.nan,
            include_groups=False),
    })
    agg["process"] = process
    ## Mean effect per module of the process, which is what makes the figure
    ## comparable between cell lines whose networks hold different numbers of
    ## apoptotic modules.
    agg["mean_delta_per_module"] = agg["net_delta_final"] / max(len(readout_is_process), 1)
    return agg


def responsive_readouts(pairs, mods, process, membership_col):
    """Turn the question round: one row per PROCESS MODULE, not per knockout.

    The knockout ranking answers "what should I target"; this answers "which
    apoptotic / proliferative module is actually controllable, and by what".
    A module that no knockout can move is a sink the network holds fixed, and
    that is as much a result as a responsive one - so every process module gets
    a row, including the inert ones.

    `best_*` name the single strongest knockout for that module, which is what
    makes the table usable directly: it names the specific perturbation to run
    through compare_maboss.py next.
    """
    members = mods.loc[mods[membership_col]]
    if members.empty:
        return pd.DataFrame()

    sub = pairs[pairs["readout"].isin(set(members["Module_ID"])) & ~pairs["is_self"]]
    rows = []
    for mod_id, grp in sub.groupby("readout"):
        sig = grp[grp["sig_final"]]
        sig_t = grp[grp["sig_max"]]
        ## Strongest knockout at the attractor; falls back to the strongest
        ## transient one when nothing survives to the end, so a module that is
        ## only transiently controllable still names its controller.
        best = (sig.loc[sig["delta_final"].abs().idxmax()] if len(sig)
                else (sig_t.loc[sig_t["delta_max"].abs().idxmax()] if len(sig_t) else None))
        rows.append({
            "readout": mod_id,
            "n_knockouts_moving_it": int(len(sig)),
            "n_knockouts_moving_it_transient": int(len(sig_t)),
            "n_knockouts_reaching_it": int(grp["reachable"].sum()),
            "max_abs_delta_final": float(sig["delta_final"].abs().max()) if len(sig) else 0.0,
            "max_abs_delta_max": float(sig_t["delta_max"].abs().max()) if len(sig_t) else 0.0,
            "n_up": int((sig["delta_final"] > 0).sum()),
            "n_down": int((sig["delta_final"] < 0).sum()),
            "best_knockout": str(best["perturbed"]) if best is not None else "",
            "best_knockout_display": str(best["perturbed_display"]) if best is not None else "",
            "best_delta_final": float(best["delta_final"]) if best is not None else 0.0,
            "best_delta_max": float(best["delta_max"]) if best is not None else 0.0,
            "best_hop": float(best["hop"]) if best is not None else np.nan,
        })

    out = pd.DataFrame(rows)
    ## Report the genes that actually put the module in the process under the
    ## annotation in force - the GO matches when GO is driving, the curated ones
    ## otherwise. Showing the curated column next to a GO-selected module would
    ## sometimes be blank, which reads as "no reason given".
    gene_col = (f"{process}_go_genes" if f"{process}_go_genes" in members.columns
                else f"{process}_genes")
    keep = ["Module_ID", "Display", "Cluster Function", "Module_HGNCs", "P_ON_wt",
            "in_degree", "out_degree", gene_col,
            f"{process}_tier", f"{process}_purity", f"{process}_n_genes",
            "stem_profile"]
    keep += [c for c in ("GO_Function",) if c in members.columns]
    out = out.merge(members[keep], left_on="readout", right_on="Module_ID", how="left") \
             .drop(columns="Module_ID")
    out["process"] = process
    return out.sort_values("max_abs_delta_final", ascending=False)


def top_pairs(pairs, mods, process, membership_col, n=15):
    """The single strongest (knockout -> process module) effects, unaggregated.

    The aggregate tables can hide a very large single effect inside a knockout
    whose other readouts cancel it out, so the raw pairs are reported too -
    these are the concrete statements the write-up can make ("knocking out X
    drives module Y from P(ON)=0.5 to 0.0").
    """
    members = set(mods.loc[mods[membership_col], "Module_ID"])
    sub = pairs[pairs["readout"].isin(members) & ~pairs["is_self"] & pairs["sig_final"]].copy()
    if sub.empty:
        return pd.DataFrame()
    cols = ["perturbed", "perturbed_display", "perturbed_function", "readout",
            "readout_display", "readout_function", "P_ON_wt", "P_ON_perturbed",
            "delta_final", "delta_max", "t_delta_max", "hop", "path_sign_agrees",
            "floor_final"]
    sub = sub[[c for c in cols if c in sub.columns]]
    sub["process"] = process
    return sub.reindex(sub["delta_final"].abs().sort_values(ascending=False).index).head(n)


def attach_knockout_labels(agg, mods):
    """Put the knocked-out module's own labels back on the aggregate rows.

    Without this, a ranking is a column of bare IDs like `18.37`, which say
    nothing about what was knocked out. Module IDs are also not comparable
    between cell lines (each comes from that line's own STEM clustering), so
    the gene list is the only thing that makes a cross-network comparison
    possible at all.
    """
    keep = ["Module_ID", "Display", "Cluster Function", "Module_HGNCs",
            "P_ON_wt", "in_degree", "out_degree",
            "apoptosis_core", "proliferation_core"]
    m = mods[keep].rename(columns={
        "Display": "perturbed_display",
        "Cluster Function": "perturbed_function",
        "Module_HGNCs": "perturbed_genes",
        "P_ON_wt": "perturbed_P_ON_wt",
        "apoptosis_core": "perturbed_is_apoptosis",
        "proliferation_core": "perturbed_is_proliferation",
    })
    return agg.merge(m, left_index=True, right_on="Module_ID", how="left") \
              .rename(columns={"Module_ID": "perturbed"})


## ---------------------------------------------------------------------------
## Figures - deliberately mirroring compare_maboss.py's save_compare_plots()
## ---------------------------------------------------------------------------

def _label(mod, mods_idx):
    """Plot label for one module: `ID (Display)` when the Display differs.

    Same rule as compare_maboss.py's _label, so figures from the two scripts
    label their axes identically.
    """
    disp = mods_idx.get(mod)
    if isinstance(disp, str) and disp and disp != mod:
        return f"{mod} ({disp})"
    return str(mod)


def plot_pair(traj, stem, cl, knockout, readouts, mods_idx, outdir,
              tag, width=10.0, height=6.0):
    """Wild-type vs knocked-out trajectories for a chosen set of readouts.

    One colour per module, solid = wild-type and dashed = knocked out, so each
    pair reads as one module rather than two unrelated lines - compare_maboss.py's
    convention. The curves come straight out of the .npz; nothing is re-simulated.
    """
    import matplotlib
    matplotlib.use("Agg")                       # headless: write files, no display
    import matplotlib.pyplot as plt

    if knockout not in traj["perturbed"]:
        return None
    k = traj["perturbed"].index(knockout)
    midx = {m: i for i, m in enumerate(traj["modules"])}
    sel = [m for m in readouts if m in midx]
    if not sel:
        return None

    t = traj["times"]
    ## `height or 6.0` rather than a plain default: callers pass --fig-height
    ## straight through, and that is None unless the user set it, which would
    ## otherwise reach figsize as None and fail inside matplotlib.
    fig, ax = plt.subplots(figsize=(width, height or 6.0))
    colours = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(sel), 1)))
    for colour, mod in zip(colours, sel):
        j = midx[mod]
        ax.plot(t, traj["wt"][:, j], color=colour, lw=1.6, label=_label(mod, mods_idx))
        ax.plot(t, traj["traj"][k, :, j], color=colour, lw=1.6, ls="--")

    ax.set_xlabel("time")
    ax.set_ylabel("P(module = ON)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{cl} ({CELL_LINE_NAMES.get(cl, cl)}): {_label(knockout, mods_idx)} "
                 f"knocked OFF\n{tag}")

    ## Two legends stacked down the right margin, both anchored to the top so the
    ## module list grows downwards into empty margin rather than colliding with
    ## the arms key - the fix compare_maboss.py carries for the same reason.
    arms = [plt.Line2D([], [], color="0.3", lw=1.6, label="wild-type"),
            plt.Line2D([], [], color="0.3", lw=1.6, ls="--", label="knockout")]
    arms_legend = ax.legend(handles=arms, loc="upper left",
                            bbox_to_anchor=(1.0, 1.0), fontsize=8)
    ax.add_artist(arms_legend)
    mods_legend = ax.legend(loc="upper left", bbox_to_anchor=(1.0, 0.86),
                            fontsize=8, title="module")

    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, f"{stem}_{_safe(knockout)}_{_safe(tag)}_traj.png")
    ## A legend re-added with add_artist is not picked up by the tight bbox on
    ## its own - name both explicitly or long module labels get sliced off.
    fig.savefig(p, bbox_inches="tight", dpi=150,
                bbox_extra_artists=(mods_legend, arms_legend))
    plt.close(fig)
    return p


def plot_heatmap(traj, stem, cl, knockout, readouts, mods_idx, outdir,
                 tag, width=10.0, height=None):
    """Delta over time for the chosen readouts, as a heatmap.

    The view that scales: many more rows are legible here than as lines. Rows
    are ordered by size of final change, largest at the top.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if knockout not in traj["perturbed"]:
        return None
    k = traj["perturbed"].index(knockout)
    midx = {m: i for i, m in enumerate(traj["modules"])}
    sel = [m for m in readouts if m in midx]
    if not sel:
        return None

    cols = [midx[m] for m in sel]
    D = traj["traj"][k][:, cols] - traj["wt"][:, cols]      # time x module
    order = np.argsort(-np.abs(D[-1]))                      # by final change
    sel = [sel[i] for i in order]
    D = D[:, order].T                                       # module x time

    labelled = len(sel) <= 40                  # labels stop being legible past this
    per_row = 0.16 if labelled else 0.06       # only reserve height for real labels
    auto_height = max(3.0, min(per_row * len(sel) + 1.5, 18.0))
    fig, ax = plt.subplots(figsize=(width, height or auto_height))

    ## Colour scaled to the data and symmetric about zero. A fixed +-1 scale is
    ## the honest one but renders uniformly white whenever a knockout shifts
    ## things by a tenth - the normal case - so the range is stated on the
    ## colourbar instead of being assumed.
    lim = max(float(np.abs(D).max()), 0.02)
    ## shading="nearest" takes cell CENTRES, so the coordinate arrays match the
    ## data shape - with "flat" they would each have to be one longer.
    im = ax.pcolormesh(traj["times"], np.arange(len(sel)), D,
                       cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    ax.set_xlabel("time")
    ax.set_ylabel("module (largest final change at top)")
    ax.invert_yaxis()
    if labelled:
        ax.set_yticks(np.arange(len(sel)))
        ax.set_yticklabels([_label(m, mods_idx) for m in sel], fontsize=7)
    else:
        ax.set_yticks([])
    ax.set_title(f"{cl}: change in P(module = ON), {_label(knockout, mods_idx)} OFF "
                 f"- wild-type\n{tag}")
    fig.colorbar(im, ax=ax, label=f"Δ P(ON), scale ±{lim:.3f}", fraction=0.03, pad=0.02)

    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, f"{stem}_{_safe(knockout)}_{_safe(tag)}_heatmap.png")
    fig.savefig(p, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return p


def plot_tradeoff_scatter(combined, stem, cl, outdir, width=8.0, height=7.0, label_n=8):
    """Every knockout placed by what it does to BOTH processes at once.

    x = net change across proliferation modules, y = net change across apoptosis
    modules, both summed over significant pairs only. The quadrant that matters
    is upper-left (apoptosis up, proliferation down); the origin is where the
    great majority of knockouts sit, which is itself the result - most single
    knockouts do nothing to either process. Only the outliers are labelled,
    since labelling 200 points at the origin is a black smudge.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if combined is None or combined.empty:
        return None
    x = combined["pro_net_delta_final"].astype(float)
    y = combined["apo_net_delta_final"].astype(float)

    fig, ax = plt.subplots(figsize=(width, height))
    ## Shade the therapeutically interesting quadrant rather than only drawing
    ## axes through zero - the point of the figure is which quadrant a knockout
    ## lands in, and that reads faster as an area than as a sign convention.
    lim = max(float(np.abs(x).max()), float(np.abs(y).max()), 0.5) * 1.15
    ax.axhspan(0, lim, xmin=0, xmax=0.5, color="tab:green", alpha=0.07, zorder=0)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.scatter(x, y, s=28, c="tab:blue", alpha=0.65, edgecolors="none", zorder=3)

    ## Label the knockouts furthest from the origin in either process. Two
    ## knockouts landing on the same point is common - a pair that shares an
    ## out-component moves the same readouts by the same amount - so labels are
    ## stacked when their coordinates collide rather than overprinted, which
    ## otherwise renders as one unreadable word.
    dist = (x ** 2 + y ** 2) ** 0.5
    used = {}
    for mod in dist.sort_values(ascending=False).head(label_n).index:
        disp = combined.at[mod, "perturbed_display"]
        text = str(disp) if isinstance(disp, str) and disp and disp != mod else str(mod)
        key = (round(float(x[mod]), 3), round(float(y[mod]), 3))   # same point?
        n = used.get(key, 0)
        used[key] = n + 1
        ax.annotate(text, (x[mod], y[mod]), fontsize=8,
                    xytext=(5, 5 + 11 * n), textcoords="offset points")

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("net Δ P(ON) summed over PROLIFERATION modules")
    ax.set_ylabel("net Δ P(ON) summed over APOPTOSIS modules")
    ax.set_title(f"{cl} ({CELL_LINE_NAMES.get(cl, cl)}): every module knockout,\n"
                 f"apoptosis vs proliferation effect "
                 f"(shaded = apoptosis up, proliferation down)")
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, f"{stem}_tradeoff_scatter.png")
    fig.savefig(p, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return p


def _safe(s):
    """Filename-safe form of a module ID or tag (`18.37` -> `18_37`)."""
    return "".join(c if c.isalnum() else "_" for c in str(s)).strip("_")


## ---------------------------------------------------------------------------
## Driver
## ---------------------------------------------------------------------------

def process_model(cl, stem, args):
    """Run the whole annotate -> recalibrate -> aggregate -> plot chain on one model."""
    print(f"\n{'=' * 72}\n{cl}  ({CELL_LINE_NAMES.get(cl, cl)})  -  {stem}\n{'=' * 72}")

    ## GO annotation is used when annotate_go.R has produced it, and the fallback
    ## is announced rather than silent - which annotation a number came from is
    ## exactly the thing that must not be ambiguous.
    gene_flags, module_labels = (None, None) if args.no_go else \
        load_go_annotation(args.go_dir, cl)
    mods = annotate_modules(load_modules(stem, args.sweep_dir),
                            gene_flags, module_labels)
    pairs = load_pairs(stem, args.sweep_dir)
    traj = load_traj(stem, args.sweep_dir) if not args.no_figures else None

    src = mods["annotation_source"].iloc[0]
    print(f"modules: {len(mods)}   pairs: {len(pairs)}   "
          f"process annotation: {src}"
          + ("" if module_labels else "  (no GO module labels found)"))
    if src == "curated":
        print(f"  note: no GO annotation in {args.go_dir} - "
              f"run gabi-maboss/scripts/annotate_go.R to replace the "
              f"hand-curated gene sets")

    for proc in ("apoptosis", "proliferation"):
        n_core = int(mods[f"{proc}_core"].sum())
        n_cur = int(mods[f"{proc}_curated"].sum())
        n_ext = int(mods[f"{proc}_extended"].sum())
        n_cf = int(mods[f"{proc}_cf"].sum())
        line = (f"  {proc:14s} used={n_core:3d}  curated={n_cur:3d}  "
                f"extended={n_ext:3d}  legacy-cluster-function={n_cf:3d}")
        ## When both annotations exist, say how far apart they are. A large
        ## disagreement is the point of having re-annotated, not a warning sign.
        if f"{proc}_go" in mods:
            go, cur = mods[f"{proc}_go"], mods[f"{proc}_curated"]
            line += (f"   [GO adds {int((go & ~cur).sum())}, "
                     f"drops {int((cur & ~go).sum())}]")
        print(line)
    print(f"  dual-process (BIRC5-type): {int(mods['dual_process'].sum())}")

    ## --- significance recalibration -------------------------------------
    floors = empirical_floors(pairs)
    pairs, pooled_final, pooled_max = apply_significance(pairs, floors)
    sweep_flag = int(pairs["above_noise"].sum()) if "above_noise" in pairs else -1
    print(f"\nnoise recalibration")
    print(f"  empirical floor (|delta_final|, {NULL_PERCENTILE:g}th pct of unreachable):"
          f" median {float(floors['floor_final'].median()):.4f}, "
          f"range {float(floors['floor_final'].min()):.4f}-{float(floors['floor_final'].max()):.4f}")
    print(f"  pooled fallback floor: {pooled_final:.4f}")
    print(f"  pairs flagged: sweep's above_noise {sweep_flag} -> recalibrated "
          f"{int(pairs['sig_final'].sum())} of {len(pairs)} "
          f"({int(pairs['reachable'].sum())} reachable)")

    ## --- process-level aggregation --------------------------------------
    tables, readouts_tbl, pairs_tbl = {}, {}, {}
    for proc, col in [("apoptosis", "apoptosis_core"),
                      ("proliferation", "proliferation_core")]:
        agg = aggregate_by_process(pairs, mods, proc, col)
        if agg.empty:
            continue
        agg = attach_knockout_labels(agg, mods)
        agg.insert(0, "cell_line", cl)
        agg.insert(1, "model", stem)
        ## Ranked on the attractor effect first, then the transient one, so a
        ## knockout whose effect washes out entirely still sorts above the
        ## knockouts that never did anything at all.
        tables[proc] = agg.sort_values(["sum_abs_delta_final", "sum_abs_delta_max"],
                                       ascending=False)

        r = responsive_readouts(pairs, mods, proc, col)
        if not r.empty:
            r.insert(0, "cell_line", cl)
            readouts_tbl[proc] = r
        t = top_pairs(pairs, mods, proc, col)
        if not t.empty:
            t.insert(0, "cell_line", cl)
            pairs_tbl[proc] = t

    ## --- the trade-off view ---------------------------------------------
    ## One row per knockout carrying BOTH processes side by side. This is the
    ## table the therapeutic question is actually asked of: apoptosis up AND
    ## proliferation down is the quadrant that matters, and it cannot be read
    ## off either single-process ranking.
    combined = None
    if "apoptosis" in tables and "proliferation" in tables:
        a = tables["apoptosis"].set_index("perturbed")
        p = tables["proliferation"].set_index("perturbed")
        shared = ["cell_line", "model", "perturbed_display", "perturbed_function",
                  "perturbed_genes", "perturbed_P_ON_wt", "in_degree", "out_degree"]
        combined = a[shared].copy()
        for src, pre in [(a, "apo"), (p, "pro")]:
            for c in ["n_sig", "n_readouts", "net_delta_final", "sum_abs_delta_final",
                      "max_abs_delta_final", "mean_delta_per_module", "min_hop_sig"]:
                combined[f"{pre}_{c}"] = src[c]
        ## The score: how far into the "apoptosis up, proliferation down"
        ## quadrant a knockout sits. Per-module means, so networks with
        ## different numbers of apoptotic modules stay comparable.
        combined["tradeoff_score"] = (combined["apo_mean_delta_per_module"]
                                      - combined["pro_mean_delta_per_module"])
        combined = combined.sort_values("tradeoff_score", ascending=False)
        combined.index.name = "perturbed"

    ## --- write the tables -----------------------------------------------
    os.makedirs(args.query_dir, exist_ok=True)
    written = []
    for proc, tbl in tables.items():
        p = os.path.join(args.query_dir, f"{stem}_knockout_vs_{proc}.csv")
        tbl.to_csv(p, index=False)
        written.append(p)
    for proc, tbl in readouts_tbl.items():
        p = os.path.join(args.query_dir, f"{stem}_{proc}_module_responsiveness.csv")
        tbl.to_csv(p, index=False)
        written.append(p)
    for proc, tbl in pairs_tbl.items():
        p = os.path.join(args.query_dir, f"{stem}_{proc}_top_pairs.csv")
        tbl.to_csv(p, index=False)
        written.append(p)
    if combined is not None:
        p = os.path.join(args.query_dir, f"{stem}_tradeoff.csv")
        combined.to_csv(p)
        written.append(p)
    p = os.path.join(args.query_dir, f"{stem}_module_process_annotation.csv")
    mods.to_csv(p, index=False)
    written.append(p)
    p = os.path.join(args.query_dir, f"{stem}_noise_floors.csv")
    floors.to_csv(p)
    written.append(p)

    ## --- print the headline rankings ------------------------------------
    for proc, tbl in tables.items():
        print(f"\ntop knockouts by total effect on {proc.upper()} "
              f"({int(mods[proc + '_core'].sum())} readout modules):")
        show = tbl.head(args.top)[[
            "perturbed", "perturbed_display", "perturbed_function", "n_sig",
            "net_delta_final", "sum_abs_delta_final", "max_abs_delta_final",
            "min_hop_sig", "perturbed_genes"]].copy()
        show["perturbed_genes"] = show["perturbed_genes"].fillna("").str.slice(0, 40)
        show["perturbed_function"] = show["perturbed_function"].fillna("").str.slice(0, 26)
        print(show.to_string(index=False, float_format=lambda v: f"{v:7.3f}"))

        ## The transient ranking, printed separately because it selects a
        ## different set of knockouts - ones whose effect the network absorbs.
        tr = tbl[tbl["sum_abs_delta_max"] > 0].sort_values("sum_abs_delta_max",
                                                           ascending=False)
        print(f"\n  ...and by largest TRANSIENT effect on {proc} "
              f"(delta_max; effect may recover by the attractor):")
        show = tr.head(args.top)[["perturbed", "perturbed_display", "n_sig_max",
                                  "sum_abs_delta_max", "max_abs_delta_max",
                                  "sum_abs_delta_final", "t_at_max_swing"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:7.3f}"))

    ## Which specific process modules move, and what moves them - the readout
    ## side of the same database.
    for proc, tbl in readouts_tbl.items():
        ## Interpretable readouts first: a result on a node that IS the gene
        ## outranks the same result on a node where the gene is one of fifteen.
        ## Sorting by tier before magnitude is the practical form of "focus on
        ## the modules whose process membership we actually know".
        rank = {"single_gene": 0, "clean": 1, "mixed": 2, "none": 3}
        tbl = tbl.assign(_r=tbl[f"{proc}_tier"].map(rank).fillna(3)) \
                 .sort_values(["_r", "max_abs_delta_final"],
                              ascending=[True, False]).drop(columns="_r")
        readouts_tbl[proc] = tbl
        print(f"\nmost responsive {proc.upper()} modules, most interpretable first "
              f"(tier: single_gene > clean > mixed):")
        show = tbl.head(args.top)[[
            "readout", "Display", f"{proc}_tier", f"{proc}_purity",
            "Module_HGNCs", "P_ON_wt", "n_knockouts_moving_it",
            "max_abs_delta_final", "best_knockout_display", "best_delta_final",
            "best_hop"]].copy()
        show["Module_HGNCs"] = show["Module_HGNCs"].fillna("").str.slice(0, 30)
        print(show.to_string(index=False, float_format=lambda v: f"{v:7.3f}"))
        counts = tbl[f"{proc}_tier"].value_counts().to_dict()
        print(f"   tiers: {counts}")

    if combined is not None:
        print(f"\ntop knockouts by APOPTOSIS-UP / PROLIFERATION-DOWN trade-off:")
        show = combined.head(args.top)[[
            "perturbed_display", "perturbed_function", "apo_net_delta_final",
            "pro_net_delta_final", "tradeoff_score", "perturbed_genes"]].copy()
        show["perturbed_genes"] = show["perturbed_genes"].fillna("").str.slice(0, 40)
        show["perturbed_function"] = show["perturbed_function"].fillna("").str.slice(0, 26)
        print(show.to_string(float_format=lambda v: f"{v:7.3f}"))

    ## --- figures ---------------------------------------------------------
    ## Filed as <fig-dir>/<function>/<cell line>/, the same tree
    ## make_followup_figures.py writes into, so every figure in the project sits
    ## under the process it is about and then the network it came from. A flat
    ## folder of 45+ files sorted by model stem buried the biology.
    def fig_out(kind):
        return os.path.join(args.fig_dir, kind, cl)

    figs = []
    ## The overview first: it needs no trajectories, only the aggregate table,
    ## so it is drawn even under --no-figures' sibling case where the .npz is
    ## missing (an interrupted sweep, or a copy of the CSVs alone).
    if combined is not None and not args.no_figures:
        p = plot_tradeoff_scatter(combined, stem, cl, fig_out("overview"))
        if p:
            figs.append(p)
    if traj is not None:
        mods_idx = dict(zip(mods["Module_ID"], mods["Display"]))
        apo_readouts = list(mods.loc[mods["apoptosis_core"], "Module_ID"])
        pro_readouts = list(mods.loc[mods["proliferation_core"], "Module_ID"])

        ## Which knockouts get drawn: the strongest for each process, plus the
        ## strongest trade-off - chosen from the recalibrated tables rather than
        ## picked by hand, so the figure selection is reproducible.
        picks = []
        if "apoptosis" in tables:
            for m in tables["apoptosis"].head(args.n_figs)["perturbed"]:
                picks.append((m, "apoptosis readouts", apo_readouts))
        if "proliferation" in tables:
            for m in tables["proliferation"].head(args.n_figs)["perturbed"]:
                picks.append((m, "proliferation readouts", pro_readouts))
        if combined is not None and len(combined):
            top_trade = combined.index[0]
            picks.append((top_trade, "apoptosis and proliferation readouts",
                          apo_readouts + pro_readouts))

        ## Folder is taken from which readout set the figure is about, so the
        ## two single-process figures file under their process and the combined
        ## one under `overview`.
        folder_for = {"apoptosis readouts": "apoptosis",
                      "proliferation readouts": "proliferation",
                      "apoptosis and proliferation readouts": "overview"}

        seen = set()
        for knockout, tag, readouts in picks:
            key = (knockout, tag)
            if key in seen:
                continue
            seen.add(key)
            outdir = fig_out(folder_for.get(tag, "overview"))
            for fn in (plot_pair, plot_heatmap):
                p = fn(traj, stem, cl, knockout, readouts, mods_idx, outdir, tag)
                if p:
                    figs.append(p)

    for p in written + figs:
        print(f"-> {p}")
    return tables, combined, mods, pairs


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Process-level (apoptosis / proliferation) report over "
                    "sweep_maboss.py's knockout database.")
    parser.add_argument("--sweep-dir", default=SWEEP_DIR,
                        help=f"folder holding the sweep tables (default: {SWEEP_DIR})")
    parser.add_argument("--query-dir", default=QUERY_DIR,
                        help="where the derived tables are written")
    parser.add_argument("--fig-dir", default=FIG_DIR,
                        help="where the figures are written")
    parser.add_argument("--go-dir", default=GO_DIR,
                        help=f"annotate_go.R's output folder; process membership "
                             f"comes from GO when this exists (default: {GO_DIR})")
    parser.add_argument("--no-go", action="store_true",
                        help="ignore the GO annotation even if present and use the "
                             "hand-curated gene sets, for a like-for-like comparison "
                             "against the first version of this report")
    parser.add_argument("--top", type=int, default=12,
                        help="how many rows of each ranking to print (default 12)")
    parser.add_argument("--n-figs", type=int, default=2,
                        help="how many top knockouts per process get figures (default 2)")
    parser.add_argument("--no-figures", action="store_true",
                        help="skip the .npz load and the plots (tables only)")
    args = parser.parse_args(argv)

    all_trade = []
    for cl, stem in MODELS:
        if not os.path.isfile(os.path.join(args.sweep_dir, f"{stem}_sweep_OFF.csv")):
            print(f"{stem}: SKIPPED - no sweep table in {args.sweep_dir}")
            continue
        _, combined, _, _ = process_model(cl, stem, args)
        if combined is not None:
            all_trade.append(combined.reset_index())

    ## One concatenated trade-off table across the three networks. Safe to
    ## concatenate ONLY because every row carries its own cell_line and the gene
    ## list - module IDs come from each line's own STEM clustering and are not
    ## comparable between them, so the ID column must never be joined across.
    if all_trade:
        p = os.path.join(args.query_dir, "all_models_tradeoff.csv")
        pd.concat(all_trade, ignore_index=True).to_csv(p, index=False)
        print(f"\n-> {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
