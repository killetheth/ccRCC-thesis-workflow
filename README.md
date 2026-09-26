# Systems Medicine of Renal Cancer Drug Resistance — thesis workflow

Code for the MSc Bioinformatics and Computational Genomics dissertation
*Systems Medicine of Renal Cancer Drug Resistance: Towards New Diagnostics and
Therapy* (Rachael Adamson, Queen's University Belfast, 2026; supervisors
Dr Ian Overton and Dr Laurence Calzone).

The workflow turns time-course expression data for four cell lines (A498,
CAKI-1, HUVEC, UMRC2), each in hypoxia and normoxia (AH, AN, CH, CN, HH, HN,
UH, UN), into module-level causal networks with Gabi, converts the three
hypoxic ccRCC networks (AH, CH, UH) into Boolean models, and simulates them in
MaBoSS, including a knockout of every module in turn.

This repository holds only the scripts needed to reproduce the thesis: the
pipeline in the thesis's Table 1 and the analysis scripts behind the reported
numbers. Data are not included (see [Inputs](#inputs)).

## Layout

The folder layout matches the one the scripts were run in, because each script
works from paths relative to its own folder (`../data/`, `../output/`,
`../../data/` from a branch subfolder).

```
stem-gabi/scripts/           Stage 1: temporal profiles -> modules
gabi-maboss/scripts/         Stages 2-3: networks -> Boolean models -> simulation
  reactome/ humannet/ huri/  one scaffold per folder (steps 1-4)
  query/                     analysis scripts behind the reported numbers
```

## Workflow

Run the steps in this order. "Section" is the thesis Methods section that
describes each step.

| Stage | Section | Script | What it does |
|---|---|---|---|
| 1 | 2.1 | `stem-gabi/scripts/1-separate_genes.R` | Splits STEM gene tables into one file per temporal profile |
| 1 | 2.1 | `stem-gabi/scripts/2-genes_in_profile.R` | Counts significant genes per profile; builds each dataset's gene-profile table |
| 1 | 2.1 | `stem-gabi/scripts/3-cyto_module.R` (+ `cyto_module_common.R`) | Contracts genes into modules; writes gene-to-module tables and module mean expression |
| 2 | 2.2 | `gabi-maboss/scripts/<scaffold>/1-undir_scaff_<scaffold>.R` (+ `scaffold_common.R`) | Builds the undirected module scaffold from the prior-knowledge network |
| 2 | 2.3 | `gabi-maboss/scripts/<scaffold>/2-gabi_<scaffold>.R` (driver: `run_all_gabi_<scaffold>.sh`) | Orients the scaffold with Gabi (ClNetScOn) from module mean expression |
| 2 | 2.4 | `gabi-maboss/scripts/<scaffold>/3-sign_export_<scaffold>.R` | Adds edge confidence, unresolved-edge flags and activation/inhibition signs; exports GraphML |
| 2 | 2.5 | `gabi-maboss/scripts/<scaffold>/4-export_cyto_<scaffold>.R` | Exports each network back out of Cytoscape |
| 2 | 2.5 | `gabi-maboss/scripts/5-sort_modules.R` | Degree, Source/Relay/Sink roles and potential hubs for every module (Table 7, Table B2) |
| — | 2.6 | *manual* | In Cytoscape, File → Export → Table to File on the chosen network's **edge table** (`.csv`) |
| 3 | 2.6 | `gabi-maboss/scripts/6-cyto_to_neko.py` | Edge table → signed SIF with MaBoSS-safe names, plus a name map back to module IDs |
| 3 | 2.6 | `gabi-maboss/scripts/7-sif_to_bnet.py` | SIF → Boolean model (`.bnet`) with NeKo |
| 3 | 2.7 | `gabi-maboss/scripts/8-run_maboss.py` (+ `maboss_common.py`) | Simulates a model; attractor probabilities and fixed points |
| 3 | 2.7 | `gabi-maboss/scripts/compare_maboss.py` (batch: `run_all_maboss.sh`) | Paired WT and perturbed arms; per-module change in p(ON) |
| 3 | 2.7 | `gabi-maboss/scripts/trace_cascade.py` | Path length and sign agreement from the perturbed module to each module that moved |
| 3 | 2.8 | `gabi-maboss/scripts/sweep_maboss.py` | Knocks out every module in turn; writes the pair table and trajectories |

`<scaffold>` is `reactome`, `humannet` or `huri`. Only the Reactome networks go
on to Stage 3.

### Analysis scripts (`gabi-maboss/scripts/query/`)

| Script | Produces |
|---|---|
| `scaffold_comparison.py` | Orientation rates per scaffold and dataset (Table 6, Table B1) |
| `scaffold_overlap.py` | Edge overlap between Reactome and HumanNet, and the paired orientation comparison (Section 3.3) |
| `family_agreement.py` | Backbone vs bidirectional model agreement (Section 3.7, Table B3) |
| `sweep_process_report.py` | The sweep's per-knockout empirical noise floors (Section 2.8) and the process-level proliferation and apoptosis tables behind Section 3.8 |

## Inputs

Not included in this repository. Place them where the scripts expect them:

- **Expression data and STEM output** — pre-processed time-course data for the
  eight datasets (Dewey, 2019), clustered with STEM v1.3.13; the Cytoscape
  module exports that `3-cyto_module.R` reads.
- **Reactome FI network** — `gabi-maboss/data/FIsInGene_04142025_with_annotations.txt`
  (Reactome Functional Interactions, April 2025 release).
- **HumanNet v3 functional network** — `gabi-maboss/data/HumanNet-FN.tsv.gz`.
- **HuRI** — `gabi-maboss/data/HI-union.tsv` and `gabi-maboss/data/Lit-BM.tsv`.

Each script's header states its exact inputs and outputs.

## Software

Versions used for the thesis (Table 3):

- **R 4.6.0** — igraph 2.3.2, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0,
  stringr 1.6.0, doParallel 1.0.17, FNN 1.1.4.1, RCy3 2.30.1 (Cytoscape steps
  only), org.Hs.eg.db 3.22.0 and AnnotationDbi 1.72.0 (HuRI branch only)
- **ClNetScOn (Gabi) 1.0.1** — from the Overton group; not distributed here
- **Python 3.12** — pandas 2.3.3, numpy 2.4.0, maboss (pyMaBoSS) 0.8.16,
  nekomata 1.1.3 (imported as `neko`)
- **MaBoSS 2.6.6** engine on `PATH` for `8-run_maboss.py`, `compare_maboss.py`
  and `sweep_maboss.py`
- **Cytoscape 3.10.4**, **STEM 1.3.13**

`6-cyto_to_neko.py`, `trace_cascade.py` and the `query/` scripts need only
the Python standard library, numpy and pandas.

## Running notes

- Most R scripts pin their working directory with
  `setwd("C:/ccRCC-causal-networks/...")`. Change that line to your own copy's
  location, or run each script from its own folder with the line removed.
  The three `2-gabi_<scaffold>.R` scripts work out their own folder and run
  anywhere; they were run on a 12-core Linux workstation, one run at a time,
  because each uses all 12 cores.
- Steps 3–4 push to and pull from a running Cytoscape session through RCy3.
- Module IDs look numeric (`18.30`), so read the output tables as text
  (`dtype=str` in pandas, `colClasses = "character"` in R) or IDs such as
  `18.30` become `18.3`.

## Licence

No licence is granted. The code is shared for assessment and reproducibility
of the thesis; please contact the author before reusing it.
