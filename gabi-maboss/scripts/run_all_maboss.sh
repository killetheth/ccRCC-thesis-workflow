#!/usr/bin/env bash
## Run a batch of compare_maboss.py perturbations back to back, one log per run.
##
## Each entry in the `runs` table below is one simulation: a model stem, the
## perturbation flags for that run, and a short tag used to name its log. To
## test different modules, edit the table - nothing else in this script needs
## touching. The tag is free text; keep it short and filename-safe, since it
## becomes the log filename.
##
## Logging mirrors the gabi drivers: every run writes its full stdout/stderr to
## logs/maboss/<model_stem>/<tag>.log, and this driver adds a high-level
## START/DONE/FAIL line per run to logs/maboss/run_all.log so batch progress can
## be watched without tailing every per-run log.
##
## Unlike the gabi drivers this does NOT abort the batch on a failed run - the
## simulations are independent, so a model that fails to build shouldn't cost
## you the rest of the queue. Failures are recorded as FAIL in run_all.log and
## summarised at the end, and the script exits non-zero if any run failed.
##
## Usage (from gabi-maboss/scripts/, inside the activated `maboss` conda env):
##   ./run_all_maboss.sh
##
## Resume after a crash without redoing finished runs by naming the tag to
## restart at (check run_all.log for the last START with no matching DONE):
##   START_TAG=ch_cebpb_on ./run_all_maboss.sh
##
## Other overrides:
##   COMMON_FLAGS="--free-inputs --csv" ./run_all_maboss.sh   # different defaults
##   DRY_RUN=1 ./run_all_maboss.sh                            # print, don't run

set -uo pipefail

## ---------------------------------------------------------------------------
## Runs to perform: "<model stem> | <perturbation flags> | <tag>"
##
## The model stem is the .bnet filename without its extension, resolved against
## MODEL_DIR below. The perturbation flags are whatever distinguishes this run
## from the others (--compare / --mutate); shared flags live in COMMON_FLAGS so
## they don't have to be repeated on every line.
## ---------------------------------------------------------------------------
runs=(
  # --- AH: the A498 JUN hub - does the UMRC2 JUN result generalise across cell lines?
  #     Also a direct activator of A498's ABCB1 module (6.48).
  "AH_reactome_backbone_subnet_drop_unsigned_edge   | --compare 18.48=ON            | ah_jun_on"
  "AH_reactome_bi_subnet_drop_unsigned_edge         | --compare 18.48=ON            | ah_bi_jun_on"
  # --- AH: NRF2 knockout
  "AH_reactome_backbone_subnet_drop_unsigned_edge   | --compare 22.3_12=OFF         | ah_nrf2_off"
  "AH_reactome_bi_subnet_drop_unsigned_edge         | --compare 22.3_12=OFF         | ah_bi_nrf2_off"
  # --- AH: the insulated-APAF1 control (expected to move little beyond its island)
  "AH_reactome_backbone_subnet_drop_unsigned_edge   | --compare 18.14=ON            | ah_apaf1_on"
  "AH_reactome_bi_subnet_drop_unsigned_edge         | --compare 18.14=ON            | ah_bi_apaf1_on"
  # --- AH: the A498 CREB1/MYC hub - counterpart of ch_creb1_off, highest out-degree
  #     in A498 (40).
  "AH_reactome_backbone_subnet_drop_unsigned_edge   | --compare 18.3_12=OFF         | ah_creb1_off"
  "AH_reactome_bi_subnet_drop_unsigned_edge         | --compare 18.3_12=OFF         | ah_bi_creb1_off"
  # --- CH: XPO1/PKM hub knockout - selinexor's target; 21 of its 29 out-edges are
  #     inhibitory, so removing it should de-repress widely.
  "CH_reactome_backbone_subnet_edge                 | --compare 21.36=OFF           | ch_xpo1_off"
  "CH_reactome_bi_subnet_edge                       | --compare 21.36=OFF           | ch_bi_xpo1_off"
  # --- CH: HIF1A knockout - targeted test of the HIF1A -| ABCB1 edge. Small
  #     out-component (8 modules), so expect a local result, not a cascade.
  "CH_reactome_backbone_subnet_edge                 | --compare 6.45=OFF            | ch_hif1a_off"
  "CH_reactome_bi_subnet_edge                       | --compare 6.45=OFF            | ch_bi_hif1a_off"
  # --- CH: CEBPB forced on - four pro-death outputs
  "CH_reactome_backbone_subnet_edge                 | --compare 26.8=ON             | ch_cebpb_on"
  "CH_reactome_bi_subnet_edge                       | --compare 26.8=ON             | ch_bi_cebpb_on"
  # --- CH: remove the sole activator of BCL2 and BCL-xL
  "CH_reactome_backbone_subnet_edge                 | --compare 18.36=OFF           | ch_creb1_off"
  "CH_reactome_bi_subnet_edge                       | --compare 18.36=OFF           | ch_bi_creb1_off"
  # --- CH: does CEBPB activation still kill with the CREB1 input already gone?
  "CH_reactome_backbone_subnet_edge                 | --compare 26.8=ON --mutate 18.36=OFF | ch_cebpb_on_creb1_bg"
  "CH_reactome_bi_subnet_edge                       | --compare 26.8=ON --mutate 18.36=OFF | ch_bi_cebpb_on_creb1_bg"
  # --- UH: JUN lever - also the hypoxia-axis run
  "UH_reactome_backbone_subnet_edge                 | --compare 18.37=ON            | uh_jun_on"
  "UH_reactome_bi_subnet_edge                       | --compare 18.37=ON            | uh_bi_jun_on"
  # --- UH: HSP90AA1 knockout, de-represses the APAF1 module
  "UH_reactome_backbone_subnet_edge                 | --compare 21.14=OFF           | uh_hsp90aa1_off"
  "UH_reactome_bi_subnet_edge                       | --compare 21.14=OFF           | uh_bi_hsp90aa1_off"
  )

## Flags applied to every run in the table above.
##
## MaBoSS samples independent stochastic trajectories, so --threads splits the
## sample_count over cores and scales near-linearly. Runs stay sequential (one
## run at a time, each using the whole box) rather than running several models
## in parallel on one core each - same total work, but you get each result as
## it finishes instead of all seven at the end.
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
COMMON_FLAGS="${COMMON_FLAGS:---free-inputs --noise-floor --csv --plot --top 25 --threads $THREADS}"

## ---------------------------------------------------------------------------

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pipeline_root="$(cd "$script_dir/.." && pwd)"
model_dir="${MODEL_DIR:-$pipeline_root/output/for_maboss}"
log_root="$pipeline_root/logs/maboss"
master_log="$log_root/run_all.log"
mkdir -p "$log_root"

## The engine binary is separate from the Python package - a missing MaBoSS on
## PATH fails every run identically, so say so once up front rather than seven
## times in seven logs.
if ! command -v MaBoSS >/dev/null 2>&1; then
  echo "ERROR: the MaBoSS engine binary is not on PATH." >&2
  echo "       Activate the conda env first:  conda activate maboss" >&2
  exit 1
fi

start_tag="${START_TAG:-}"
started=true
[[ -n "$start_tag" ]] && started=false

failed=()
completed=0

echo "[$(date '+%F %T')] BATCH START (${#runs[@]} runs queued)" | tee -a "$master_log"

for entry in "${runs[@]}"; do
  IFS='|' read -r stem flags tag <<< "$entry"
  # Trim the padding whitespace the table uses for readability.
  stem="$(echo "$stem" | xargs)"
  flags="$(echo "$flags" | xargs)"
  tag="$(echo "$tag" | xargs)"

  if [[ "$started" != true ]]; then
    if [[ "$tag" == "$start_tag" ]]; then
      started=true
    else
      echo "[$(date '+%F %T')] SKIP  ${tag} (resuming at ${start_tag})" | tee -a "$master_log"
      continue
    fi
  fi

  model="$model_dir/${stem}.bnet"
  run_log="$log_root/${stem}/${tag}.log"
  mkdir -p "$(dirname "$run_log")"

  if [[ ! -f "$model" ]]; then
    echo "[$(date '+%F %T')] FAIL  ${tag} - no such model: ${model}" | tee -a "$master_log"
    failed+=("$tag (missing model)")
    continue
  fi

  echo "[$(date '+%F %T')] START ${tag} - ${stem} ${flags}" | tee -a "$master_log"

  if [[ -n "${DRY_RUN:-}" ]]; then
    echo "         would run: python 8-run_maboss.py $model $flags $COMMON_FLAGS"
    continue
  fi

  # Record what produced the log, so a log read months later is self-describing.
  {
    echo "### ${tag}"
    echo "### $(date '+%F %T')"
    echo "### python compare_maboss.py ${model} ${flags} ${COMMON_FLAGS}"
    echo
  } > "$run_log"

  if (cd "$script_dir" && python compare_maboss.py "$model" $flags $COMMON_FLAGS) >> "$run_log" 2>&1; then
    echo "[$(date '+%F %T')] DONE  ${tag} -> ${run_log}" | tee -a "$master_log"
    completed=$((completed + 1))
  else
    echo "[$(date '+%F %T')] FAIL  ${tag} - see ${run_log}" | tee -a "$master_log"
    failed+=("$tag")
  fi
done

echo "[$(date '+%F %T')] BATCH END - ${completed} completed, ${#failed[@]} failed" | tee -a "$master_log"

if (( ${#failed[@]} > 0 )); then
  for f in "${failed[@]}"; do
    echo "  FAILED: $f" | tee -a "$master_log"
  done
  exit 1
fi
