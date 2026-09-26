#!/usr/bin/env bash
## Run the HumanNet gabi orient over every threshold x cell line on this single
## 12-core box. 2-gabi_humannet.R claims all 12 cores via
## gabi.scaffold.master(cores = 12), so runs are strictly sequential (never in
## parallel) to avoid CPU contention.
##
## Each run writes its own verbose log to
## logs/gabi/humannet/<CL>/<CL>_gabi_d<THR>.log; this driver adds a high-level
## START/DONE timestamp per run to run_all.log so you can watch batch progress
## without tailing every per-run log.
##
## Usage: ./run_all_gabi_humannet.sh   (run from gabi-maboss/scripts/humannet/,
##        e.g. inside tmux)
##
## Resume after a crash without redoing finished runs by naming the cell line
## and/or threshold to restart at (check run_all.log + the per-run log for the
## last START with no matching DONE):
##   START_CL=HN ./run_all_gabi_humannet.sh              # resume at cell line HN
##   START_THR=0.2 START_CL=HN ./run_all_gabi_humannet.sh # resume at 0.2 / HN
##
## The thresholds to run come from THRESHOLDS (space separated), so an
## alternative cutoff needs no edit here. The matching scaffold must exist
## first: build it with THR=<value> Rscript 1-undir_scaff_humannet.R.
##   THRESHOLDS="0" ./run_all_gabi_humannet.sh            # the d = 0 comparison
##   THRESHOLDS="0 0.2" ./run_all_gabi_humannet.sh        # both, in that order
## Note START_THR selects a resume point within THRESHOLDS; naming a value that
## is not in the list is an error rather than a silent no-op.

set -euo pipefail

read -r -a thresholds <<< "${THRESHOLDS:-0.2}"
cell_lines=(AH AN CH CN HH HN UH UN)

start_thr="${START_THR:-${thresholds[0]}}"
start_cl="${START_CL:-${cell_lines[0]}}"
started=false

## Guard against the silent no-op: if START_THR or START_CL names something not
## in the list, the loop below would run nothing and still report success.
printf '%s\n' "${thresholds[@]}" | grep -qx -- "$start_thr" \
  || { echo "START_THR='$start_thr' is not in THRESHOLDS='${thresholds[*]}'" >&2; exit 1; }
printf '%s\n' "${cell_lines[@]}" | grep -qx -- "$start_cl" \
  || { echo "START_CL='$start_cl' is not a known cell line" >&2; exit 1; }

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pipeline_root="$(cd "$script_dir/../.." && pwd)"
master_log="$pipeline_root/logs/gabi/humannet/run_all.log"
mkdir -p "$(dirname "$master_log")"

for thr in "${thresholds[@]}"; do
  for cl in "${cell_lines[@]}"; do
    if [[ "$started" != true ]]; then
      if [[ "$thr" == "$start_thr" && "$cl" == "$start_cl" ]]; then
        started=true
      else
        continue
      fi
    fi
    echo "[$(date '+%F %T')] START ${cl} d${thr}" | tee -a "$master_log"
    (cd "$script_dir" && Rscript 2-gabi_humannet.R "$cl" "$thr")
    echo "[$(date '+%F %T')] DONE  ${cl} d${thr}" | tee -a "$master_log"
  done
done

echo "[$(date '+%F %T')] All runs complete." | tee -a "$master_log"
