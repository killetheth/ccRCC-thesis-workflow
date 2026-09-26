#!/usr/bin/env bash
## Run the Reactome gabi orient over every threshold x cell line on this single
## 12-core box. 2-gabi_reactome.R claims all 12 cores via
## gabi.scaffold.master(cores = 12), so runs are strictly sequential (never in
## parallel) to avoid CPU contention.
##
## Reactome/HGNC counterpart of ../humannet/run_all_gabi_humannet.sh: it drives
## the d0 (no-threshold) orient and writes to the reactome/ log subtree, so it
## never clobbers the HumanNet batch.
##
## Each run writes its own verbose log to
## logs/gabi/reactome/<CL>/<CL>_gabi_d<THR>.log; this driver adds a high-level
## START/DONE timestamp per run to run_all.log so you can watch batch progress
## without tailing every per-run log.
##
## Usage: ./run_all_gabi_reactome.sh   (run from gabi-maboss/scripts/reactome/,
##        e.g. inside tmux)
##
## Resume after a crash without redoing finished runs by naming the cell line
## and/or threshold to restart at (check run_all.log + the per-run log for the
## last START with no matching DONE):
##   START_CL=HN ./run_all_gabi_reactome.sh             # resume at cell line HN
##   START_THR=0 START_CL=HN ./run_all_gabi_reactome.sh  # resume at 0 / HN

set -euo pipefail

thresholds=(0)
cell_lines=(AH AN CH CN HH HN UH UN)

start_thr="${START_THR:-${thresholds[0]}}"
start_cl="${START_CL:-${cell_lines[0]}}"
started=false

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pipeline_root="$(cd "$script_dir/../.." && pwd)"
master_log="$pipeline_root/logs/gabi/reactome/run_all.log"
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
    (cd "$script_dir" && Rscript 2-gabi_reactome.R "$cl" "$thr")
    echo "[$(date '+%F %T')] DONE  ${cl} d${thr}" | tee -a "$master_log"
  done
done

echo "[$(date '+%F %T')] All runs complete." | tee -a "$master_log"
