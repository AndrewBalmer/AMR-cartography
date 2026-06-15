#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
RSCRIPT="${RSCRIPT:-Rscript}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-4000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$FARM_OUT/logs"

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/recomputed_install_poolr.%J.out" \
  -eo "$FARM_OUT/logs/recomputed_install_poolr.%J.err" \
  "$RSCRIPT '$SCRIPT_DIR/install_r_poolr.R'"
