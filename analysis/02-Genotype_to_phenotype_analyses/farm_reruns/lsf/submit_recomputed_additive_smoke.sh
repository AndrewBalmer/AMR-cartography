#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/additive/chunks}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$CHUNK_DIR" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/additive_smoke.%J.out" \
  -eo "$FARM_OUT/logs/additive_smoke.%J.err" \
  "$PYTHON '$SCRIPT_DIR/run_additive_chunk.py' \
    --data-dir '$DATA_DIR' \
    --out-dir '$CHUNK_DIR' \
    --start 0 \
    --end 1 \
    --chunk-label smoke \
    --force"
