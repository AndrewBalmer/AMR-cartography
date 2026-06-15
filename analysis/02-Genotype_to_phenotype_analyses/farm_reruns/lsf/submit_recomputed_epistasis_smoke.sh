#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_DIR="${INTERACTION_DIR:-$FARM_OUT/interactions}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/epistasis/chunks}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$INTERACTION_DIR" "$CHUNK_DIR" "$FARM_OUT/logs"
"$PYTHON" "$SCRIPT_DIR/generate_corrected_epistasis_interactions.py" \
  --data-dir "$DATA_DIR" \
  --out-dir "$INTERACTION_DIR"

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/epistasis_smoke.%J.out" \
  -eo "$FARM_OUT/logs/epistasis_smoke.%J.err" \
  "$PYTHON '$SCRIPT_DIR/run_epistasis_chunk.py' \
    --data-dir '$DATA_DIR' \
    --interaction-file '$INTERACTION_DIR/corrected_epistasis_interactions.csv' \
    --out-dir '$CHUNK_DIR' \
    --start 0 \
    --end 5 \
    --chunk-label smoke \
    --force"
