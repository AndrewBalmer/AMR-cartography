#!/usr/bin/env bash
set -euo pipefail

# Required:
#   export PROJECT_ROOT=/path/to/AMR-cartography
# Optional:
#   export PYTHON=/path/to/python
#   export FARM_OUT=/path/to/output

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/original_logic_rebuild}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_DIR="${INTERACTION_DIR:-$FARM_OUT/interactions}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/chunks}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$INTERACTION_DIR" "$CHUNK_DIR" "$FARM_OUT/logs"

"$PYTHON" "$SCRIPT_DIR/generate_corrected_epistasis_interactions.py" \
  --data-dir "$DATA_DIR" \
  --out-dir "$INTERACTION_DIR"

bsub -q normal -M 8000 -R "select[mem>8000] rusage[mem=8000]" \
  -oo "$FARM_OUT/logs/smoke.%J.out" \
  -eo "$FARM_OUT/logs/smoke.%J.err" \
  "$PYTHON '$SCRIPT_DIR/run_epistasis_chunk.py' \
    --data-dir '$DATA_DIR' \
    --interaction-file '$INTERACTION_DIR/corrected_epistasis_interactions.csv' \
    --out-dir '$CHUNK_DIR' \
    --start 0 \
    --end 5 \
    --chunk-label smoke \
    --force"
