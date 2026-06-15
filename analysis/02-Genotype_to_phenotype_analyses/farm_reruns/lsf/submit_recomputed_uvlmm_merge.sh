#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/uvlmm/chunks}"
MERGED_DIR="${MERGED_DIR:-$FARM_OUT/uvlmm/merged}"
THRESHOLD_DIR="${THRESHOLD_DIR:-$FARM_OUT/thresholds}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

ADD_MEFF=$("$PYTHON" - "$THRESHOLD_DIR/recomputed_meff.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["additive_meff"])
PY
)

mkdir -p "$MERGED_DIR" "$THRESHOLD_DIR" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/uvlmm_merge.%J.out" \
  -eo "$FARM_OUT/logs/uvlmm_merge.%J.err" \
  "$PYTHON '$SCRIPT_DIR/merge_exact_unilmm_chunks.py' \
    --chunk-dir '$CHUNK_DIR' \
    --out-dir '$MERGED_DIR' \
    --threshold-dir '$THRESHOLD_DIR' \
    --galwey-meff '$ADD_MEFF'"
