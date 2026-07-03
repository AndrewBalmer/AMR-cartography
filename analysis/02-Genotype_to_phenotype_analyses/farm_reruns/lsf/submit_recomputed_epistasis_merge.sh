#!/usr/bin/env bash
set -euo pipefail

# Merge the completed recomputed epistasis observed/permutation chunks.
# This uses the recomputed 4,052-interaction Galwey meff and the original
# lowest-minimum permutation threshold policy, then applies the lower-bound
# effect filter to define supported interactions.

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/epistasis/chunks}"
MERGED_DIR="${MERGED_DIR:-$FARM_OUT/epistasis/merged}"
THRESHOLD_DIR="${THRESHOLD_DIR:-$FARM_OUT/thresholds}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

EPI_MEFF=$("$PYTHON" - "$THRESHOLD_DIR/recomputed_meff.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["epistasis_meff"])
PY
)

mkdir -p "$MERGED_DIR" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/epistasis_merge.%J.out" \
  -eo "$FARM_OUT/logs/epistasis_merge.%J.err" \
  "$PYTHON '$SCRIPT_DIR/merge_epistasis_chunks.py' \
    --chunk-dir '$CHUNK_DIR' \
    --out-dir '$MERGED_DIR' \
    --galwey-meff '$EPI_MEFF' \
    --threshold-mode lowest-min-p \
    --expected-observed-interactions 4052 \
    --expected-permutation-rows 405200 \
    --expected-permutations 100 \
    --fail-on-non-ok"
