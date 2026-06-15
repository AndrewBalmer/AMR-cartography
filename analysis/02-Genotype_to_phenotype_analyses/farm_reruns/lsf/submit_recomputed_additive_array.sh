#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/additive/chunks}"
CHUNK_SIZE="${CHUNK_SIZE:-5}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
MAX_CONCURRENT="${MAX_CONCURRENT:-100}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$CHUNK_DIR" "$FARM_OUT/logs"
TOTAL=$("$PYTHON" - "$DATA_DIR" <<'PY'
from pathlib import Path
import pandas as pd
import sys
print(pd.read_csv(Path(sys.argv[1]) / "S.pneumo_map_dummy_gen_test_markers.csv", nrows=0).shape[1])
PY
)
JOBS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "Submitting $JOBS additive chunks for $TOTAL markers"

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -J "recomp_add[1-$JOBS]%$MAX_CONCURRENT" \
  -oo "$FARM_OUT/logs/additive.%I.%J.out" \
  -eo "$FARM_OUT/logs/additive.%I.%J.err" \
  "$PYTHON '$SCRIPT_DIR/run_additive_chunk.py' \
    --data-dir '$DATA_DIR' \
    --out-dir '$CHUNK_DIR' \
    --array-index \${LSB_JOBINDEX} \
    --chunk-size '$CHUNK_SIZE'"
