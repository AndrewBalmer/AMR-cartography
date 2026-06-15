#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_FILE="${INTERACTION_FILE:-$FARM_OUT/interactions/corrected_epistasis_interactions.csv}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/epistasis/chunks}"
CHUNK_SIZE="${CHUNK_SIZE:-25}"
N_PERMUTATIONS="${N_PERMUTATIONS:-100}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$CHUNK_DIR" "$FARM_OUT/logs"
TOTAL=$("$PYTHON" - "$INTERACTION_FILE" <<'PY'
import pandas as pd, sys
print(len(pd.read_csv(sys.argv[1])))
PY
)
JOBS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "Submitting $N_PERMUTATIONS recomputed epistasis permutation arrays, $JOBS chunks each"

for PERM in $(seq 1 "$N_PERMUTATIONS"); do
  bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
    -J "recomp_epi_perm_${PERM}[1-$JOBS]" \
    -oo "$FARM_OUT/logs/epistasis_perm_${PERM}.%I.%J.out" \
    -eo "$FARM_OUT/logs/epistasis_perm_${PERM}.%I.%J.err" \
    "$PYTHON '$SCRIPT_DIR/run_epistasis_permutation_chunk.py' \
      --data-dir '$DATA_DIR' \
      --interaction-file '$INTERACTION_FILE' \
      --out-dir '$CHUNK_DIR' \
      --array-index \${LSB_JOBINDEX} \
      --chunk-size '$CHUNK_SIZE' \
      --permutation-index '$PERM'"
done
