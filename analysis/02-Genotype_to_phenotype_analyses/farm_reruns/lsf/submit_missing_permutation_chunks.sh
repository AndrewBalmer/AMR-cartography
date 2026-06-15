#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/original_logic_rebuild}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_FILE="${INTERACTION_FILE:-$FARM_OUT/interactions/corrected_epistasis_interactions.csv}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/chunks}"
CHUNK_SIZE="${CHUNK_SIZE:-25}"
N_PERMUTATIONS="${N_PERMUTATIONS:-100}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
ON_FIT_ERROR="${ON_FIT_ERROR:-pvalue-one}"
FIT_TIMEOUT_SECONDS="${FIT_TIMEOUT_SECONDS:-300}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"
MISSING_FILE="${MISSING_FILE:-$FARM_OUT/missing_permutation_chunks.tsv}"

mkdir -p "$CHUNK_DIR" "$FARM_OUT/logs"

"$PYTHON" - "$CHUNK_DIR" "$N_PERMUTATIONS" "$CHUNK_SIZE" "$INTERACTION_FILE" "$MISSING_FILE" <<'PY'
from pathlib import Path
import math
import sys

import pandas as pd

chunk_dir = Path(sys.argv[1])
n_permutations = int(sys.argv[2])
chunk_size = int(sys.argv[3])
interaction_file = Path(sys.argv[4])
missing_file = Path(sys.argv[5])

total = len(pd.read_csv(interaction_file))
n_chunks = math.ceil(total / chunk_size)
missing: list[tuple[int, int]] = []
for perm in range(1, n_permutations + 1):
    for array_index in range(1, n_chunks + 1):
        expected = chunk_dir / f"epistasis_perm_p_values.chunk_perm_{perm:04d}_array_{array_index:06d}.csv"
        if not expected.exists():
            missing.append((perm, array_index))

missing_file.write_text(
    "permutation_index\tarray_index\n"
    + "".join(f"{perm}\t{array_index}\n" for perm, array_index in missing)
)
print(f"missing_chunks={len(missing)} n_chunks={n_chunks} total_interactions={total}")
PY

tail -n +2 "$MISSING_FILE" | while IFS=$'\t' read -r PERM ARRAY_INDEX; do
  [ -n "$PERM" ] || continue
  bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
    -J "corr_epi_perm_missing_${PERM}_${ARRAY_INDEX}" \
    -oo "$FARM_OUT/logs/perm_missing_${PERM}_${ARRAY_INDEX}.%J.out" \
    -eo "$FARM_OUT/logs/perm_missing_${PERM}_${ARRAY_INDEX}.%J.err" \
    "$PYTHON '$SCRIPT_DIR/run_epistasis_permutation_chunk.py' \
      --data-dir '$DATA_DIR' \
      --interaction-file '$INTERACTION_FILE' \
      --out-dir '$CHUNK_DIR' \
      --array-index '$ARRAY_INDEX' \
      --chunk-size '$CHUNK_SIZE' \
      --permutation-index '$PERM' \
      --on-fit-error '$ON_FIT_ERROR' \
      --fit-timeout-seconds '$FIT_TIMEOUT_SECONDS'"
done
