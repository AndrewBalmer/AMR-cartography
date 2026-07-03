#!/usr/bin/env bash
# Head-node-safe verification: low-rank epistatic mvLMM engine vs exact limix dense
# engine. Fits LMMs, so it is submitted through LSF and must never be run on a
# farm head node directly (see README "Verification runs" exception).
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_FILE="${INTERACTION_FILE:-$FARM_OUT/interactions/corrected_epistasis_interactions.csv}"
OBSERVED_FILE="${OBSERVED_FILE:-$FARM_OUT/epistasis/merged/corrected_epistasis_p_values.csv}"
OUT_DIR="${OUT_DIR:-$FARM_OUT/validation/engine_equivalence}"
N_LOWEST="${N_LOWEST:-10}"
N_RANDOM="${N_RANDOM:-10}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$OUT_DIR" "$FARM_OUT/logs"

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/verify_engine.%J.out" \
  -eo "$FARM_OUT/logs/verify_engine.%J.err" \
  "$PYTHON '$SCRIPT_DIR/verify_epistasis_engine_equivalence.py' \
    --data-dir '$DATA_DIR' \
    --interaction-file '$INTERACTION_FILE' \
    --observed-file '$OBSERVED_FILE' \
    --out-dir '$OUT_DIR' \
    --n-lowest $N_LOWEST \
    --n-random $N_RANDOM"
