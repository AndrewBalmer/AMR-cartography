#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
INTERACTION_DIR="${INTERACTION_DIR:-$FARM_OUT/interactions}"
THRESHOLD_DIR="${THRESHOLD_DIR:-$FARM_OUT/thresholds}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$INTERACTION_DIR" "$THRESHOLD_DIR" "$FARM_OUT/logs"

BSUB_DEPENDENCY_ARGS=()
if [[ -n "${DEPENDENCY:-}" ]]; then
  BSUB_DEPENDENCY_ARGS=(-w "$DEPENDENCY")
fi

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/recomputed_meff.%J.out" \
  -eo "$FARM_OUT/logs/recomputed_meff.%J.err" \
  "${BSUB_DEPENDENCY_ARGS[@]}" <<EOF
set -euo pipefail
"$PYTHON" "$SCRIPT_DIR/generate_corrected_epistasis_interactions.py" \\
  --data-dir "$DATA_DIR" \\
  --out-dir "$INTERACTION_DIR"
"$RSCRIPT" "$SCRIPT_DIR/compute_recomputed_meff.R" \\
  "$DATA_DIR" \\
  "$INTERACTION_DIR/corrected_epistasis_interactions.csv" \\
  "$THRESHOLD_DIR/recomputed_meff.json"
EOF
