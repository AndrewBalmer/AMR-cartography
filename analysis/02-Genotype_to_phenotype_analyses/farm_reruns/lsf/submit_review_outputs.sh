#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/original_logic_rebuild}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
OLD_DIR="${OLD_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1}"
SUPPORT_DIR="${SUPPORT_DIR:-$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis}"
MERGED_DIR="${MERGED_DIR:-$FARM_OUT/merged}"
OUTPUT_DIR="${OUTPUT_DIR:-$FARM_OUT/manuscript_outputs}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$MERGED_DIR" "$OUTPUT_DIR" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/review_outputs.%J.out" \
  -eo "$FARM_OUT/logs/review_outputs.%J.err" \
  "$PYTHON '$SCRIPT_DIR/merge_exact_unilmm.py' \
    --legacy-dir '$OLD_DIR' \
    --new-dir '$DATA_DIR' \
    --out-dir '$MERGED_DIR' && \
   $PYTHON '$SCRIPT_DIR/build_corrected_evidence_table.py' \
    --legacy-dir '$OLD_DIR' \
    --new-dir '$DATA_DIR' \
    --support-dir '$SUPPORT_DIR' \
    --analysis-out '$MERGED_DIR' \
    --output-dir '$OUTPUT_DIR'"
