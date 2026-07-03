#!/usr/bin/env bash
set -euo pipefail

# Build the manuscript-review outputs after additive, uvLMM, epistasis, and
# threshold merges are complete. This writes into the recomputed FARM_OUT review
# directory only; it intentionally does not copy files into manuscript/.

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112}"
OLD_DIR="${OLD_DIR:-$PROJECT_ROOT/AMRC-repo-files/pythonProject1}"
SUPPORT_DIR="${SUPPORT_DIR:-$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis}"
ADDITIVE_DIR="${ADDITIVE_DIR:-$FARM_OUT/additive/merged}"
UV_DIR="${UV_DIR:-$FARM_OUT/uvlmm/merged}"
EPISTASIS_DIR="${EPISTASIS_DIR:-$FARM_OUT/epistasis/merged}"
THRESHOLD_DIR="${THRESHOLD_DIR:-$FARM_OUT/thresholds}"
OUTPUT_DIR="${OUTPUT_DIR:-$FARM_OUT/manuscript_outputs}"
ORIGINAL_LOGIC_PUBLIC="${ORIGINAL_LOGIC_PUBLIC:-$PROJECT_ROOT/farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

ADD_MEFF=$("$PYTHON" - "$THRESHOLD_DIR/recomputed_meff.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["additive_meff"])
PY
)
ADD_THRESHOLD=$("$PYTHON" - "$THRESHOLD_DIR/additive_threshold_summary.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["adjusted_threshold"])
PY
)

mkdir -p "$OUTPUT_DIR" "$FARM_OUT/validation" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/recomputed_review_outputs.%J.out" \
  -eo "$FARM_OUT/logs/recomputed_review_outputs.%J.err" \
  "$PYTHON '$SCRIPT_DIR/build_recomputed_summary_report.py' \
    --farm-out '$FARM_OUT' \
    --original-logic-public '$ORIGINAL_LOGIC_PUBLIC' \
    --out-dir '$OUTPUT_DIR' && \
   $PYTHON '$SCRIPT_DIR/build_corrected_evidence_table.py' \
    --legacy-dir '$OLD_DIR' \
    --new-dir '$DATA_DIR' \
    --support-dir '$SUPPORT_DIR' \
    --analysis-out '$EPISTASIS_DIR' \
    --additive-dir '$ADDITIVE_DIR' \
    --uv-dir '$UV_DIR' \
    --epistasis-dir '$EPISTASIS_DIR' \
    --threshold-file '$THRESHOLD_DIR/recomputed_thresholds.json' \
    --output-dir '$OUTPUT_DIR' \
    --additive-threshold '$ADD_THRESHOLD' \
    --additive-galwey-meff '$ADD_MEFF' \
    --analysis-label 'Recomputed 170-marker analysis' && \
   $PYTHON '$SCRIPT_DIR/build_recomputed_summary_report.py' \
    --farm-out '$FARM_OUT' \
    --original-logic-public '$ORIGINAL_LOGIC_PUBLIC' \
    --out-dir '$OUTPUT_DIR' && \
   $PYTHON '$SCRIPT_DIR/validate_recomputed_outputs.py' \
    --farm-out '$FARM_OUT'"
