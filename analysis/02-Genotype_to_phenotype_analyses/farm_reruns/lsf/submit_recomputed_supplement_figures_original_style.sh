#!/usr/bin/env bash
set -euo pipefail

# Regenerate S17/S18/S19 using the original R-style visual grammar and the
# recomputed 170-marker outputs. This is the publication-review figure path;
# the older Python plotting script is diagnostic only.

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
RSCRIPT="${RSCRIPT:-Rscript}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds}"
OUTPUT_DIR="${OUTPUT_DIR:-$FARM_OUT/manuscript_outputs/supplement_figures/original_style}"
FIGURES="${FIGURES:-S17,S18,S19}"
INSTALL_DEPS="${INSTALL_DEPS:-true}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$OUTPUT_DIR" "$FARM_OUT/logs"

if [[ "$INSTALL_DEPS" == "true" ]]; then
  PREP_CMD="$RSCRIPT '$SCRIPT_DIR/install_r_figure_packages.R' && "
else
  PREP_CMD=""
fi

bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/recomputed_supplement_figures_original_style.%J.out" \
  -eo "$FARM_OUT/logs/recomputed_supplement_figures_original_style.%J.err" \
  "${PREP_CMD}FIGURES='$FIGURES' $RSCRIPT '$SCRIPT_DIR/generate_recomputed_supplement_figures_original_style.R' '$PROJECT_ROOT' '$FARM_OUT' '$OUTPUT_DIR'"
