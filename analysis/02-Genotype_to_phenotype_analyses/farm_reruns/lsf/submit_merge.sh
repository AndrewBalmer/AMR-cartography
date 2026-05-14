#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the checked-out repo path}"
PYTHON="${PYTHON:-python}"
FARM_OUT="${FARM_OUT:-$PROJECT_ROOT/farm_outputs/corrected_epistasis}"
CHUNK_DIR="${CHUNK_DIR:-$FARM_OUT/chunks}"
MERGED_DIR="${MERGED_DIR:-$FARM_OUT/merged}"
QUEUE="${QUEUE:-normal}"
MEM_MB="${MEM_MB:-8000}"
SCRIPT_DIR="$PROJECT_ROOT/analysis/02-Genotype_to_phenotype_analyses/farm_reruns"

mkdir -p "$MERGED_DIR" "$FARM_OUT/logs"
bsub -q "$QUEUE" -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
  -oo "$FARM_OUT/logs/merge.%J.out" \
  -eo "$FARM_OUT/logs/merge.%J.err" \
  "$PYTHON '$SCRIPT_DIR/merge_epistasis_chunks.py' \
    --chunk-dir '$CHUNK_DIR' \
    --out-dir '$MERGED_DIR'"
