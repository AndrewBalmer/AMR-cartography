# Recomputed 170-marker mvLMM/epistasis workflow

This directory contains the reproducible workflow for the revised
*Streptococcus pneumoniae* PBP genotype-to-phenotype analysis. The scripts are
ordinary Python and R entry points and can be run on a workstation, lab server,
or computing cluster.

The workflow has two roles:

1. Original-logic validation: reproduce the original 157-marker manuscript
   output frame and evidence counts as a quality-control check.
2. Primary recomputed analysis: rerun the additive mvLMM, exact uvLMM,
   epistatic mvLMM, permutation thresholds, manuscript tables, and supporting
   figures on the corrected 170-marker panel.

Older R/Rmd scripts elsewhere in `analysis/` are retained as historical
manuscript-generation scripts. The corrected 170-marker analysis is implemented
here.

## Compute requirements

The additive and epistatic mvLMM steps are computationally intensive. On a
single workstation or server, run them in chunks and monitor CPU, memory, and
wall time. On a shared cluster, adapt the same chunk commands to the local job
scheduler.

## Inputs and environment

Large input bundles and analysis outputs are not committed. The expected data files
are listed in `required_data_manifest.txt`.

Set these variables before running the workflow:

```bash
export PROJECT_ROOT=/path/to/AMR-cartography
export PYTHON=$PROJECT_ROOT/.venv-mvlmm/bin/python
export DATA_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112
export OLD_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1
export SUPPORT_DIR="$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis"
export RESULTS_DIR=$PROJECT_ROOT/analysis_outputs/recomputed_170_thresholds
export N_PERMUTATIONS=100
```

Create a Python environment with the required Limix stack:

```bash
python3 -m venv .venv-mvlmm
source .venv-mvlmm/bin/activate
pip install -r analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/requirements.txt
```

If a suitable conda environment or system Python already exists, set `PYTHON`
to that interpreter instead.

## Original-logic validation

Run this quality-control workflow before using revised manuscript outputs. It
checks that the code reproduces the original manuscript counting frame on the
old snapshot, then verifies that the corrected 170-marker panel contains the old
157 markers plus the 13 added markers.

```bash
export RESULTS_DIR=$PROJECT_ROOT/analysis_outputs/original_logic_rebuild

$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/compare_additive_mvlmm.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --out-dir "$RESULTS_DIR/additive"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/validate_original_logic.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --support-dir "$SUPPORT_DIR" \
  --added-marker-file "$RESULTS_DIR/additive/added_markers.csv" \
  --out-dir "$RESULTS_DIR/validation"
```

The historical public table used by these checks is read from a frozen,
checksum-verified fixture in the repository, not from a branch pointer:

- `recomputed_170_workflow/fixtures/Supplementary_File_1_historical_354rows.csv`
- SHA256 `0fdd239d438dcbab610e3bf05f2611def83f2bd6672704753914bca0a4ec06f5`
- extracted from the bioRxiv release commit `c7890d7a76ebe611f0ad6e0001d0dcf8a03bb572`

Override it with `--historical-public-supplement` if needed. The checksum is
verified on every load, so an altered or truncated fixture fails loudly rather
than silently changing the comparison.

The validation checks include:

- old Supplementary File 1 has 354 component-expanded rows;
- old evidence counts are Very Strong 1, Strong 5, Moderate 82, Weak 105,
  Weak/No Evidence 161;
- old multi-method rows and positions are 88 and 81 in the
  component-expanded frame;
- the old 157 fitted additive markers are a strict subset of the corrected 170;
- the corrected-minus-old set is exactly the 13 added markers;
- original-rule epistasis candidates are 3,542 for the old panel and 4,052 for
  the corrected panel;
- all 510 corrected-only epistasis pairs involve at least one added marker;
- raw and Galwey-adjusted p-value columns are not conflated.

## Primary recomputed 170-marker analysis

Use this workflow for the revised analysis. It recomputes effective-test counts
and permutation thresholds on the 170-marker / 4,052-interaction universe.

```bash
export RESULTS_DIR=$PROJECT_ROOT/analysis_outputs/recomputed_170_thresholds
```

### Effective-test counts

Generate the corrected epistasis interaction metadata, then compute Galwey
effective-test counts with the original R logic:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/generate_corrected_epistasis_interactions.py \
  --data-dir "$DATA_DIR" \
  --out-dir "$RESULTS_DIR/interactions"
```

The generator applies the original four-cell rule: all four parent genotype
cells must include at least 1% of isolates. It does not hash-deduplicate marker
products. The corrected 170-marker panel must produce 4,052 candidate
interactions.

```bash
Rscript analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/compute_recomputed_meff.R \
  "$DATA_DIR" \
  "$RESULTS_DIR/interactions/corrected_epistasis_interactions.csv" \
  "$RESULTS_DIR/thresholds/recomputed_meff.json"
```

The output is written to:

- `$RESULTS_DIR/thresholds/recomputed_meff.json`

### Additive mvLMM

Run additive observed chunks, then the 100 null/random-phenotype repeats, then
merge them:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/run_additive_chunk.py --help
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/merge_additive_chunks.py --help
```

Expected merged shapes:

- 170 observed additive p-value rows;
- 680 observed additive effect-size rows;
- 17,000 additive null p-value rows;
- 100 complete additive null repeats.

### Exact uvLMM

Run and merge all 170 x 6 marker-drug tests:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/run_exact_unilmm_chunk.py --help
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/merge_exact_unilmm_chunks.py --help
```

The uvLMM results are used for the display/comparison column only. They are not
counted as a fifth evidence stream.

### Epistatic mvLMM

Using the interaction file generated above, run observed chunks, then the 100
permutation repeats, then merge them:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/run_epistasis_chunk.py --help
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/run_epistasis_permutation_chunk.py --help
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/merge_epistasis_chunks.py --help
```

Expected merged shapes:

- 4,052 observed epistasis p-value rows;
- 8,104 observed epistasis effect-size rows;
- 405,200 epistasis permutation p-value rows;
- 100 complete epistasis permutation repeats.

## Threshold and evidence policy

Additive and epistasis thresholds are estimated separately.

For each null family, the scripts compute the minimum p-value within each of the
100 repeats and use the lowest of those 100 minima as the analysis threshold.
The exact threshold values are stored in `$RESULTS_DIR/thresholds/`.

Primary evidence gates:

- additive mvLMM evidence: recomputed Galwey-adjusted additive threshold and
  `Joint_effsize >= 1`;
- epistasis evidence: recomputed Galwey-adjusted epistasis threshold and
  `joint_effect_size - joint_effect_size_se >= 1`;
- exact uvLMM display support: `pv20_adj_galwey <= 0.001`;
- evidence categories remain 0 = Weak/No Evidence, 1 = Weak, 2 = Moderate,
  3 = Strong, 4 = Very Strong.

Raw `pv20` and adjusted `pv20_adj_galwey` are retained as separate columns.

## Build manuscript outputs

After all observed, null, and merge steps finish:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/build_recomputed_summary_report.py \
  --results-dir "$RESULTS_DIR" \
  --original-logic-public "$PROJECT_ROOT/analysis_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv" \
  --out-dir "$RESULTS_DIR/manuscript_outputs"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/build_corrected_evidence_table.py \
  --legacy-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --support-dir "$SUPPORT_DIR" \
  --analysis-out "$RESULTS_DIR/epistasis/merged" \
  --additive-dir "$RESULTS_DIR/additive/merged" \
  --uv-dir "$RESULTS_DIR/uvlmm/merged" \
  --epistasis-dir "$RESULTS_DIR/epistasis/merged" \
  --threshold-file "$RESULTS_DIR/thresholds/recomputed_thresholds.json" \
  --output-dir "$RESULTS_DIR/manuscript_outputs" \
  --analysis-label "Recomputed 170-marker analysis"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/build_recomputed_manuscript_summary.py \
  --results-dir "$RESULTS_DIR" \
  --historical-public "$PROJECT_ROOT/analysis_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv"

Rscript analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/generate_recomputed_supplement_figures_original_style.R \
  "$PROJECT_ROOT" "$RESULTS_DIR" \
  "$RESULTS_DIR/manuscript_outputs/supplement_figures/original_style"
```

Primary output locations:

- `$RESULTS_DIR/manuscript_outputs/Supplementary_File_1.csv`
- `$RESULTS_DIR/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
- `$RESULTS_DIR/manuscript_outputs/corrected_rerun_manuscript_summary.md`
  (the shipped results tree carries this file under its earlier name,
  `corrected_rerun_manuscript_audit.md`; the content is the same artefact)
- `$RESULTS_DIR/manuscript_outputs/recomputed_170_manuscript_statistics.md`
- `$RESULTS_DIR/manuscript_outputs/manuscript_summary/recomputed_170_manuscript_result_summary.md`
- `$RESULTS_DIR/manuscript_outputs/manuscript_summary/recomputed_170_table3_top20.csv`
- `$RESULTS_DIR/manuscript_outputs/manuscript_summary/recomputed_170_manuscript_replacement_text.md`
- `$RESULTS_DIR/manuscript_outputs/supplement_figures/original_style/`

Run the validator after building outputs:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/validate_recomputed_outputs.py \
  --results-dir "$RESULTS_DIR"
```

The validation report is written to:

- `$RESULTS_DIR/validation/recomputed_170_validation_report.md`

## Optional engine-equivalence check

This check compares the recomputed low-rank epistatic mvLMM implementation
against the original dense Limix call on selected interactions.

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/verify_epistasis_engine_equivalence.py --help
```

Outputs are written to:

- `$RESULTS_DIR/validation/engine_equivalence/`

## Committing outputs

Do not copy files into `manuscript/` automatically. Keep recomputed outputs
under `analysis_outputs/` while they are being generated and checked. If manuscript
files need to be replaced, copy only the approved small CSV/PDF outputs in a
separate commit with a clear provenance note.
