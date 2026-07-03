# Recomputed 170-marker mvLMM/epistasis workflow

This directory contains the reproducible farm workflow for the revised
*Streptococcus pneumoniae* PBP genotype-to-phenotype analysis.

The workflow has two roles:

1. Original-logic validation: reproduce the original 157-marker manuscript
   output frame and evidence counts as a quality-control check.
2. Primary recomputed analysis: rerun the additive mvLMM, exact uvLMM,
   epistatic mvLMM, permutation thresholds, manuscript tables, and supporting
   figures on the corrected 170-marker panel.

Older R/Rmd scripts elsewhere in `analysis/` are retained as historical
manuscript-generation scripts. The corrected 170-marker farm workflow is
implemented here.

## Head-node safety

Do not run analysis-scale Python or R commands directly on farm head nodes. Use
the LSF submitters in `lsf/` for model fits, permutation jobs, large merges, and
figure rebuilds that read many result files.

Suitable head-node tasks are editing, `git`, `bsub`, `bjobs`, small log checks,
and file transfer. Submit these through LSF:

- `run_additive_chunk.py`
- `run_exact_unilmm_chunk.py`
- `run_epistasis_chunk.py`
- `run_epistasis_permutation_chunk.py`
- `merge_additive_chunks.py`
- `merge_epistasis_chunks.py`
- full R figure rebuilds
- any ad hoc loop over all chunk or permutation outputs

For a single problematic chunk, use an interactive LSF session or a small LSF
job rather than running the model fit on a head node.

## Inputs and environment

Large input bundles and farm outputs are not committed. The expected data files
are listed in `required_data_manifest.txt`.

Set these variables before running the submitters:

```bash
export PROJECT_ROOT=/path/to/AMR-cartography
export PYTHON=$PROJECT_ROOT/.venv-mvlmm/bin/python
export DATA_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112
export OLD_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1
export SUPPORT_DIR="$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis"
export FARM_OUT=$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds
export N_PERMUTATIONS=100
```

Optional LSF settings:

```bash
export QUEUE=normal
export MEM_MB=8000
export CHUNK_SIZE=25
export MAX_CONCURRENT=100
```

Create a Python environment with the farm-compatible Limix stack:

```bash
python3 -m venv .venv-mvlmm
source .venv-mvlmm/bin/activate
pip install -r analysis/02-Genotype_to_phenotype_analyses/farm_reruns/requirements-farm.txt
```

If a suitable farm module or conda environment already exists, set `PYTHON` to
that interpreter instead.

## Original-logic validation

Run this quality-control workflow before using revised manuscript outputs. It
checks that the code reproduces the original manuscript counting frame on the
old snapshot, then verifies that the corrected 170-marker panel contains the old
157 markers plus the 13 added markers.

```bash
export FARM_OUT=$PROJECT_ROOT/farm_outputs/original_logic_rebuild

$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/compare_additive_mvlmm.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/additive"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/validate_original_logic.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --support-dir "$SUPPORT_DIR" \
  --added-marker-file "$FARM_OUT/additive/added_markers.csv" \
  --out-dir "$FARM_OUT/validation"
```

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
export FARM_OUT=$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds
```

### Effective-test counts

Compute Galwey effective-test counts with the original R logic:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_install_poolr.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_meff.sh
```

The output is written to:

- `$FARM_OUT/thresholds/recomputed_meff.json`

### Additive mvLMM

Submit the additive smoke test, observed chunks, and 100 null/permutation
repeats:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_smoke.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_permutation_array.sh
```

Merge after all jobs finish:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_merge.sh
```

Expected merged shapes:

- 170 observed additive p-value rows;
- 680 observed additive effect-size rows;
- 17,000 additive null p-value rows;
- 100 complete additive null repeats.

### Exact uvLMM

Run and merge all 170 x 6 marker-drug tests:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_uvlmm_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_uvlmm_merge.sh
```

The uvLMM results are used for the display/comparison column only. They are not
counted as a fifth evidence stream.

### Epistatic mvLMM

Generate corrected interaction metadata:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/generate_corrected_epistasis_interactions.py \
  --data-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/interactions"
```

The generator applies the original four-cell rule: all four parent genotype
cells must include at least 1% of isolates. It does not hash-deduplicate marker
products. The corrected 170-marker panel must produce 4,052 candidate
interactions.

Submit the smoke test, observed chunks, and 100 permutation repeats:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_smoke.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_permutation_array.sh
```

Merge after all jobs finish:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_merge.sh
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
The exact threshold values are stored in `$FARM_OUT/thresholds/`.

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

After all observed, null, and merge jobs finish:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_review_outputs.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_supplement_figures_original_style.sh
```

Primary output locations:

- `$FARM_OUT/manuscript_outputs/Supplementary_File_1.csv`
- `$FARM_OUT/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
- `$FARM_OUT/manuscript_outputs/corrected_rerun_manuscript_summary.md`
- `$FARM_OUT/manuscript_outputs/recomputed_170_manuscript_statistics.md`
- `$FARM_OUT/manuscript_outputs/manuscript_summary/recomputed_170_manuscript_result_summary.md`
- `$FARM_OUT/manuscript_outputs/manuscript_summary/recomputed_170_table3_top20.csv`
- `$FARM_OUT/manuscript_outputs/manuscript_summary/recomputed_170_manuscript_replacement_text.md`
- `$FARM_OUT/manuscript_outputs/supplement_figures/original_style/`

Run the validator after building outputs:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/validate_recomputed_outputs.py \
  --farm-out "$FARM_OUT"
```

The validation report is written to:

- `$FARM_OUT/validation/recomputed_170_validation_report.md`

## Structure context outputs

Structural annotations are generated separately from the evidence scoring and
are used only for biological interpretation.

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/annotate_structure_context.py \
  --farm-out "$FARM_OUT"
```

Outputs are written under:

- `$FARM_OUT/manuscript_outputs/structure_context/`

## Optional engine-equivalence check

This LSF job compares the recomputed low-rank epistatic mvLMM implementation
against the original dense Limix call on selected interactions.

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_verify_epistasis_engine.sh
```

Outputs are written to:

- `$FARM_OUT/validation/engine_equivalence/`

## Committing outputs

Do not copy files into `manuscript/` automatically. Keep recomputed outputs
under `farm_outputs/` while they are being generated and checked. If manuscript
files need to be replaced, copy only the approved small CSV/PDF outputs in a
separate commit with a clear provenance note.
