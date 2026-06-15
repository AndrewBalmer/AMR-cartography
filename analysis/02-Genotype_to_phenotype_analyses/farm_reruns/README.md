# Exact Original-Logic mvLMM/Epistasis Rebuild

This folder contains the reproducible farm workflow for the corrected 170-marker
PBP analysis. The rebuild policy is exact historical replication first: the
scripts must reproduce the old manuscript outputs and counting frame on the old
snapshot before the same implemented logic is applied to the corrected
170-marker panel.

Current corrected manuscript-facing outputs should be treated as non-final until
the golden validation passes, the 4,052-candidate epistasis universe has been
rerun, and rebuilt output tables have been reviewed.

## Roles

- GitHub branch: source of truth for code, scripts, manifests, and small final
  manuscript-facing tables/audits.
- `AMRC-repo-files/`: ignored data/results bundle. Transfer it separately with
  `rsync`; do not commit it.
- Farm repo/data directory: authoritative source for heavy analysis outputs
  after the full epistasis run starts.
- Laptop repo: lightweight checks, manuscript writing, final review.

## 1. Clone the branch on the farm

```bash
git clone https://github.com/AndrewBalmer/AMR-cartography.git
cd AMR-cartography
git checkout rerun-corrected-mvlmm-epistasis
```

## 2. Create the Python environment

Use the local farm Python/conda setup that works on your farm account. The
analysis expects `limix==3.0.4` and compatible `glimix-core`/`numpy-sugar`.

```bash
python -m venv .venv-mvlmm
source .venv-mvlmm/bin/activate
pip install -r analysis/02-Genotype_to_phenotype_analyses/farm_reruns/requirements-farm.txt
```

If the farm already has a working Limix module/conda env, use that instead and
set `PYTHON=/path/to/python` before submitting jobs.

## 3. Transfer ignored data from laptop to farm

From the laptop, run something like:

```bash
rsync -av --progress /Users/ab69/AMR-cartography/AMRC-repo-files/ \
  USER@FARM_HOST:/farm/path/AMR-cartography/AMRC-repo-files/
```

For the first run, transfer the full `AMRC-repo-files/` folder. It is about
1.3 GB locally, and transferring the full bundle avoids missing hidden
dependencies. The minimum expected files are listed in
`required_data_manifest.txt`.

## 4. Set farm environment variables

On the farm:

```bash
export PROJECT_ROOT=/farm/path/AMR-cartography
export PYTHON=$PROJECT_ROOT/.venv-mvlmm/bin/python
export FARM_OUT=$PROJECT_ROOT/farm_outputs/original_logic_rebuild
export DATA_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112
export OLD_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1
export SUPPORT_DIR="$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis"
```

Optional defaults:

```bash
export CHUNK_SIZE=25
export QUEUE=normal
export MEM_MB=8000
export N_PERMUTATIONS=100
```

## 5. Run the golden original-logic validation

Before trusting any corrected output, create the additive old-vs-new comparison
files. This writes `added_markers.csv` and reports raw-vs-adjusted p-value
counts for audit.

```bash
mkdir -p "$FARM_OUT/additive"
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/compare_additive_mvlmm.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/additive"
```

Then run the validation harness:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/validate_original_logic.py \
  --old-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --support-dir "$SUPPORT_DIR" \
  --added-marker-file "$FARM_OUT/additive/added_markers.csv" \
  --out-dir "$FARM_OUT/validation"
```

This must assert:

- old Supplementary File 1 has 354 rows with evidence counts
  VS 1 / Strong 5 / Moderate 82 / Weak 105 / Weak-No 161;
- old multi-method rows/positions are 88/81 in the component-expanded frame;
- old fitted additive markers are a strict subset of the corrected 170;
- corrected-minus-old is exactly the 13 added markers;
- old epistasis candidates are 3,542 and corrected candidates are 4,052;
- all 510 corrected-only epistasis pairs involve one of the 13 added markers;
- old adjusted-vs-raw additive and uvLMM counts differ as expected.

## 6. Generate corrected interaction metadata

Then generate the corrected interaction metadata:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/generate_corrected_epistasis_interactions.py \
  --data-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/interactions"
```

This creates:

- `corrected_epistasis_interactions.csv`
- `corrected_epistasis_interaction_summary.csv`

The generator starts from all `170 choose 2` marker pairs and keeps only pairs
where all four parent genotype cells are `>= 1%` of isolates. It does not
hash-deduplicate products. The corrected panel must produce 4,052 candidates.

## 7. Submit smoke test

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_smoke.sh
```

Check the smoke output before submitting the full run:

```bash
ls -lh "$FARM_OUT/chunks"
tail -n 50 "$FARM_OUT"/logs/smoke.*.err
tail -n 50 "$FARM_OUT"/logs/smoke.*.out
```

Expected smoke outputs:

- `epistasis_p_values.chunk_smoke.csv`
- `epistasis_effect_sizes.chunk_smoke.csv`

## 8. Submit full corrected epistasis array

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_epistasis_array.sh
```

Each LSF array task runs one chunk of interaction tests. Existing chunk files are
skipped unless `--force` is passed manually to the runner.

## 9. Submit corrected permutation arrays

```bash
export N_PERMUTATIONS=100
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_permutation_array.sh
```

Permutation outputs are retained for audit/sensitivity checks. Exact manuscript
outputs use the historical adjusted epistasis threshold `0.0007620121` with
Galwey meff `39`, matching the original implemented manuscript logic.

## 10. Merge farm outputs

After the observed and permutation jobs finish:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_merge.sh
```

Or run directly:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_epistasis_chunks.py \
  --chunk-dir "$FARM_OUT/chunks" \
  --out-dir "$FARM_OUT/merged"
```

Merged outputs include:

- `corrected_epistasis_p_values.csv`
- `corrected_epistasis_effect_sizes.csv`
- `corrected_epistasis_permutation_p_values.csv`
- `corrected_epistasis_permutation_minima.csv`
- `corrected_epistasis_marker_support.csv`
- `corrected_epistasis_summary.json`

Epistasis support requires both:

- `pv20_adj_galwey <= 0.0007620121`;
- `joint_effect_size - joint_effect_size_se >= 1`.

The merge script asserts 4,052 observed interactions by default. If it sees the
older 4,563 hash-deduplicated universe, it should fail rather than build final
tables.

## 11. Rebuild review outputs, not manuscript/ directly

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_exact_unilmm.py \
  --legacy-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/merged"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/build_corrected_evidence_table.py \
  --legacy-dir "$OLD_DIR" \
  --new-dir "$DATA_DIR" \
  --support-dir "$SUPPORT_DIR" \
  --analysis-out "$FARM_OUT/merged" \
  --output-dir "$FARM_OUT/manuscript_outputs"
```

This writes:

- `Supplementary_File_1.csv`: the legacy manuscript/thesis
  component-level supplementary table, rebuilt with corrected additive mvLMM,
  exact uvLMM, and corrected epistasis evidence.
- `Supplementary_File_1_corrected_marker_level.csv`: the 170-row
  corrected marker-level audit table.
- `corrected_rerun_manuscript_audit.md`: row counts, evidence
  counts, thresholds, and corrected epistasis merge diagnostics.

Only after review should these files be deliberately copied over
`manuscript/` outputs in one provenance-noted commit.

## 12. Copy small final outputs back if needed

Do not copy the whole farm output directory back into Git by default. Copy only
the final small manuscript-facing outputs and any summary files needed for
review.

```bash
rsync -av USER@FARM_HOST:/farm/path/AMR-cartography/farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv \
  /Users/ab69/AMR-cartography/manuscript/
```

## Notes

- Additive mvLMM evidence uses historical adjusted p-value logic:
  `pv20_adj_galwey <= 0.000588` and `Joint_effsize >= 1`.
- Exact uvLMM display counts use `pv20_adj_galwey <= 0.001` for all 170
  markers: 157 historical marker results plus exact reruns for the 13
  corrected-panel additions.
- Raw `pv20` and adjusted `pv20_adj_galwey` are kept as separate columns.
- Once the full farm run starts, do not regenerate heavy analysis outputs on the
  laptop. The farm outputs become authoritative.
