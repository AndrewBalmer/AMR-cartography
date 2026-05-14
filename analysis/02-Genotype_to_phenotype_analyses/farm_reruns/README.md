# Corrected mvLMM/Epistasis Farm Rerun

This folder contains the reproducible farm workflow for the corrected 170-marker
PBP analysis. The goal is to keep `main` as the previous manuscript-analysis
baseline, while this branch runs the corrected additive/uvLMM/epistasis analyses
on the farm.

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
export FARM_OUT=$PROJECT_ROOT/farm_outputs/corrected_epistasis
export DATA_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112
```

Optional defaults:

```bash
export CHUNK_SIZE=25
export QUEUE=normal
export MEM_MB=8000
export N_PERMUTATIONS=100
```

## 5. Generate corrected interaction metadata

First create the additive old-vs-new comparison files. This also writes
`added_markers.csv`, which is useful if the exact added-marker uvLMM ever needs
to be rerun.

```bash
mkdir -p "$FARM_OUT/additive"
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/compare_additive_mvlmm.py \
  --old-dir "$PROJECT_ROOT/AMRC-repo-files/pythonProject1" \
  --new-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/additive"
```

Then generate the corrected interaction metadata:

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/generate_corrected_epistasis_interactions.py \
  --data-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/interactions"
```

This creates:

- `corrected_epistasis_interactions.csv`
- `corrected_epistasis_interaction_summary.csv`

The generator starts from all `170 choose 2` marker pairs and removes constant,
parent-identical, and duplicate interactions.

## 6. Submit smoke test

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

## 7. Submit full corrected epistasis array

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_epistasis_array.sh
```

Each LSF array task runs one chunk of interaction tests. Existing chunk files are
skipped unless `--force` is passed manually to the runner.

## 8. Submit corrected permutation arrays

```bash
export N_PERMUTATIONS=100
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_permutation_array.sh
```

The corrected epistasis threshold is calculated from corrected-panel
permutations, not from the previous manuscript threshold.

## 9. Merge farm outputs

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

## 10. Rebuild final manuscript-facing table

```bash
$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_exact_unilmm.py \
  --legacy-dir "$PROJECT_ROOT/AMRC-repo-files/pythonProject1" \
  --new-dir "$DATA_DIR" \
  --out-dir "$FARM_OUT/merged"

$PYTHON analysis/02-Genotype_to_phenotype_analyses/farm_reruns/build_corrected_evidence_table.py \
  --legacy-dir "$PROJECT_ROOT/AMRC-repo-files/pythonProject1" \
  --new-dir "$DATA_DIR" \
  --support-dir "$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis" \
  --analysis-out "$FARM_OUT/merged" \
  --manuscript-dir "$PROJECT_ROOT/manuscript"
```

This writes:

- `manuscript/Supplementary_File_1.csv`
- `manuscript/Supplementary_File_1_corrected_marker_level.csv`
- `manuscript/corrected_rerun_manuscript_audit.md`

## 11. Copy small final outputs back if needed

Do not copy the whole farm output directory back into Git by default. Copy only
the final small manuscript-facing outputs and any summary files needed for
review.

```bash
rsync -av USER@FARM_HOST:/farm/path/AMR-cartography/manuscript/Supplementary_File_1.csv \
  /Users/ab69/AMR-cartography/manuscript/
```

## Notes

- The additive mvLMM threshold is fixed at `pv20 < 0.0009078488974311251`.
- Exact uvLMM evidence is available for all 170 markers: 157 historical marker
  results plus exact reruns for the 13 corrected-panel additions.
- Once the full farm run starts, do not regenerate heavy analysis outputs on the
  laptop. The farm outputs become authoritative.
