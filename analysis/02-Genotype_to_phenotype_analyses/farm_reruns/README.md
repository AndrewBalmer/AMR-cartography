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

## Head Node Safety

Do not run analysis-scale Python or R commands directly on `farm22-head*` or
`farm22-pam-01`. The farm head nodes are only for editing, `git`, `bsub`,
`bjobs`, small log checks, and file transfer. Model fits, full merges over chunk
directories, large CSV scans, row-count sweeps across permutation chunks, and
missing-chunk reruns must be submitted through LSF.

In particular, do not run these directly on a head node:

- `run_additive_chunk.py`
- `run_exact_unilmm_chunk.py`
- `run_epistasis_chunk.py`
- `run_epistasis_permutation_chunk.py`
- `merge_additive_chunks.py`
- `merge_epistasis_chunks.py`
- ad hoc Python loops over all `epistasis_perm_p_values.chunk_*.csv`
- full R/Rmd manuscript rebuilds or any script that reads all permutation chunks

Use the `lsf/submit_*.sh` scripts below. If a single chunk needs debugging, use
`bsub` or an interactive LSF session rather than running the Python command
directly on the login/head node. This matters even for "just one chunk": a
pathological LMM fit can burn CPU on the head node long enough to trigger
Arbiter penalties.

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
files. This step is lightweight, but if the head node is under Arbiter penalty
or the data filesystem is slow, submit it via LSF rather than running it
directly.

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

Do not run `merge_epistasis_chunks.py` directly on a head node. It reads all
observed, effect-size, and permutation chunk files and can trigger head-node CPU
or I/O penalties.

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
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_review_outputs.sh
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

## Missing Or Stalled Permutation Chunks

If permutation chunks stall or are missing, do not rerun them directly on a head
node. Use the missing-chunk LSF submitter:

```bash
export ON_FIT_ERROR=pvalue-one
export FIT_TIMEOUT_SECONDS=300
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_missing_permutation_chunks.sh
```

This scans for missing expected permutation chunk files, writes the missing list
to `$FARM_OUT/missing_permutation_chunks.tsv`, and submits one LSF job per
missing chunk. With `ON_FIT_ERROR=pvalue-one` and `FIT_TIMEOUT_SECONDS=300`, a
pathological fit becomes an explicit conservative p=1 permutation row instead of
hanging indefinitely. If a chunk still cannot complete, record any synthetic
conservative fallback in `$FARM_OUT/synthetic_conservative_permutation_chunks.tsv`
before merging.

## 12. Copy small final outputs back if needed

Do not copy the whole farm output directory back into Git by default. Copy only
the final small manuscript-facing outputs and any summary files needed for
review.

```bash
rsync -av USER@FARM_HOST:/farm/path/AMR-cartography/farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv \
  /Users/ab69/AMR-cartography/manuscript/
```

## Recomputed 170-marker primary analysis

Use this workflow for a revised manuscript that claims the analysis was rerun on
the corrected 170-marker universe with recomputed thresholds. It is separate
from `original_logic_rebuild/`, which remains the historical-threshold
comparability audit.

Set the recomputed output root:

```bash
export PROJECT_ROOT=/farm/path/AMR-cartography
export PYTHON=$PROJECT_ROOT/.venv-mvlmm/bin/python
export FARM_OUT=$PROJECT_ROOT/farm_outputs/recomputed_170_thresholds
export DATA_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1-additive-production-20260507-150112
export OLD_DIR=$PROJECT_ROOT/AMRC-repo-files/pythonProject1
export SUPPORT_DIR="$PROJECT_ROOT/AMRC-repo-files/Streptococcus pneumoniae analysis"
export N_PERMUTATIONS=100
export MAX_CONCURRENT=200
```

Submit jobs through LSF only:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_install_poolr.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_meff.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_smoke.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_permutation_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_uvlmm_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_smoke.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_array.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_permutation_array.sh
```

`submit_recomputed_install_poolr.sh` installs or verifies R `poolr` in the user
R library. It is required for the Galwey effective-test calculation and should
be allowed to finish before submitting `submit_recomputed_meff.sh`. To let LSF
enforce that ordering, submit the meff step with a dependency on the install job,
for example:

```bash
DEPENDENCY='done(123456)' \
  bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_meff.sh
```

After all observed and permutation jobs finish:

```bash
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_additive_merge.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_uvlmm_merge.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_epistasis_merge.sh
bash analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/submit_recomputed_review_outputs.sh
```

Primary recomputed thresholds use the original lowest-min-p policy: compute the
minimum p-value within each of the 100 null/permutation repeats, then use the
lowest of those 100 minima as the analysis threshold. The scripts store both raw
and Galwey-adjusted threshold equivalents in `$FARM_OUT/thresholds/`.

The array wrappers use LSF's `%` concurrency limit. `MAX_CONCURRENT` defaults to
100 for observed arrays and 200 for permutation arrays; lower it if the cluster
is busy or IT asks for gentler scheduling. The permutation wrappers submit one
combined LSF array and map each task to a `(permutation, chunk)` pair, rather
than submitting 100 separate arrays.

Primary recomputed outputs are valid only if
`$FARM_OUT/validation/recomputed_170_validation_report.md` reports all checks
passing. The validator requires complete additive, uvLMM, and epistasis outputs,
zero synthetic chunks, and zero non-ok fit rows.

## Notes

- `original_logic_rebuild/` uses the historical literal thresholds for
  comparability with the implemented manuscript pipeline.
- `recomputed_170_thresholds/` uses recomputed 170-marker thresholds. Additive
  evidence requires the recomputed Galwey-adjusted mvLMM threshold and
  `Joint_effsize >= 1`; epistasis evidence requires the recomputed
  Galwey-adjusted epistasis threshold and
  `joint_effect_size - joint_effect_size_se >= 1`.
- Exact uvLMM display counts use `pv20_adj_galwey <= 0.001` for all 170
  markers. uvLMM remains a display/comparison column and is not counted as a
  fifth evidence stream.
- Raw `pv20` and adjusted `pv20_adj_galwey` are kept as separate columns.
- Once the full farm run starts, do not regenerate heavy analysis outputs on the
  laptop. The farm outputs become authoritative.
