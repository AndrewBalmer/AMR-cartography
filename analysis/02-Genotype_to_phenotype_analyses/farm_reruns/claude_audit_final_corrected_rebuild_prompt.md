# Claude Prompt: Audit Final Corrected Original-Logic Rebuild

You are auditing the `AMR-cartography` repo after a corrected mvLMM/epistasis rerun. Work surgically and read-only unless explicitly asked otherwise. Do not modify code, manuscript files, or outputs. Produce a detailed audit report that checks whether the corrected branch faithfully implements the original manuscript logic, with only the intended addition of the validated 170-marker panel and rerun epistasis outputs.

The user is worried about silent changes in scientific logic, p-value thresholds, counting frames, manuscript-facing numbers, and whether any result changes are real. Your job is to independently verify Codex's work and identify any mistakes, inconsistencies, stale outputs, changed assumptions, or places where the new output should not be trusted.

## Branch And Scope

- Repo path: `/nfs/users/nfs_a/ab69/AMR-cartography`
- Current branch should be: `rerun-corrected-mvlmm-epistasis`
- This must not modify or overwrite `main`.
- `main` is the old/preprint/manuscript baseline.
- Compare current branch against `main` and against `origin/rerun-corrected-mvlmm-epistasis` where useful.
- Audit the whole branch relative to `main`, not only the latest local changes. Run this first and treat it as the source of truth:

```bash
git log --oneline main..HEAD
```

The scientific/workflow branch history relative to `main` is expected to include these commits, newest first, plus possibly a later documentation-only commit that adds this audit prompt:
  - `03f756e Document head-node-safe farm rerun workflow`
  - `5f60c79 Allow timed conservative permutation reruns`
  - `6ac2d87 Implement original-logic corrected rerun workflow`
  - `e0d68f5 Count non-ok fit rows in epistasis merge`
  - `4fad94f Add fit timeout for permutation repairs`
  - `a0fc176 Tolerate numerical failures in permutation cleanup`
  - `758ea36 Exclude smoke chunks from epistasis merge`
  - `a445481 Use configured Python in LSF array wrappers`
  - `0393e0b Ignore farm data symlinks`
  - `89d77fc Add farm rerun workflow for corrected epistasis`
- Current branch is expected to be ahead of origin by these local-only commits:
  - `6ac2d87 Implement original-logic corrected rerun workflow`
  - `5f60c79 Allow timed conservative permutation reruns`
  - `03f756e Document head-node-safe farm rerun workflow`

Run lightweight checks on the head node only. Do not run model fits, all-chunk scans, or full merges directly on the head node. If you need heavy processing, write the command that should be submitted through LSF instead of running it.

## Important Output Locations

Corrected review outputs are intentionally not in `manuscript/`. They are here:

- `farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv`
- `farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
- `farm_outputs/original_logic_rebuild/manuscript_outputs/corrected_rerun_manuscript_audit.md`
- `farm_outputs/original_logic_rebuild/validation/original_logic_validation_report.md`
- `farm_outputs/original_logic_rebuild/merged/corrected_epistasis_summary.json`
- `farm_outputs/original_logic_rebuild/synthetic_conservative_permutation_chunks.tsv`

Old baseline public table should be read from:

```bash
git show main:manuscript/Supplementary_File_1.csv
```

Do not compare marker-level counts to component-expanded manuscript counts as if they are the same frame.

## Source Changes To Audit

Please review all files changed on the branch relative to `main`:

```bash
git diff --name-status main..HEAD
```

Key intended changes:

### `.gitignore`

- Adds `.venv*/` so local Python environments do not appear as dirty files.

### `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/common.py`

Intended changes:

- Defines historical thresholds/constants:
  - additive mvLMM evidence threshold: `HISTORICAL_ADDITIVE_THRESHOLD = 0.000588`
  - uvLMM display threshold: `HISTORICAL_UV_THRESHOLD = 0.001`
  - epistasis threshold: `HISTORICAL_EPISTASIS_THRESHOLD = 0.0007620121`
  - additive/uvLMM Galwey meff: `28.0`
  - epistasis Galwey meff: `39.0`
  - epistasis lower-bound effect filter: `joint_effect_size - joint_effect_size_se >= 1`
- Adds `add_galwey_adjusted_pvalues`, preserving raw `pv20` separately.
- Adds `four_cell_interaction_metadata`, which should implement the original R `model.matrix(Y ~ .^2)` pairwise interaction rule with all four genotype cells `>= 1%` of isolates.
- Does not hash-deduplicate interaction products for exact manuscript replication.

Audit questions:

- Are raw and adjusted p-value columns always distinct?
- Is the four-cell rule exactly the original historical rule?
- Is the `>= 1%` comparison correct, given `n=3620` and one percent `36.2`, so minimum integer count is `37`?

### `generate_corrected_epistasis_interactions.py`

Intended behavior:

- Reads `S.pneumo_map_dummy_gen_test_markers.csv` from the corrected data dir.
- Generates all `170 choose 2` parent marker pairs.
- Keeps only pairs where all four genotype cells pass the original 1% rule.
- Does not remove duplicates/complements/hash-identical products.
- Asserts default expected count `4052`.
- Writes:
  - `corrected_epistasis_interactions.csv`
  - `corrected_epistasis_interaction_summary.csv`
  - `corrected_epistasis_interaction_summary.json`

Audit questions:

- Does it really produce 4,052 interactions on the corrected 170-marker panel?
- Does it reproduce old 3,542 interactions when run/compared against the old epistasis universe?
- Are the old 3,542 unordered pairs a strict subset of the corrected 4,052?
- Are the 510 corrected-only pairs all involving one of the 13 added markers?

### `compare_additive_mvlmm.py`

Intended behavior:

- Compares old and corrected additive mvLMM marker panels.
- Uses Galwey-adjusted p-values, not raw p-values, for threshold counts.
- Reports raw-vs-adjusted counts separately.
- Writes `added_markers.csv`.

Audit questions:

- Does it show old `157`, corrected `170`, shared `157`, added `13`, dropped `0`?
- Does `added_markers.csv` exactly match corrected-minus-old fitted markers?
- Does it avoid using raw p-values as adjusted p-values?

### `merge_exact_unilmm.py`

Intended behavior:

- Combines old 157 uvLMM tests with exact reruns for the 13 added markers.
- Writes `pv20` and `pv20_adj_galwey` separately.
- Counts significant uvLMM drugs using `pv20_adj_galwey <= 0.001`.
- Expected marker-drug rows: `170 * 6 = 1020`.
- Expected support rows: `170`.

Audit questions:

- Are uvLMM support counts adjusted-p based?
- Are exact added-marker effect sizes merged correctly if present?
- Are there any old raw-p or `< 0.001` raw comparisons still driving display counts?

### `merge_epistasis_chunks.py`

Intended behavior:

- Fails unless observed rows equal `4052`.
- Merges observed p-values and effect sizes.
- Computes `pv20_adj_galwey = pv20 * 39`.
- Computes `joint_effect_size = sqrt(D1^2 + D2^2)`.
- Computes `joint_effect_size_se = sqrt(SE_D1^2 + SE_D2^2)`.
- Defines support as:
  - `pv20_adj_galwey <= 0.0007620121`
  - and `joint_effect_size - joint_effect_size_se >= 1`
- Pivots supported interactions to both parent markers and counts marker-level epistasis support.
- Writes:
  - `corrected_epistasis_p_values.csv`
  - `corrected_epistasis_effect_sizes.csv`
  - `corrected_epistasis_permutation_p_values.csv`
  - `corrected_epistasis_permutation_minima.csv`
  - `corrected_epistasis_supported_interactions.csv`
  - `corrected_epistasis_marker_support.csv`
  - `corrected_epistasis_summary.json`

Actual merged summary reported:

```json
{
  "support_rule": "pv20_adj_galwey <= threshold AND joint_effect_size - joint_effect_size_se >= lower_bound",
  "galwey_meff": 39.0,
  "historical_epistasis_threshold": 0.0007620121,
  "support_effect_lower_bound": 1.0,
  "observed_interactions": 4052,
  "effect_rows": 8104,
  "permutation_rows": 405200,
  "permutations": 100,
  "raw_permutation_threshold": 0.0002065143234331,
  "adjusted_permutation_threshold": 0.0080540586138909,
  "observed_non_ok_fit_rows": 0,
  "permutation_non_ok_fit_rows": 25,
  "alpha": 0.05,
  "pvalue_threshold_only_interactions": 3093,
  "support_interactions": 1924,
  "markers_with_epistasis_support": 131
}
```

Audit questions:

- Is the support table based on p-value plus lower-bound effect, not p-value alone?
- Are both parents counted exactly once per supported interaction?
- Are the permutation thresholds only audit/sensitivity outputs, not manuscript evidence thresholds?
- Do the 25 non-ok permutation rows all come from the documented conservative fallback?
- Could the 25 p=1 rows affect any manuscript-facing results? They should not, because exact manuscript outputs use the historical literal epistasis threshold, not the recomputed permutation threshold.

### `build_corrected_evidence_table.py`

Intended behavior:

- Writes review outputs to `farm_outputs/original_logic_rebuild/manuscript_outputs` by default, not `manuscript/`.
- Uses additive evidence:
  - `pv20_adj_galwey <= 0.000588`
  - and `joint_effect_size >= 1`
- Uses uvLMM display counts:
  - `pv20_adj_galwey <= 0.001`
- Uses epistasis evidence from `corrected_epistasis_marker_support.csv`, already filtered by p-value and lower-bound effect.
- Builds two frames:
  - component-expanded legacy/public Supplementary File 1
  - marker-level audit table
- Should preserve the original manuscript evidence categories:
  - `0 = Weak/No Evidence`
  - `1 = Weak`
  - `2 = Moderate`
  - `3 = Strong`
  - `4 = Very Strong`

Audit questions:

- Does it preserve the original public/component-expanded counting frame?
- Does it avoid labeling raw `pv20` as `pv20_adj_galwey`?
- Does the public `Adj. p-value` column contain Galwey-adjusted additive p-values?
- Does it avoid silently copying outputs into `manuscript/`?
- Does it reproduce old 354-row counts on old data, if run in old-data mode or through the validation harness?

### `validate_original_logic.py`

Intended behavior:

- Golden validation harness.
- Should pass before corrected results are trusted.
- Validates old `main` Supplementary File 1:
  - rows `354`
  - evidence counts: Very Strong `1`, Strong `5`, Moderate `82`, Weak `105`, Weak/No Evidence `161`
  - multi-method rows `88`
  - multi-method positions `81`
  - zero slash-containing public substitution rows
- Validates marker panels:
  - old fitted additive markers `157`
  - corrected fitted additive markers `170`
  - old is a strict subset of corrected
  - corrected-minus-old is `13`
  - each added marker is binary, non-missing, non-duplicate, non-complement, and passes `Present > 1%` / `Absent > 1%`
- Validates epistasis candidates:
  - old count table rows `3542`
  - old epistasis matrix columns `3542`
  - old epistasis p-value rows `3542`
  - corrected original-rule candidates `4052`
  - corrected-only pairs `510`
  - all corrected-only pairs involve added markers
- Validates adjusted-vs-raw scoring:
  - old additive adjusted evidence hits `78`
  - old additive raw-threshold hits `84`
  - old uvLMM raw marker-drug tests `127`
  - old uvLMM adjusted marker-drug tests `66`
  - old uvLMM raw significant markers `61`
  - old uvLMM adjusted significant markers `40`

Audit questions:

- Are all assertions correct and meaningful?
- Is it accidentally relying on changed branch outputs rather than old `main` where it should use old data?
- Does it use component-expanded manuscript frame for manuscript claims?

### LSF scripts and head-node safety

Changed or added:

- `lsf/submit_smoke.sh`
- `lsf/submit_epistasis_array.sh`
- `lsf/submit_permutation_array.sh`
- `lsf/submit_merge.sh`
- `lsf/submit_review_outputs.sh`
- `lsf/submit_missing_permutation_chunks.sh`

Intended behavior:

- Defaults should use `farm_outputs/original_logic_rebuild`.
- Model fits, full merges, review-output rebuilds, and missing permutation chunks should be submitted via LSF.
- `submit_permutation_array.sh` accepts:
  - `ON_FIT_ERROR`
  - `FIT_TIMEOUT_SECONDS`
- `submit_missing_permutation_chunks.sh` scans missing chunks and submits one LSF job per missing permutation chunk, with conservative fallback controls.
- README should clearly warn not to run heavy Python/R commands directly on `farm22-head*`.

Audit questions:

- Are the documented commands head-node safe?
- Are there any remaining instructions telling users to run heavy all-chunk operations directly?
- Are the new LSF scripts robust to paths with spaces, especially `SUPPORT_DIR`?
- Is quoting safe in the `bsub` command strings?

## Actual Workflow Run By Codex

Codex ran the workflow on the rerun branch, with outputs under `farm_outputs/original_logic_rebuild`.

High-level sequence:

1. Tidied dirty files:
   - Restored tracked `manuscript/Supplementary_File_1.csv`.
   - Moved non-final dirty manuscript outputs to:
     - `farm_outputs/original_logic_rebuild/tidied_invalid_manuscript_outputs/Supplementary_File_1.pre_tidy.csv`
     - `farm_outputs/original_logic_rebuild/tidied_invalid_manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
     - `farm_outputs/original_logic_rebuild/tidied_invalid_manuscript_outputs/corrected_rerun_manuscript_audit.md`
2. Committed source changes:
   - `6ac2d87`
   - `5f60c79`
   - `03f756e`
3. Generated corrected interaction metadata:
   - `4052` interactions
   - `170` additive markers
   - `14365` possible pairs
   - `10313` failed four-cell filter
   - one-percent count `36.2`
   - min retained cell count `37`
4. Ran smoke test via LSF:
   - 5 fits
   - all `fit_status=ok`
5. Ran full observed epistasis via LSF:
   - 163 chunks
   - `4052` rows
   - all `fit_status=ok`
6. Ran 100 permutation arrays via LSF:
   - expected `16300` permutation chunk files
   - expected `405200` permutation rows
7. Some permutation chunks stalled.
   - 131 missing chunks were resubmitted with timeout/fallback controls.
   - One chunk still wedged: permutation `9`, array `146`.
   - Codex wrote one documented synthetic conservative p=1 chunk with 25 rows:
     - `farm_outputs/original_logic_rebuild/synthetic_conservative_permutation_chunks.tsv`
     - reason: LSF rerun wedged and direct timeout rerun exceeded outer timeout
8. Merged epistasis outputs:
   - observed rows `4052`
   - effect rows `8104`
   - permutation rows `405200`
   - permutations `100`
   - observed non-ok rows `0`
   - permutation non-ok rows `25`
9. Rebuilt uvLMM support:
   - `1020` marker-drug rows
   - `170` marker support rows
10. Built review manuscript outputs:
   - `farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv`
   - `farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
   - `farm_outputs/original_logic_rebuild/manuscript_outputs/corrected_rerun_manuscript_audit.md`
11. Ran validation harness against rebuilt outputs:
   - passed.
12. No outputs were deliberately copied into `manuscript/`.

## Headline Calculations To Recheck

All manuscript headline comparisons below use the original public/component-expanded Supplementary File 1 frame unless explicitly stated otherwise.

Old baseline from `main:manuscript/Supplementary_File_1.csv`:

- Public rows: `354`
- Evidence counts:
  - Very Strong `1`
  - Strong `5`
  - Moderate `82`
  - Weak `105`
  - Weak/No Evidence `161`
- Any evidence (`Very Strong`, `Strong`, `Moderate`, `Weak`):
  - rows `193`
  - unique component substitutions `172`
  - unique PBP positions `147`
  - proportion of 285 variable positions: `147 / 285 = 51.6%`
  - positions by PBP:
    - 1A `50`
    - 2B `31`
    - 2X `66`
- Multi-method evidence (`Very Strong`, `Strong`, `Moderate`):
  - rows `88`
  - unique component substitutions `88`
  - unique PBP positions `81`
  - positions by PBP:
    - 1A `18`
    - 2B `24`
    - 2X `39`
  - component rows by PBP:
    - 1A `20`
    - 2B `26`
    - 2X `42`
  - effect size range: `0.330–3.079`
  - additive p-threshold plus epistasis: `79 / 88 = 89.8%`
  - only Very Strong row: PBP2X `I371T`, beta `3.079`, sig interactions `36`, mv/uv `Yes/Yes (6)`

Corrected rebuilt public/component-expanded output:

- Public rows: `394`
- Evidence counts:
  - Very Strong `1`
  - Strong `5`
  - Moderate `109`
  - Weak `121`
  - Weak/No Evidence `158`
- Any evidence:
  - rows `236`
  - unique component substitutions `215`
  - unique PBP positions `179`
  - proportion of 285 variable positions: `179 / 285 = 62.8%`
  - positions by PBP:
    - 1A `60`
    - 2B `50`
    - 2X `69`
- Multi-method evidence:
  - rows `115`
  - unique component substitutions `115`
  - unique PBP positions `106`
  - positions by PBP:
    - 1A `24`
    - 2B `42`
    - 2X `40`
  - component rows by PBP:
    - 1A `26`
    - 2B `46`
    - 2X `43`
  - effect size range: `0.374–3.187`
  - additive p-threshold plus epistasis: `106 / 115 = 92.2%`
  - only Very Strong row: PBP2X `I371T`, beta `3.187`, sig interactions `40`, mv/uv `Yes/Yes (6)`

Corrected marker-level audit frame:

- Marker rows: `170`
- Evidence counts:
  - Very Strong `1`
  - Strong `5`
  - Moderate `78`
  - Weak `68`
  - Weak/No Evidence `18`
- Markers with `>=2` evidence streams: `84`
- `>=2` by PBP:
  - PBP1A `18`
  - PBP2B `25`
  - PBP2X `41`
- Evidence stream totals:
  - single-sub evidence `1`
  - cluster evidence `21`
  - mvLMM evidence `90`
  - epistasis evidence `131`

Set changes in multi-method component substitutions:

- Gained unique multi-method component substitutions: `31`
- Of these, `30/31` are components of the 13 added markers.
- One gained component is not itself part of the 13 added markers: `PBP2X A347S`
- Lost unique multi-method component substitutions: `4`
  - `1A A550P`
  - `1A N546G`
  - `2B D578E`
  - `2B T569K`
- The four lost rows remain Weak in the rebuilt output because beta is now just under 1:
  - `1A A550P`: beta `0.993`, adj p `3.038e-11`, sig interactions `33`, total `1`
  - `1A N546G`: beta `0.993`, adj p `3.038e-11`, sig interactions `33`, total `1`
  - `2B D578E`: beta `0.964`, adj p `<1e-16`, sig interactions `57`, total `1`
  - `2B T569K`: beta `0.964`, adj p `<1e-16`, sig interactions `57`, total `1`

Gained unique multi-method component substitutions:

```text
1A E397V
1A H503N
1A L421I
1A N517D
1A T392S
1A T495I
1A V408L
1A Y497H
2B A533G
2B A619G
2B D555E
2B D625G
2B D641E
2B G597A
2B G597P
2B H585Q
2B K556Q
2B L566V
2B L609T
2B M581V
2B N606D
2B Q565S
2B Q567E
2B S473T
2B S480A
2B S640T
2B T595G
2B T630N
2B V542I
2B V574I
2X A347S
```

## Required Audit Tasks

Please produce a report with findings first, ordered by severity. Include exact file/line references where possible. At minimum, check:

1. **Branch safety**
   - Confirm current branch is `rerun-corrected-mvlmm-epistasis`.
   - Confirm `main` was not modified.
   - Confirm rebuilt outputs live under `farm_outputs/original_logic_rebuild`, not `manuscript/`.

2. **Code change audit**
   - Review every file changed by commits `6ac2d87`, `5f60c79`, `03f756e`.
   - Identify any unintended scientific logic changes beyond adding the 13 validated markers and rerunning epistasis on the corrected 4,052-candidate universe.

3. **Original logic replication**
   - Verify old 354-row public table counts from `main`.
   - Verify old 88/81 multi-method calculation uses component-expanded frame.
   - Verify old epistasis candidate universe is 3,542 and corrected original-rule universe is 4,052.

4. **P-value scale audit**
   - Confirm additive evidence uses Galwey-adjusted p-values, not raw p-values.
   - Confirm uvLMM display counts use adjusted p-values.
   - Confirm no column named `pv20_adj_galwey` contains raw `pv20`.

5. **Epistasis support audit**
   - Confirm support requires both historical adjusted p threshold and lower-bound effect filter.
   - Confirm p-value-only counts are not used as manuscript evidence.
   - Confirm parent-marker support pivoting is correct.

6. **Permutation audit**
   - Confirm 100 permutations and 405,200 permutation rows.
   - Confirm 25 non-ok rows are exactly the documented conservative p=1 fallback for permutation 9 array 146.
   - Confirm those rows cannot make manuscript support less conservative or incorrectly add significance.
   - Confirm recomputed raw/adjusted permutation thresholds are not silently used as manuscript thresholds.

7. **Headline calculation audit**
   - Independently recalculate all old-vs-new headline values listed above.
   - Confirm component-expanded and marker-level frames are never mixed.
   - Confirm per-PBP splits and gained/lost substitution lists.

8. **Head-node safety audit**
   - Confirm README no longer instructs users to run heavy scripts directly on head nodes.
   - Confirm LSF scripts cover smoke, observed, permutations, merge, review outputs, and missing permutation chunks.
   - Identify any remaining commands likely to trigger Arbiter penalties.

9. **Stale output audit**
   - Identify any old invalid manuscript-facing files still present or likely to be confused with final outputs.
   - Check whether `farm_outputs/original_logic_rebuild/tidied_invalid_manuscript_outputs` is clearly ignored/non-final.

10. **Manuscript implication audit**
    - State which manuscript claims would change:
      - any-evidence substitutions/positions/proportion
      - multi-method substitutions/positions/per-PBP split
      - additive+epistatic overlap
      - effect-size range
      - Supplementary File 1 row/evidence counts
      - Table 3 top rows
    - State which claims appear stable:
      - PBP2X I371T remains sole Very Strong row
      - Strong count remains 5
      - broad conclusion of substantial epistasis remains true

## Expected Deliverable

Write a markdown audit report. Do not change files unless explicitly instructed. The report should include:

- Critical findings, if any.
- Major concerns or uncertainties.
- Minor issues/nits.
- Exact commands or scripts used for independent recalculation.
- Recomputed old-vs-new headline table.
- A clear recommendation:
  - safe to use rebuilt review outputs;
  - use with caveats;
  - or do not trust until specific fixes/reruns are done.
