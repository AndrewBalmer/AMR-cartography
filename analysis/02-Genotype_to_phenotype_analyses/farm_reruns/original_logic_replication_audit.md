# Original Logic Replication Audit

Date: 2026-06-15

Branch audited: `rerun-corrected-mvlmm-epistasis`

Purpose: determine whether the corrected rerun pipeline exactly replicates the original manuscript logic, except for including the corrected marker panel.

## Executive Verdict

The 13 added markers are valid and should not have been excluded from the additive marker panel.

However, the current corrected rerun pipeline does not exactly replicate the original manuscript logic. The current manuscript-facing corrected counts should not be used as final results until the deviations below are fixed and the pipeline is rerun or rebuilt.

The largest deviations are:

1. Additive mvLMM evidence is currently scored on raw `pv20`; the original table scored Galwey-adjusted `pv20_adj_galwey`.
2. uvLMM display/support counts are currently scored on raw exact `pv20`; the original comparison table used Galwey-adjusted uvLMM p-values.
3. Corrected epistasis interaction candidates do not use the original four-cell sample-size filter.
4. Corrected epistasis support currently uses p-value significance only; the original support file required significance and `joint_effect_size - joint_effect_size_se >= 1`.
5. Corrected epistasis currently removes duplicate interaction vectors by hash; the original Rmd did not remove these duplicates after the four-cell filter.

## Validity Of The 13 Added Markers

The marker-panel correction itself is needed.

Verified facts:

- Old fitted additive panel: 157 markers.
- Corrected fitted additive panel: 170 markers.
- Old 157 markers are a strict subset of the corrected 170.
- The 13 added markers exactly match `farm_outputs/corrected_epistasis/additive/added_markers.csv`.
- The 13 added markers are all present in the corrected `S.pneumo_map_dummy_gen_test_markers.csv`.
- All 13 are binary, have no missing values, have no exact duplicate marker column, and have no exact inverse/complement marker column in the corrected 170-marker test matrix.
- The original additive filter is `Present > one_percent & Absent > one_percent`, where `one_percent = nrow(dummy_gen) * 0.01`; for 3,620 isolates this means `Present > 36.2` and `Absent > 36.2`.
- The smallest present count among the 13 added markers is 82, so all 13 pass the original additive inclusion rule.
- The corrected relatedness matrix adds 26 columns relative to the old relatedness matrix; exactly 13 of those pass the original additive test-marker frequency filter. The other 13 have present counts from 1 to 22 and correctly remain excluded from the additive test panel.

Source logic:

- Additive count filter: `30-mvLMM-and-heritability.Rmd`, lines 255-266.
- Test marker output: `S.pneumo_map_dummy_gen_test_markers.csv`.

Added marker counts:

| Marker | Components | Present | Absent | Passes additive 1 percent rule |
|---|---:|---:|---:|---|
| `PBP1A_D473_S_PBP1A_L640_V` | 2 | 165 | 3455 | yes |
| `PBP1A_E397_V_PBP1A_L421_I` | 2 | 222 | 3398 | yes |
| `PBP1A_S458_D_PBP1A_Y487_F` | 2 | 82 | 3538 | yes |
| `PBP1A_T392_S_PBP1A_V408_L` | 2 | 219 | 3401 | yes |
| `PBP1A_T495_I_PBP1A_Y497_H_PBP1A_H503_N_PBP1A_N517_D` | 4 | 533 | 3087 | yes |
| `PBP2B_A533_G_PBP2B_V542_I_PBP2B_D555_E_PBP2B_K556_Q_PBP2B_L566_V_PBP2B_V574_I_PBP2B_M581_V_PBP2B_H585_Q_PBP2B_T595_G_PBP2B_G597_A` | 10 | 207 | 3413 | yes |
| `PBP2B_A619_G_PBP2B_S640_T_PBP2B_D641_E` | 3 | 188 | 3432 | yes |
| `PBP2B_D561_E_PBP2B_Q565_A_PBP2B_L566_I_PBP2B_S582_A` | 4 | 230 | 3390 | yes |
| `PBP2B_D625_G_PBP2B_T630_N` | 2 | 710 | 2910 | yes |
| `PBP2B_G597_P_PBP2B_N606_D_PBP2B_L609_T` | 3 | 234 | 3386 | yes |
| `PBP2B_Q565_S_PBP2B_Q567_E` | 2 | 206 | 3414 | yes |
| `PBP2B_S473_T_PBP2B_S480_A` | 2 | 474 | 3146 | yes |
| `PBP2X_V537_I_PBP2X_L546_V_PBP2X_S574_T` | 3 | 512 | 3108 | yes |

## Original Pipeline Logic

### 1. Marker Matrix Generation

Original implemented logic:

1. Build dummy amino-acid variables from PBP transpeptidase-domain positions.
2. Remove invariant variables.
3. Remove sensitive/reference amino acids by filtering cases where the amino acid equals the reference amino acid.
4. Merge/remove exactly identical and exactly inverse variables via correlation `1` or `-1`.
5. Save the merged post-correlation matrix as `S.pneumo_map_dummy_gen_relatedness_matrix.csv`.
6. Count present and absent isolates per merged marker.
7. Keep additive test markers only where present and absent counts are both greater than 1 percent of isolates.
8. Save those markers as `S.pneumo_map_dummy_gen_test_markers.csv`.

Replication status:

- Corrected 170-marker panel is consistent with this logic.
- Old 157 test markers are a strict subset of corrected 170 test markers.
- No evidence found that the 13 added markers violate an original exclusion rule.

### 2. Additive mvLMM Fit

Original model logic in `31-mvLMM-heritability-and-epistatic-mvLMM.py`:

- Candidate marker: one column from `S.pneumo_map_dummy_gen_test_markers.csv`.
- Phenotype: map dimensions `D1`, `D2`.
- Relatedness: `linear_kinship()` on the relatedness matrix with the tested marker dropped.
- Trait covariance: `A = map_coords.cov()`.
- Output raw p-values to `mvLMM_p_values_normal_pneumo_low_freq_vars.csv`.
- Output candidate effects to `mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv`.

Current corrected fitted outputs:

- Corrected output directory has 170 p-value rows and 680 effect-size rows.
- Shared 157 markers are unchanged at the marker-vector level.
- Shared-marker p-values are extremely concordant between old and corrected fitted outputs.

Important caveat:

- The committed `farm_reruns/` scripts do not contain the additive mvLMM runner that produced `pythonProject1-additive-production-20260507-150112`. The fitted output can be audited, but the additive production code path is not fully captured in `farm_reruns/`.

### 3. Additive mvLMM Threshold And Evidence

Original downstream table logic:

- The Rmd computes `pv20_adj_galwey = pv20 * galwey_meff`.
- In the old additive file, `pv20_adj_galwey / pv20 = 28`.
- The manuscript table reads `Sub_effect_sizes_mv_pneumo.csv`, which contains `pv20_adj_galwey`.
- Evidence is `pv20_adj_galwey_mv_LMM <= 0.000588 & Joint_effsize_mv_LMM >= 1`.

Source lines:

- Galwey adjustment: `30-mvLMM-and-heritability.Rmd`, lines 402-416.
- Final evidence rule: `32-Table-of-subs-and-overlap.Rmd`, line 554.

Current corrected logic:

- `build_corrected_evidence_table.py` reads raw `pv20` from `mvLMM_p_values_normal_pneumo_low_freq_vars.csv`.
- It labels this raw value as `pv20_adj_galwey_mv_LMM`.
- It applies `pv20 < 0.0009078488974311251`.
- It applies `joint_effect_size >= 1`.

Status: not an exact replication.

This may be a defensible corrected statistical threshold, but it is not the same implemented manuscript logic. Exact replication requires using the same p-value scale and threshold derivation/threshold application as the original pipeline.

Observed old-scale difference:

- Old adjusted-p rule `pv20_adj_galwey <= 0.000588`: 118 additive markers.
- Old raw-p rule `pv20 <= 0.000588`: 128 additive markers.
- Old adjusted-p plus joint effect `>=1`: 78 markers.
- Old raw-p plus joint effect `>=1`: 84 markers.

Those are not equivalent.

### 4. uvLMM Logic

Original downstream comparison logic:

- The Rmd computes `pv20_adj_galwey = pv20 * galwey_meff` for uvLMM.
- The uvLMM significance display/count uses `pv20_adj_galwey <= 0.001`.
- uvLMM does not contribute to the four-stream `Evidence` score.
- uvLMM does contribute to the display column `Sig. mvLMM/uvLMM (No. drugs)` and to manuscript comparison claims.

Source lines:

- `30-mvLMM-and-heritability.Rmd`, lines 1189-1304.
- Supplementary table join/display: `32-Table-of-subs-and-overlap.Rmd`, lines 715-809.

Current corrected logic:

- `merge_exact_unilmm.py` and `build_corrected_evidence_table.py` count raw exact uvLMM `pv20 < 0.001`.

Status: not an exact replication for uvLMM display/comparison claims.

Observed old-scale difference:

- Old raw uvLMM `pv20 < 0.001`: 127 marker-drug tests, 61 markers with at least one significant drug.
- Old Galwey-adjusted uvLMM `pv20_adj_galwey <= 0.001`: 66 marker-drug tests, 40 markers with at least one significant drug.

### 5. Epistasis Candidate Generation

Original implemented logic:

1. Build all pairwise interaction terms from the post-merge marker matrix.
2. Keep interaction columns only where all four genotype combinations have at least 1 percent of isolates:
   - `AA_Combination_0_0 >= one_percent`
   - `AA_Combination_0_1 >= one_percent`
   - `AA_Combination_1_0 >= one_percent`
   - `AA_Combination_1_1 >= one_percent`
3. Save the kept interaction matrix as `S.pneumo_map_test_markers_incl_2nd_order_epistatic.csv`.
4. Do not remove duplicate interaction vectors by hash after this filter.

Source lines:

- Candidate generation and four-cell filter: `30-mvLMM-and-heritability.Rmd`, lines 312-353.
- Old fitted epistasis output has exactly 3,542 p-value rows, matching the old 3,542-column candidate file.

Current corrected logic:

- `generate_corrected_epistasis_interactions.py` starts from 170 choose 2 additive marker pairs.
- It removes constant interaction terms.
- It removes interaction terms equal to either parent.
- It removes duplicate interaction vectors by hash.
- It does not apply the original four-cell genotype-combination filter.

Status: not an exact replication.

Quantified impact:

- Old-rule candidate generation on the old panel gives exactly 3,542 interactions, matching the old file.
- Old-rule candidate generation on the corrected 170-marker panel gives 4,052 interactions.
- The old 3,542 unordered parent pairs are a strict subset of those 4,052 corrected old-rule pairs.
- The 510 new old-rule pairs all involve at least one of the 13 added markers.
- Current corrected candidate generator produced 4,563 interactions.
- Overlap between current 4,563 and old-rule corrected 4,052 is only 2,412.
- Current generator includes 2,151 pairs the original four-cell logic would exclude.
- Current generator omits 1,640 pairs the original four-cell logic would include, mostly because of duplicate-vector removal and parent-identical logic.

### 6. Epistasis Model Fit

Original model logic:

- Candidate marker: one interaction column from `S.pneumo_map_test_markers_incl_2nd_order_epistatic.csv`.
- Parent covariates: the two parent marker columns.
- Relatedness: `linear_kinship()` on the post-merge relatedness matrix with the two parent columns dropped.
- Trait covariance/design: original code uses `A = mod_matrix` for observed epistasis.
- Output raw p-values and candidate effects.

Current corrected model:

- `run_epistasis_chunk.py` uses the parent markers as covariates and drops parent columns from the relatedness features.
- `common.py` uses a low-rank implementation intended to be algebraically equivalent to the historical model form.

Status: model form appears intended to match, but it cannot be considered a faithful corrected rerun while the candidate universe differs.

Required validation before final use:

- After regenerating the old-rule 4,052 interaction universe, run a small overlap comparison on old data or a fixed subset to confirm the low-rank implementation reproduces the historical Limix `scan()` outputs within numerical tolerance.

### 7. Epistasis Threshold And Support

Original implemented downstream logic:

- Epistasis p-values are Galwey-adjusted with `galwey_meff <- 39`.
- Significance label uses `pv20_adj_galwey <= 0.0007620121`.
- The support file used by the evidence table is not simply p-significant interactions.
- It is built from:
  - `pv20_adj_galwey <= 0.0007620121`
  - `joint_effect_size - joint_effect_size_se >= 1`
- It then pivots each significant interaction to both parent markers and counts interactions per marker as `PBP_count`.
- The evidence table uses `PBP_count_mv_LMM_epistatic >= 1`.

Source lines:

- Epistasis Galwey factor and threshold: `30-mvLMM-and-heritability.Rmd`, lines 794-803.
- Effect-filtered support file construction: `30-mvLMM-and-heritability.Rmd`, lines 905-944.
- Evidence table epistasis rule: `32-Table-of-subs-and-overlap.Rmd`, line 555.

Current corrected logic:

- `merge_epistasis_chunks.py` derives a raw p-value permutation threshold from raw `pv20`.
- It marks `significant_epistasis_perm = pv20 <= threshold`.
- It counts marker support from p-significant interactions only.
- It does not apply `joint_effect_size - joint_effect_size_se >= 1`.
- It does not use Galwey-adjusted epistasis p-values.

Status: not an exact replication.

Quantified on current invalid 4,563-interaction universe:

- Raw p-significant interactions: 2,721.
- Raw p-significant and `joint_effect_size >= 1`: 2,363.
- Raw p-significant and `joint_effect_size - joint_effect_size_se >= 1`: 1,870.
- Markers with p-significant support only: 155.
- Markers with p-significant plus old lower-bound effect support: 150.

These counts are not final because they are based on the wrong interaction universe, but they show that the support rule matters.

### 8. Evidence Table Reconstruction

Original evidence table logic:

- Full join single-substitution rows to additive mvLMM rows.
- Split compound additive `Location` rows with `separate_rows(Location, sep="/")`.
- Derive `Loci` from `Location`.
- Full join clustering by `PBP` and `Loci`.
- Fill missing `Location` from `Loci`.
- Split epistasis support rows with `separate_rows(Location, sep="/")`.
- Full join epistasis support by `PBP` and `Location`.
- Four evidence streams:
  - single substitution: `median_phenotypic_distance >= 1`
  - clustering: `Confidence == "Strong"`
  - additive mvLMM: adjusted p threshold and `Joint_effsize_mv_LMM >= 1`
  - epistasis: `PBP_count_mv_LMM_epistatic >= 1`, where `PBP_count` came from the already effect-filtered epistasis support file
- Total score:
  - 0: Weak/No Evidence
  - 1: Weak
  - 2: Moderate
  - 3: Strong
  - 4: Very Strong
- Supplementary file is sorted by descending `total`, then descending additive joint effect size.

Current corrected builder status:

- Component-expanded frame: mostly matches original intent.
- Evidence category mapping: matches.
- Single-substitution evidence: matches.
- Clustering evidence: likely matches, although clustering de-duplication should use explicit confidence ordering rather than alphabetical sort.
- Additive evidence: not exact because raw p-value scale is used.
- Epistasis evidence: not exact because upstream epistasis support is p-only and candidate universe is wrong.
- uvLMM display: not exact because raw p-value threshold is used.

## Current Outputs That Should Be Treated As Invalid For Manuscript Results

Do not use the following as final manuscript-facing corrected counts until the exact-replication fixes are made:

- `manuscript/Supplementary_File_1.csv`
- `manuscript/Supplementary_File_1_corrected_marker_level.csv`
- `manuscript/corrected_rerun_manuscript_audit.md`
- `farm_outputs/corrected_epistasis/merged/corrected_epistasis_summary.json`
- Current corrected epistasis marker-support counts
- Current multi-method evidence counts such as 141, 134, 121, or 101

Those values are products of a pipeline that does not replicate the original implemented logic.

## Required Exact-Replication Fixes

Before regenerating the manuscript or supplement, the corrected pipeline should be changed so that it has a "golden old-data test" and a corrected-panel run.

### Golden Old-Data Test

Run the rewritten Python/R pipeline on the old 157-marker data and require it to reproduce:

- Old additive test marker count: 157.
- Old epistasis candidate interaction count: 3,542.
- Old `Supplementary_File_1.csv` total rows: 354.
- Old evidence counts:
  - Very Strong: 1
  - Strong: 5
  - Moderate: 82
  - Weak: 105
  - Weak/No Evidence: 161
- Old multi-method rows: 88.
- Old multi-method unique positions: 81.
- Old any-evidence rows: 193.
- Old any-evidence unique substitutions: 172.
- Old any-evidence unique positions: 147.

If the old-data test fails, the corrected-panel outputs should not be used.

### Corrected Marker Panel

Keep the corrected 170-marker panel. The 13 added markers pass the original inclusion rules.

Expected corrected old-rule epistasis candidate count:

- 4,052 interactions.
- Old 3,542 unordered parent pairs must be a subset.
- The 510 new old-rule pairs must all involve at least one added marker.

### Additive Evidence

The corrected table must use the same p-value scale as the original table:

- Compute corrected `pv20_adj_galwey` on the corrected marker panel.
- Apply the corrected threshold on the same adjusted-p scale, or explicitly document that the thresholding method has changed.
- If the threshold is intentionally changed to a raw permutation threshold, that is not "same original logic" and should be treated as a scientific-method change.

### uvLMM Display

The corrected `Sig. mvLMM/uvLMM (No. drugs)` column must use the same uvLMM p-value scale as the original:

- Compute uvLMM `pv20_adj_galwey`.
- Count significant drugs using adjusted `<= 0.001`, unless a deliberate method change is documented.

### Epistasis Candidates

Replace the current candidate generator with the old candidate-generation logic:

- Generate all pairwise products.
- Keep only pairs where all four parent genotype combinations are at least 1 percent of isolates.
- Do not hash-deduplicate interaction vectors if the aim is exact original replication.
- Preserve a deterministic parent-pair key so old/new overlap can be audited as unordered pairs.

### Epistasis Support

Rebuild epistasis support exactly as before:

- Use the corrected old-rule interaction universe.
- Run observed epistasis and permutations on that universe.
- Compute corrected p-values on the same scale as the original.
- Mark support rows only where:
  - corrected epistasis p threshold is met
  - `joint_effect_size - joint_effect_size_se >= 1`
- Pivot significant interactions to both parent markers.
- Count interactions per parent marker as `PBP_count`.
- Evidence stream is `PBP_count >= 1`.

### Final Evidence Table

Only after upstream fixes:

- Rebuild the component-expanded supplementary table.
- Rebuild the marker-level audit table separately.
- Assert final shapes/counts.
- Recompute all manuscript summary statistics from the rebuilt component-expanded table.

## Minimal Independent Verification Commands

These checks were used for this audit and should be repeated by an independent reviewer.

```bash
cd /nfs/users/nfs_a/ab69/AMR-cartography
git status --short --branch
```

```bash
.venv-mvlmm/bin/python - <<'PY'
from pathlib import Path
import pandas as pd
root = Path('/nfs/users/nfs_a/ab69/AMR-cartography')
old = root/'AMRC-repo-files/_rerun_snapshots/mvLMM-pre-rerun-20260505-154523/pythonProject1'
new = root/'AMRC-repo-files/pythonProject1-additive-production-20260507-150112'
old_m = pd.read_csv(old/'S.pneumo_map_dummy_gen_test_markers.csv')
new_m = pd.read_csv(new/'S.pneumo_map_dummy_gen_test_markers.csv')
added = set(pd.read_csv(root/'farm_outputs/corrected_epistasis/additive/added_markers.csv')['marker'])
print(old_m.shape, new_m.shape)
print(set(old_m.columns).issubset(set(new_m.columns)))
print(set(new_m.columns) - set(old_m.columns) == added)
PY
```

```bash
.venv-mvlmm/bin/python - <<'PY'
from pathlib import Path
import pandas as pd, numpy as np
root = Path('/nfs/users/nfs_a/ab69/AMR-cartography')
old = root/'AMRC-repo-files/_rerun_snapshots/mvLMM-pre-rerun-20260505-154523/pythonProject1'
new = root/'AMRC-repo-files/pythonProject1-additive-production-20260507-150112'

def old_rule_pairs(path):
    cols = list(pd.read_csv(path/'S.pneumo_map_dummy_gen_test_markers.csv', nrows=0).columns)
    rel = pd.read_csv(path/'S.pneumo_map_dummy_gen_relatedness_matrix.csv', usecols=cols)
    x = rel.to_numpy(dtype=np.uint8)
    threshold = len(rel) * 0.01
    out = set()
    for i, a_name in enumerate(cols):
        a = x[:, i]
        for j in range(i + 1, len(cols)):
            b = x[:, j]
            n11 = int((a & b).sum())
            n10 = int((a & (1 - b)).sum())
            n01 = int(((1 - a) & b).sum())
            n00 = len(a) - n11 - n10 - n01
            if min(n00, n01, n10, n11) >= threshold:
                out.add(frozenset((a_name, cols[j])))
    return out

old_pairs = old_rule_pairs(old)
new_pairs = old_rule_pairs(new)
added = set(pd.read_csv(root/'farm_outputs/corrected_epistasis/additive/added_markers.csv')['marker'])
print(len(old_pairs), len(new_pairs), old_pairs.issubset(new_pairs), len(new_pairs - old_pairs))
print(all(bool(set(pair) & added) for pair in (new_pairs - old_pairs)))
PY
```

Expected output: old pairs `3542`, corrected old-rule pairs `4052`, old subset `True`, new-only `510`, all new-only pairs involve added markers `True`.

## Bottom Line

The correct next step is not manuscript regeneration.

The correct next step is to repair the corrected rerun pipeline so that it exactly reproduces the old implemented logic on old data, then apply that same logic to the corrected 170-marker panel.

