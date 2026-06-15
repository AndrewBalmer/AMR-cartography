# Claude Independent Verification And Action Plan

Date captured: 2026-06-15

Source: user-provided Claude verification text.

## Verbatim Claude Report

```text
ndent verification complete — I worked read-only against the repo (the only write was to my own memory directory, outside the repo). I reproduced every numeric claim in original_logic_replication_audit.md exactly from files on disk. Below is my report.

Independent Verification Report — Original Logic Replication Audit
Date: 2026-06-15 · Branch: rerun-corrected-mvlmm-epistasis · Method: read-only, pandas/numpy via .venv-mvlmm, CSV parsers for all counts.

Repo state (Task 1)
pwd = /nfs/users/nfs_a/ab69/AMR-cartography; branch tracks origin/rerun-corrected-mvlmm-epistasis.
Modified (tracked): README.md, build_corrected_evidence_table.py, Supplementary_File_1.csv.
Untracked: .venv-mvlmm/, the two audit prompt/report .md files, Supplementary_File_1_corrected_marker_level.csv, corrected_rerun_manuscript_audit.md. The audit report itself is untracked (not committed).
Confirmed
The 13-marker additive-panel correction is needed and valid (Task 2). Old panel = 157, corrected = 170, both over 3,620 isolates; old is a strict subset of corrected; corrected − old is exactly the 13 markers in added_markers.csv. All 13 are binary, have zero missing values, and have no exact duplicate and no exact complement column anywhere in the 170-marker matrix. All 13 pass the original Present > 36.2 & Absent > 36.2 rule (min present = 82). Present/absent counts match the audit table row-for-row (e.g. PBP1A_S458_D_PBP1A_Y487_F 82/3538; PBP2B_D625_G_PBP2B_T630_N 710/2910).

Original additive marker-inclusion logic (Task 3). 30-mvLMM-and-heritability.Rmd:234 writes S.pneumo_map_dummy_gen_relatedness_matrix.csv after the correlation==±1 merge of identical/inverse variables (L197–230). one_percent = nrow(dummy_gen)*0.01 (L256); filter(Present > one_percent & Absent > one_percent) (L258–259); S.pneumo_map_dummy_gen_test_markers.csv written (L269). All confirmed.

Original downstream thresholds and rules.

Additive (Task 4): raw p in mvLMM_p_values_normal_pneumo_low_freq_vars.csv; downstream Sub_effect_sizes_mv_pneumo.csv carries pv20_adj_galwey; pv20_adj_galwey / pv20 = 28.0 exactly (single unique ratio). Evidence rule 32-Table-of-subs-and-overlap.Rmd:554: pv20_adj_galwey_mv_LMM <= 0.000588 & Joint_effsize_mv_LMM >= 1.
uvLMM (Task 5): pv20_adj_galwey = pv20 * galwey_meff then <= 0.001 (30…Rmd:1218,1221; count at L1292). uvLMM is not one of the four total streams (L561) but feeds the Sig. mvLMM/uvLMM (No. drugs) column and comparison claims.
Epistasis candidates (Task 6): pairwise products via model.matrix(Y ~ .^2); four-cell filter AA_Combination_{00,01,10,11} >= one_percent (L344–348); output and old fitted file both 3,542.
Epistasis support (Task 8): galwey_meff <- 39, significance pv20_adj_galwey <= 0.0007620121 (L799–804); No_sig_V1 requires both that threshold and joint_effect_size - joint_effect_size_se >= 1 (L905–907); the support file is built from No_sig_V1 (L937–944); evidence rule uses PBP_count_mv_LMM_epistatic >= 1 (32…Rmd:555).
Evidence table (Task 9): component expansion separate_rows(Location, sep="/") (32…Rmd:504); four streams (single-sub median≥1, clustering Confidence=="Strong", additive, epistasis); total→label = 0 Weak/No Evidence, 1 Weak, 2 Moderate, 3 Strong, 4 Very Strong (L599–624).
Quantified old-scale differences — all reproduced exactly:

Check	Adjusted (original)	Raw (as corrected scores it)
Additive ≤0.000588	118	128
Additive + Joint≥1	78	84
uvLMM ≤0.001 (tests / markers)	66 / 40	127 / 61
Epistasis candidate universe — all reproduced exactly: old-rule pairs = 3,542 (old panel) and 4,052 (corrected panel); old ⊂ corrected; 510 new pairs all involve an added marker. Independently, the old candidate file's parsed pairs equal the recomputed old-rule pairs — confirming the four-cell logic produced the historical file. Current generator = 4,563; overlap with old-rule-4,052 = 2,412; current-only (invalid under old rule) = 2,151; old-rule-missing-from-current = 1,640 (2412+2151=4563; 2412+1640=4052).

Epistasis support on the current 4,563 universe — all reproduced exactly: raw p-significant = 2,721; +joint≥1 = 2,363; + (joint−joint_se)≥1 = 1,870; markers by p-only = 155; markers by p + lower-bound effect = 150.

Refuted
Nothing in the audit was refuted. What is refuted is the proposition that the corrected pipeline exactly replicates the original manuscript logic — it does not. I found no error in the audit's facts, counts, or source-line citations.

Logic deviations (corrected pipeline vs original implemented logic)
Additive p-value scale & threshold. Original scores Galwey-adjusted pv20×28 at ≤0.000588; build_corrected_evidence_table.py:244,357 puts raw pv20 into the column named pv20_adj_galwey_mv_LMM and thresholds at 0.0009078 (common.py:33). On the raw scale the original gate is ≈ pv20 ≤ 2.1e-5; the corrected gate is ≤ 9.08e-4 — roughly 43× more permissive. This is a method change, not exact replication. (It may be a defensible permutation-based threshold, but it is a different scale and derivation, applied through a deceptively original-looking column name.)
uvLMM scale. merge_exact_unilmm.py and the builder count raw pv20 < 0.001; original used Galwey-adjusted ≤ 0.001.
Epistasis candidate generation. generate_corrected_epistasis_interactions.py starts from 170-choose-2, drops constant / parent-identical / hash-duplicate vectors, and omits the four-cell ≥1% filter. Result is a different universe (4,563 vs the correct-old-rule 4,052), with 2,151 pairs the original would exclude and 1,640 it would include.
Epistasis support. merge_epistasis_chunks.py:60–81 marks support by raw-pv20 permutation significance only — no Galwey adjustment and no joint_effect_size − joint_effect_size_se ≥ 1 lower-bound effect filter that the original No_sig_V1 required.
(The builder does faithfully preserve the component-expanded frame, the four-stream structure, and the 0–4 category mapping — the deviations are entirely in the upstream additive/uv/epistasis inputs.)

What must be fixed before manuscript regeneration
Golden old-data test first. Run the rewritten pipeline on the old 157-marker data and require it to reproduce: 157 additive markers, 3,542 epistasis candidates, 354-row Supplementary_File_1.csv, and the old evidence counts (1 Very Strong / 5 Strong / 82 Moderate / 105 Weak / 161 None; 88 multi-method). If it fails, do not use corrected outputs.
Keep the 170-marker panel (the 13-marker correction is verified valid).
Additive: compute corrected pv20_adj_galwey on the 170-marker panel and threshold on that adjusted scale — or explicitly document a deliberate switch to a raw permutation threshold as a method change.
uvLMM: count drugs on Galwey-adjusted ≤ 0.001.
Epistasis candidates: regenerate with the four-cell ≥1% filter and without hash-dedup (target 4,052; old 3,542 must be a subset).
Epistasis support: rebuild on that 4,052 universe with corrected-scale p-values and the joint − joint_se ≥ 1 effect lower bound before pivoting to PBP_count.
Rebuild the supplementary and marker-level tables and recompute all manuscript summary stats from the rebuilt frame.
Validity verdict
manuscript/Supplementary_File_1.csv, manuscript/Supplementary_File_1_corrected_marker_level.csv, and manuscript/corrected_rerun_manuscript_audit.md should be treated as INVALID pending rerun. The current corrected downstream counts are not final: the four logic deviations above are verified, so any multi-method / evidence-tier numbers produced by this pipeline (e.g. the 155/150 support counts, and the summary.json 2721 significant interactions over the wrong 4,563 universe) rest on logic that does not match the original manuscript. The 13-marker additive-panel correction is needed; the rest of the corrected pipeline must be brought back onto the original p-value scales and filters (or its method changes explicitly documented) before regeneration.
```

## Action Plan For Claude's Seven Fixes

### 1. Golden Old-Data Test First

Goal: build a reproducibility harness that runs the corrected table-building logic against the old 157-marker snapshot and refuses to proceed unless it reproduces the original manuscript outputs.

Tasks:

1. Create a small verification script that reads the old snapshot under `AMRC-repo-files/_rerun_snapshots/mvLMM-pre-rerun-20260505-154523/`.
2. Assert the old marker count is `157`.
3. Recompute old epistasis candidates and assert `3542`.
4. Rebuild the old component-expanded supplementary frame and assert:
   - total rows `354`
   - evidence counts `Very Strong=1`, `Strong=5`, `Moderate=82`, `Weak=105`, `Weak/No Evidence=161`
   - multi-method rows `88`
5. Save this as a fast local test/report before touching corrected outputs.

Success criterion: old-data reconstruction exactly matches the original `manuscript/Supplementary_File_1.csv` logic and counts.

### 2. Keep The 170-Marker Panel

Goal: preserve the corrected additive marker panel while proving the 13 added markers pass original inclusion rules.

Tasks:

1. Add a marker-panel validation step that compares old `157` and corrected `170`.
2. Assert old markers are a strict subset of corrected markers.
3. Assert corrected-minus-old equals `added_markers.csv`.
4. Assert all 13 added markers pass binary, no-NA, no-duplicate, no-complement, and `Present > 1%`, `Absent > 1%` checks.
5. Emit a small `added_marker_validation.csv` or markdown table.

Success criterion: the 170-marker panel is locked as the corrected panel and cannot silently drift.

### 3. Additive: Restore Original P-Value Scale Or Explicitly Declare A Method Change

Goal: separate exact original-logic replication from any possible new statistical threshold.

Tasks:

1. Compute corrected additive `pv20_adj_galwey` on the 170-marker panel.
2. Confirm the correct Galwey/meff factor for the corrected 170-marker matrix.
3. Store both raw `pv20` and adjusted `pv20_adj_galwey` with unambiguous column names.
4. Update evidence scoring to use adjusted p-values if the goal is exact original logic.
5. Add a switch or config flag only if a deliberate raw-permutation-threshold method change is desired.
6. Re-run the golden old-data test to prove old additive evidence is reproduced before applying to 170 markers.

Success criterion: additive evidence no longer places raw p-values into a column named `pv20_adj_galwey_mv_LMM`, and the evidence rule is explicit and test-covered.

### 4. uvLMM: Count Drugs On Galwey-Adjusted P-Values

Goal: make `Sig. mvLMM/uvLMM (No. drugs)` match the original display/comparison logic.

Tasks:

1. Compute uvLMM `pv20_adj_galwey` using the same marker-matrix correction factor logic as the original.
2. Count significant drugs with adjusted `<= 0.001`.
3. Preserve raw exact uvLMM values in a separate column/file for transparency.
4. Rebuild the uvLMM support table and compare old-data counts to the original `comparison_all_output.csv` behavior.

Success criterion: old-data uvLMM support reproduces adjusted old counts, and corrected 170-marker uvLMM support uses the same scale.

### 5. Epistasis Candidates: Four-Cell Filter, No Hash Dedup

Goal: regenerate the corrected epistasis candidate universe using the original candidate-selection rule.

Tasks:

1. Replace or add a strict old-logic candidate generator.
2. Generate all pairwise marker products from the corrected marker panel.
3. Keep only interactions where all four genotype cells are at least `1%` of isolates.
4. Do not hash-deduplicate interactions if the goal is exact original logic.
5. Assert:
   - old-data candidate count `3542`
   - corrected-panel candidate count `4052`
   - old unordered pairs are a subset of corrected unordered pairs
   - `510` corrected-only pairs all involve added markers
6. Write the candidate metadata and matrix with clear provenance.

Success criterion: corrected epistasis candidate universe is `4052`, not `4563`.

### 6. Epistasis Support: Correct Universe, Correct Scale, Lower-Bound Effect Filter

Goal: rebuild epistasis support exactly as the manuscript did.

Tasks:

1. Rerun observed epistasis on the 4,052 old-logic corrected candidate universe.
2. Rerun permutations on the same 4,052 universe.
3. Compute corrected p-values on the same scale as the original epistasis workflow.
4. Apply the old support filter:
   - threshold-significant
   - `joint_effect_size - joint_effect_size_se >= 1`
5. Pivot each supported interaction to both parent markers.
6. Count supported interactions per marker as `PBP_count`.
7. Use `PBP_count >= 1` as the epistasis evidence stream.
8. Re-run golden old-data test for epistasis support before trusting corrected 170-marker support.

Success criterion: epistasis evidence is no longer based on p-only support from the wrong 4,563-interaction universe.

### 7. Rebuild Tables And Recompute Manuscript Statistics

Goal: regenerate manuscript-facing outputs only after the upstream logic matches the original.

Tasks:

1. Rebuild `Supplementary_File_1.csv` in the original component-expanded frame.
2. Rebuild a separate marker-level audit table, clearly labeled as marker-level and not comparable to old `88/81`.
3. Recompute all manuscript summary statistics from the rebuilt component-expanded table:
   - any-evidence substitutions and positions
   - multi-method rows and positions
   - per-PBP splits
   - additive plus epistatic overlap
   - mvLMM-only vs uvLMM comparison claims
   - effect-size ranges
   - Table 3 top rows
4. Add assertions for expected old-data counts and sanity checks for corrected data.
5. Update the manuscript only after the rebuilt outputs pass the old-data test and corrected-panel validation.

Success criterion: every manuscript number is traceable to a rebuilt output generated by exact original logic plus the validated 170-marker panel.

