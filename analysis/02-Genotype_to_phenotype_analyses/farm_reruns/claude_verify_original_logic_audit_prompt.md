# Claude Prompt: Verify Original Logic Replication Audit

```text
I need you to independently verify a very specific audit of my AMR cartography repo. Do not trust the previous assistant's conclusions. Work from files on disk.

Repository:
- /nfs/users/nfs_a/ab69/AMR-cartography

Important constraints:
- Work read-only.
- Do not modify code, data, manuscript files, outputs, branches, or git state.
- Do not delete, reset, clean, restore, install packages, submit jobs, or rerun heavy analyses.
- Use CSV parsers for counts, not grep-only checks.
- Be surgical and precise. The key question is whether the corrected rerun pipeline exactly replicates the original manuscript logic, except for including the corrected marker panel.

Primary files to inspect:
- analysis/02-Genotype_to_phenotype_analyses/30-31-multivariate-linear-mixed-models-heritability-epistatic-analyses/30-mvLMM-and-heritability.Rmd
- analysis/02-Genotype_to_phenotype_analyses/30-31-multivariate-linear-mixed-models-heritability-epistatic-analyses/31-mvLMM-heritability-and-epistatic-mvLMM.py
- analysis/02-Genotype_to_phenotype_analyses/32-ranking-substitutions/32-Table-of-subs-and-overlap.Rmd
- analysis/02-Genotype_to_phenotype_analyses/farm_reruns/*.py
- analysis/02-Genotype_to_phenotype_analyses/farm_reruns/original_logic_replication_audit.md
- AMRC-repo-files/_rerun_snapshots/mvLMM-pre-rerun-20260505-154523/pythonProject1/
- AMRC-repo-files/pythonProject1-additive-production-20260507-150112/
- farm_outputs/corrected_epistasis/
- manuscript/Supplementary_File_1.csv
- manuscript/Supplementary_File_1_corrected_marker_level.csv
- manuscript/corrected_rerun_manuscript_audit.md

Your tasks:

1. Establish repo state
- Run `pwd`.
- Run `git status --short --branch`.
- Run `git diff --stat`.
- Identify whether the audit files are modified/untracked.

2. Verify the 13 added markers
Using Python/pandas:
- Load old markers from:
  AMRC-repo-files/_rerun_snapshots/mvLMM-pre-rerun-20260505-154523/pythonProject1/S.pneumo_map_dummy_gen_test_markers.csv
- Load corrected markers from:
  AMRC-repo-files/pythonProject1-additive-production-20260507-150112/S.pneumo_map_dummy_gen_test_markers.csv
- Load added markers from:
  farm_outputs/corrected_epistasis/additive/added_markers.csv
- Verify:
  - old marker count = 157
  - corrected marker count = 170
  - old set is a strict subset of corrected set
  - corrected minus old exactly equals the 13 added markers
  - all 13 have present and absent counts greater than 1 percent of 3,620 isolates
  - all 13 are binary, have no missing values, and have no exact duplicate or exact complement in the corrected 170-marker matrix
- Report the present/absent counts for each added marker.

3. Verify original additive marker inclusion logic
Inspect `30-mvLMM-and-heritability.Rmd` and confirm:
- after exact duplicate/inverse merging, `S.pneumo_map_dummy_gen_relatedness_matrix.csv` is written
- additive test markers are selected by `Present > one_percent & Absent > one_percent`
- `one_percent = nrow(dummy_gen) * 0.01`
- `S.pneumo_map_dummy_gen_test_markers.csv` is then written

4. Verify additive downstream threshold scale
Inspect the original Rmds and old files. Confirm or refute:
- Original additive raw p-values are in `mvLMM_p_values_normal_pneumo_low_freq_vars.csv`.
- Original downstream table used `Sub_effect_sizes_mv_pneumo.csv`.
- `Sub_effect_sizes_mv_pneumo.csv` contains `pv20_adj_galwey`.
- In the old file, `pv20_adj_galwey / pv20` is 28.
- `32-Table-of-subs-and-overlap.Rmd` defines additive evidence using `pv20_adj_galwey_mv_LMM <= 0.000588 & Joint_effsize_mv_LMM >= 1`.
- Current `build_corrected_evidence_table.py` reads raw `pv20`, labels it as `pv20_adj_galwey_mv_LMM`, and applies `pv20 < 0.0009078488974311251`.
- Decide whether this is exact original logic or a logic change.
- Quantify old differences:
  - count old adjusted `pv20_adj_galwey <= 0.000588`
  - count old raw `pv20 <= 0.000588`
  - count old adjusted plus joint `>=1`
  - count old raw plus joint `>=1`

5. Verify uvLMM display threshold scale
Inspect `30-mvLMM-and-heritability.Rmd` around the uvLMM section and `merge_exact_unilmm.py`.
Confirm or refute:
- Original uvLMM comparison uses `pv20_adj_galwey = pv20 * galwey_meff`.
- Original uvLMM significance uses `pv20_adj_galwey <= 0.001`.
- uvLMM is not part of the four-stream evidence total, but is part of `Sig. mvLMM/uvLMM (No. drugs)` and manuscript comparison claims.
- Current corrected scripts count raw uvLMM `pv20 < 0.001`.
- Quantify old difference between raw and adjusted uvLMM counts.

6. Verify original epistasis candidate-generation logic
Inspect `30-mvLMM-and-heritability.Rmd` around interaction generation.
Confirm or refute:
- Original epistasis interaction candidates are generated from pairwise products.
- Original candidate filter requires all four genotype combinations to be at least 1 percent of isolates:
  `AA_Combination_0_0 >= one_percent`,
  `AA_Combination_0_1 >= one_percent`,
  `AA_Combination_1_0 >= one_percent`,
  `AA_Combination_1_1 >= one_percent`.
- Original output `S.pneumo_map_test_markers_incl_2nd_order_epistatic.csv` has 3,542 interaction columns.
- Old fitted `mvLMM_p_values_epi_pneumo.csv` has 3,542 rows.

7. Compare corrected epistasis candidate logic
Inspect `generate_corrected_epistasis_interactions.py`.
Confirm or refute:
- It starts from 170 choose 2 marker pairs.
- It removes constant interactions.
- It removes interactions equal to either parent.
- It removes duplicate interaction vectors by hash.
- It does not apply the original four-cell genotype-combination filter.
- It therefore does not exactly replicate the original candidate-generation logic.

Using Python/pandas, recompute old-rule unordered candidate pairs:
- On the old panel: expected 3,542.
- On the corrected 170-marker panel: expected 4,052.
- Verify old 3,542 unordered pairs are a subset of corrected 4,052.
- Verify the 510 new pairs all involve one of the 13 added markers.
- Compare current generated `farm_outputs/corrected_epistasis/interactions/corrected_epistasis_interactions.csv` to the corrected old-rule 4,052:
  - current count expected 4,563
  - overlap expected 2,412
  - current extra invalid under old rule expected 2,151
  - old-rule missing from current expected 1,640

8. Verify epistasis support logic
Inspect `30-mvLMM-and-heritability.Rmd` lines that build `Sig_AA_subs_all_epistatic_interactions_S.pneumo.csv`.
Confirm or refute:
- `No_sig_V1` is filtered with `pv20_adj_galwey <= 0.0007620121`.
- `No_sig_V1` is also filtered with `joint_effect_size - joint_effect_size_se >= 1`.
- The support file is built from `No_sig_V1`, not from all p-significant rows.
- The final evidence table then uses `PBP_count_mv_LMM_epistatic >= 1`.

Inspect `merge_epistasis_chunks.py` and confirm or refute:
- Current corrected support counts p-significant interactions only.
- It does not apply `joint_effect_size - joint_effect_size_se >= 1`.
- It does not use Galwey-adjusted epistasis p-values.

Quantify on the current invalid 4,563 interaction universe:
- raw p-significant interactions
- raw p-significant and `joint_effect_size >= 1`
- raw p-significant and `joint_effect_size - joint_effect_size_se >= 1`
- markers supported by raw p-only
- markers supported by raw p plus old lower-bound effect rule

9. Verify final evidence table logic
Inspect `32-Table-of-subs-and-overlap.Rmd` and `build_corrected_evidence_table.py`.
Confirm or refute:
- original component-expanded frame uses `separate_rows(Location, sep="/")`
- original evidence streams are single-substitution, clustering, additive mvLMM, and epistasis
- original total-to-label mapping is 0 Weak/No Evidence, 1 Weak, 2 Moderate, 3 Strong, 4 Very Strong
- current builder preserves the component-expanded frame and category mapping
- current builder's upstream additive/uv/epistasis inputs are not exact original logic

10. Final answer format
Write a concise but precise report with:
- "Confirmed" items
- "Refuted" items
- "Logic deviations"
- "What must be fixed before manuscript regeneration"
- Whether `manuscript/Supplementary_File_1.csv`, `Supplementary_File_1_corrected_marker_level.csv`, and `corrected_rerun_manuscript_audit.md` should be considered final or invalid pending rerun

Be very careful with language:
- Say "the 13-marker additive-panel correction is needed" only if verified.
- Say "the current corrected downstream counts are not final" if the logic deviations are verified.
- Do not propose changing scientific methods silently. Distinguish exact original logic from potentially defensible new statistical logic.
```

