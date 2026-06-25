# Claude Prompt: Full Recomputed 170-Marker Pipeline And Manuscript/Supplement Cross-Reference Audit

You are auditing the AMR-cartography repository on branch `rerun-corrected-mvlmm-epistasis`.

Work read-only. Do not edit code, manuscript files, generated outputs, or scheduler state. Do not submit LSF jobs and do not rerun heavy analyses. Lightweight local inspection of CSV/JSON/PDF-text output is allowed.

The goal is to independently verify this claim:

> The primary revised output in `farm_outputs/recomputed_170_thresholds/` fully reran the original analysis logic on the corrected 170-marker panel, including additive mvLMM, uvLMM, epistasis, permutations, thresholds, downstream merging, Supplementary File 1, marker-level audit table, and manuscript statistics. The manuscript and old supplement have not yet been overwritten.

Then cross-reference the recomputed results against the existing rendered manuscript and supplementary information PDFs, and identify every old result statement that must change. For each old statement, decide whether the original wording remains scientifically valid after updating the number, or whether the wording/framing itself must change.

## Critical Frame Rules

- The existing manuscript source is not in the repository; use PDF text extraction from:
  - `manuscript/Manuscript.pdf`
  - `manuscript/Supplementary_Information.pdf`
- Treat PDF text line numbers as extraction aids, not authoritative source lines.
- The manuscript's original `88 changes at 81 sites` count used the public component-expanded Supplementary File 1 frame.
- Do not mix component-expanded public counts with marker-level counts.
- Use `component-expanded substitution rows` or `component-expanded substitutions` when referring to the public Supplementary File 1 frame.
- Use `markers` only for the underlying 170 binary model predictors.
- `multi-method` in the manuscript frame means `Evidence` is `Very Strong`, `Strong`, or `Moderate`.
- Evidence tiers remain:
  - `4 = Very Strong`
  - `3 = Strong`
  - `2 = Moderate`
  - `1 = Weak`
  - `0 = Weak/No Evidence`
- uvLMM is a display/comparison column, not a fifth evidence stream.
- For final revised analysis, thresholds are recomputed on the corrected 170-marker / 4,052-interaction universe using the original lowest-minimum-p logic. Do not use historical literal thresholds except when explicitly discussing the comparability output under `farm_outputs/original_logic_rebuild/`.

## Source Files To Audit

Primary recomputed output root:

- `farm_outputs/recomputed_170_thresholds/`

Validation/provenance:

- `farm_outputs/recomputed_170_thresholds/validation/recomputed_170_validation_report.md`
- `farm_outputs/recomputed_170_thresholds/thresholds/recomputed_thresholds.json`
- `farm_outputs/recomputed_170_thresholds/thresholds/recomputed_meff.json`
- `farm_outputs/recomputed_170_thresholds/thresholds/additive_threshold_summary.json`
- `farm_outputs/recomputed_170_thresholds/thresholds/uvlmm_summary.json`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_summary.json`

Primary recomputed tables:

- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/Supplementary_File_1.csv`
- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/recomputed_170_manuscript_statistics.md`
- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/surgical_result_audit/recomputed_170_surgical_result_audit.md`
- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/surgical_result_audit/recomputed_170_table3_top20.csv`
- `farm_outputs/recomputed_170_thresholds/manuscript_outputs/surgical_result_audit/recomputed_170_manuscript_replacement_text.md`

Old/preprint and historical-comparability tables:

- `manuscript/Supplementary_File_1.csv`
- `git show main:manuscript/Supplementary_File_1.csv`
- `farm_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv`

Model output tables:

- `farm_outputs/recomputed_170_thresholds/additive/merged/mvLMM_p_values_normal_pneumo_low_freq_vars_adjusted.csv`
- `farm_outputs/recomputed_170_thresholds/additive/merged/mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv`
- `farm_outputs/recomputed_170_thresholds/additive/merged/mvLMM_p_values_normal_pneumo_random_phenotype_FWAS_adjusted.csv`
- `farm_outputs/recomputed_170_thresholds/uvlmm/merged/uniLMM_exact_170_marker_drug_results.csv`
- `farm_outputs/recomputed_170_thresholds/uvlmm/merged/uniLMM_exact_170_marker_support.csv`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_p_values.csv`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_permutation_p_values.csv`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_permutation_minima.csv`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_supported_interactions.csv`
- `farm_outputs/recomputed_170_thresholds/epistasis/merged/corrected_epistasis_marker_support.csv`

Pipeline scripts to inspect for logic/provenance:

- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/common.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/run_additive_chunk.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_additive_chunks.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/run_exact_unilmm_chunk.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_exact_unilmm_chunks.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/generate_corrected_epistasis_interactions.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/run_epistasis_chunk.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/run_epistasis_permutation_chunk.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/merge_epistasis_chunks.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/build_corrected_evidence_table.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/build_recomputed_summary_report.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/build_recomputed_surgical_result_audit.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/validate_original_logic.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/validate_recomputed_outputs.py`
- `analysis/02-Genotype_to_phenotype_analyses/farm_reruns/lsf/`

Original Rmd logic to compare against:

- `analysis/02-Genotype_to_phenotype_analyses/32-ranking-substitutions/32-Table-of-subs-and-overlap.Rmd`
- `analysis/02-Genotype_to_phenotype_analyses/30-31-multivariate-linear-mixed-models-heritability-epistatic-analyses/30-mvLMM-and-heritability.Rmd`

## Expected Full-Rerun Validation Values

Independently verify these from files, not by trusting the generated Markdown:

- Current branch is `rerun-corrected-mvlmm-epistasis`.
- No file under `manuscript/` is modified by the recomputed workflow.
- No active AMR-cartography additive/uvLMM/epistasis/recomputed jobs remain in LSF.
- Additive observed p-value rows: `170`.
- Additive observed effect rows: `680`.
- Additive permutation rows: `17,000`.
- Additive permutations: `100`.
- Additive non-ok rows: `0`.
- uvLMM rows: `1,020`.
- uvLMM markers: `170`.
- uvLMM non-ok rows: `0`.
- Epistasis candidates / observed p-value rows: `4,052`.
- Epistasis observed effect rows: `8,104`.
- Epistasis permutation rows: `405,200`.
- Epistasis permutations: `100`.
- Epistasis non-ok rows: `0`.
- Additive threshold policy: `lowest permutation minimum`.
- Epistasis threshold policy: `lowest permutation minimum`.
- Additive old calibration reproduced: expected old lowest raw minimum p `0.0005875153233810342`.
- Recomputed additive Galwey meff: `28`.
- Recomputed epistasis Galwey meff: `38`.
- Marker-level table rows: `170`.

Expected recomputed thresholds:

- Additive raw threshold: `0.0001418518635803`.
- Additive Galwey-adjusted threshold: `0.0039718521802484`.
- Epistasis raw lowest-min-p threshold: `4.4080441256643215e-05`.
- Epistasis Galwey-adjusted threshold: `0.0016750567677524422`.
- Epistasis p-threshold-only interactions: `3,147`.
- Epistasis supported interactions after lower-bound effect filter: `1,926`.
- Markers with epistasis support: `131`.

Expected headline Supplementary File 1 comparison, using the public component-expanded frame:

| Metric | Old preprint | Historical-threshold corrected | Recomputed 170-marker |
|---|---:|---:|---:|
| Public rows | 354 | 394 | 394 |
| Evidence VS/S/M/W/None | 1/5/82/105/161 | 1/5/109/121/158 | 1/5/111/119/158 |
| Any-evidence rows/subs/positions | 193/172/147 | 236/215/179 | 236/215/179 |
| Any-evidence positions 1A/2B/2X | 50/31/66 | 60/50/69 | 60/50/69 |
| Multi-method rows/subs/positions | 88/88/81 | 115/115/106 | 117/117/108 |
| Multi-method positions 1A/2B/2X | 18/24/39 | 24/42/40 | 25/42/41 |
| Multi-method beta range | 0.330-3.079 | 0.374-3.187 | 0.374-3.187 |

Expected recomputed marker-level evidence counts:

- Very Strong: `1`
- Strong: `5`
- Moderate: `80`
- Weak: `66`
- Weak/No Evidence: `18`
- Marker-level multi-evidence markers: `86` of `170`

Expected gained/lost multi-method rows:

- Recomputed vs old preprint:
  - Gained: `33`
  - Lost: `4`
  - Shared rows with changed displayed values: `84`
- Recomputed vs historical-threshold corrected:
  - Gained: `2`
  - Lost: `0`
  - Shared rows with changed displayed values: `6`

Expected recomputed Table 3 top row:

- `PBP2X I371T`
- Evidence `Very Strong`
- beta Joint `3.187`
- adjusted p-value `<1e-16`
- significant interactions `40`
- total evidence `4`
- `Sig. mvLMM/uvLMM` = `Yes/Yes (6)`

Strong substitutions should remain:

- `PBP2X E320K`
- `PBP2X A491V`
- `PBP2X T490S`
- `PBP2X A369V`
- `PBP2X M343T`

## Known Manuscript PDF Result Passages To Cross-Reference

Use `pdftotext -layout` or equivalent to inspect all result-bearing text. At minimum, audit these extracted locations from `manuscript/Manuscript.pdf`:

- Lines ~505-512: four analysis streams and "multiple methods" framing.
- Lines ~515-522:
  - old `172 amino acid changes`
  - old `147 unique PBP locations`
  - old `51.6%`
  - old `88 changes at 81 sites`
  - old PBP split `18 in PBP1A, 24 in PBP2B, 39 in PBP2X`
  - old effect range `0.33-3.079`
  - old additive/epistatic percentage `89.8%`
- Lines ~525-530:
  - old subclass-specific percentage `59.08%`
  - old common-effects percentage `13.63%`
  - examples `I371T`, `N444S`, `S531Y`, `P568S`, `L609S`, `G467L/S469N/G483A/A490S`
- Lines ~549-553:
  - old claim `Multivariate mvLMMs detected 94 substitutions ... not identified by univariate LMMs`
- Lines ~574-625:
  - Table 3 caption and top-20 table rows.
- Lines ~657-659:
  - old Discussion claim `identify 88 PBP substitutions`
- Lines ~700-704:
  - old effect range `0.33-3.079`
  - old `147 of 285 variable positions (51.6%)`
  - old `81 ... positions had evidence from multiple methods`
- Lines ~729-734:
  - old subclass-specific percentage `59.08%`
  - example wording.

For each, report:

- old exact wording/value,
- recomputed replacement value,
- whether the original wording remains applicable,
- recommended revised wording,
- source file used to derive the replacement value,
- any unresolved ambiguity.

## Known Supplementary Information PDF Result Passages To Cross-Reference

At minimum, audit these extracted locations from `manuscript/Supplementary_Information.pdf`:

- Lines ~1708-1721, Supplementary Section 2.11 additive mvLMM:
  - old `157 independent tests`
  - old `97 ... associated with phenotypic change exceeding one log2 MIC unit ... p < 0.05`
  - old `100` random phenotype permutations
  - old lowest p-value `5.88 x 10-4`
  - old `89 of the 97` below threshold
  - old PBP split `25 PBP1A, 27 PBP2B, 37 PBP2X`
  - old `16 changes` with joint effect size over two MIC units
- Lines ~1724-1733, Supplementary Section 2.11 epistasis:
  - old `157 independent changes`
  - old `3542 interactions`
  - old `2129` significant interactions before lower-bound filter
  - old threshold `p < 7.62 x 10-4`
  - old `1634` interactions after lower-bound filter
  - old `138 of 285 variant PBP locations (48.4%)`
  - old `30 substitutions` with more than 40 interactions
  - old dense-connection split `1A (36.2%) and 2X (42%)`
- Lines ~1756-1768, Supplementary Section 2.12 evidence ranking:
  - old evidence definitions
  - old `172 amino acid substitutions`
  - old `147 locations (51.6%)`
  - old `88 of 172`
  - old `81 unique locations`
  - old PBP split `18, 24, 39`
  - old additive/epistatic percentage `89.8%`
- Lines ~1774-1826, Supplementary Table S14:
  - old top-20 table values and membership.
- Lines ~1845-1855, Supplementary Section 3.1 subclass/common/contrasting effects:
  - old `172 substitutions and 2129 epistatic interactions`
  - old `217 (7.76%)`
  - old `693 (24.78%)`
  - old `31 (1.11%)`
  - old `1360 (48.63%)`
  - old `296 (10.59%)`
  - examples and interpretive wording.
- Lines ~1875-1884, Supplementary Section 3.2 overlap with previous work:
  - old `77 positions`
  - old `11 positions`
  - old `22 positions`
  - old `20/22`
  - old `97 positions not highlighted by previous methods`
  - Determine which, if any, should change after the 170-marker rebuild. If the external previous-work comparison inputs are unavailable, say so explicitly.
- Lines ~1904-1917, Supplementary Section 3.3 mvLMM vs uvLMM:
  - old claim extracted as `82 substitutions were identified by either multivariate or univariate analysis`
  - parenthetical counts `36 both, 53 multivariate only, 4 univariate only`
  - note this is internally inconsistent as extracted because `36 + 53 + 4 = 93`; verify against source if possible.
  - old threshold wording `p < 0.001`
  - Recompute the comparable marker-level and/or component-expanded values from the regenerated uvLMM and mvLMM outputs, and state which frame each uses.

## Values Already Derived By Codex That Need Independent Verification

These are not authoritative until you reproduce them:

- Recomputed additive marker-level, `pv20_adj_galwey <= 0.05` and joint effect size `>= 1`: `95` markers.
- Recomputed additive marker-level, permutation threshold and joint effect size `>= 1`: `92` markers.
- Recomputed threshold-passing additive PBP split: `27 PBP1A`, `26 PBP2B`, `39 PBP2X`.
- Recomputed threshold-passing additive markers with joint effect size `> 2`: `18`.
- Recomputed multi-method public rows with both mvLMM display support and at least one supported epistatic interaction: `108 / 117 = 92.3%`.
- Recomputed marker-level multi-evidence markers with both additive mvLMM evidence and epistasis evidence: `73 / 86 = 84.9%`.
- The old `89.8%` additive-plus-epistasis percentage has not been pinned to a unique public-frame recipe. Do not quote a revised percentage unless you can reproduce the old 89.8% recipe first.
- Recomputed marker-level mvLMM vs uvLMM comparison using strict display thresholds:
  - mvLMM evidence markers: `92`
  - uvLMM-support markers: `45`
  - either mvLMM or uvLMM: `104`
  - both: `33`
  - mvLMM only: `59`
  - uvLMM only: `12`
- Recomputed public component-expanded `Sig. mvLMM/uvLMM` display counts:
  - either: `208`
  - both: `53`
  - mvLMM only: `153`
  - uvLMM only: `2`
  - This may not be the same frame as the old manuscript's `94 substitutions not identified by uvLMM`; resolve the old definition before proposing replacement prose.
- Recomputed marker-level additive subclass pattern among mvLMM evidence markers:
  - one-axis/subclass-specific effects: `62 / 92 = 67.39%`
  - same-direction common effects across axes: `9 / 92 = 9.78%`
  - opposite-sign two-axis effects: `9 / 92 = 9.78%`
  - Verify whether this is the correct denominator for replacing the old `59.08% / 13.63%` sentence. If not, derive the correct denominator from the original Rmd logic.
- Recomputed combined additive-p-threshold plus epistasis-p-threshold effect-size universe, using marker-level additive p-threshold rows plus epistasis p-threshold rows:
  - denominator `3281`
  - one-axis/subclass-specific effects `1773 / 3281 = 54.04%`
  - same-direction common effects `112 / 3281 = 3.41%`
  - positive penicillin-axis effects `262 / 3281 = 7.99%`
  - positive cephalosporin-axis effects `866 / 3281 = 26.39%`
  - positive on both axes `35 / 3281 = 1.07%`
  - negative penicillin-axis effects `1603 / 3281 = 48.86%`
  - negative cephalosporin-axis effects `336 / 3281 = 10.24%`
  - This may not exactly match the old Rmd's component-expanded denominator; verify before using.
- Recomputed epistasis supported-interaction parent markers:
  - supported parent markers: `131`
  - marker-level supported positions: `125`
  - component-expanded supported positions: `159`
  - markers with more than 40 supported interactions: `39`
  - by PBP among those `39`: PBP1A `14`, PBP2B `6`, PBP2X `19`
  - Verify which position frame matches the old `138 of 285 variant PBP locations` wording.

## Audit Tasks

1. Confirm branch, git status, and that `manuscript/` is untouched.
2. Confirm all recomputed validation checks and shape/completeness counts.
3. Confirm no old literal thresholds are used for the primary recomputed outputs.
4. Confirm thresholds are recomputed using lowest permutation minimum, not 5th percentile and not rounded `0.001`.
5. Confirm Galwey meff values come from the recomputed 170-marker and 4,052-interaction matrices.
6. Confirm additive observed and null outputs were produced under `farm_outputs/recomputed_170_thresholds/`, not copied from external paths.
7. Confirm uvLMM was rerun over all `170 x 6` marker-drug tests.
8. Confirm epistasis candidates were regenerated with the original four-cell rule and total `4,052`.
9. Confirm epistasis observed and 100 permutation outputs are complete, with no synthetic fallback chunks.
10. Confirm Supplementary File 1 rebuilt in the component-expanded public frame.
11. Confirm marker-level audit table rebuilt separately and not mixed into manuscript claims.
12. Independently recompute all headline counts from CSVs.
13. Independently recompute gained/lost multi-method rows and inspect the gained/lost lists for formatting or biological oddities.
14. Independently verify Table 3 top-20 generation logic and values.
15. Extract manuscript PDF text and supplementary PDF text, then audit all result-bearing statements listed above.
16. Search the full extracted PDFs for additional instances of old key values and result claims not listed above, including:
    - `157`, `170`, `172`, `147`, `51.6`, `88`, `81`, `89.8`, `59.08`, `13.63`, `94`, `97`, `89`, `3542`, `2129`, `1634`, `138`, `48.4`, `30`, `7.62`, `5.88`, `0.33`, `3.079`, `Table 3`, `S14`, `S15`, `S18`, `S19`.
17. For every changed passage, state whether the old wording remains applicable:
    - `number update only`
    - `wording still directionally true but needs numeric update`
    - `wording/frame must change`
    - `cannot determine from available repository files`
18. Recommend replacement wording for each manuscript/supplement passage that changes.
19. Explicitly flag values whose old recipe cannot be reproduced from the repository, especially:
    - `89.8%`
    - `59.08% / 13.63%`
    - main-text `94 substitutions not identified by uvLMM`
    - Supplementary Section 3.2 overlap with previous work
    - Supplementary Section 3.3 mvLMM vs uvLMM counts

## Desired Output

Write a Markdown audit report with these sections:

1. Overall verdict: `PASS`, `FAIL`, or `NEEDS CLARIFICATION`.
2. Full-rerun validation table: expected vs observed vs pass/fail.
3. Threshold/provenance audit.
4. Public Supplementary File 1 and marker-level count audit.
5. Table 3 / Supplementary Table S14 audit.
6. Manuscript cross-reference table:
   - PDF
   - extracted line(s)
   - old wording/value
   - recomputed value
   - wording applicability classification
   - recommended replacement text
   - source file(s)
7. Supplementary Information cross-reference table with the same columns.
8. Findings ordered by severity.
9. Residual ambiguities that Andrew must approve before manuscript revision.
10. A short list of exact manuscript/supplement passages that are safe to retain unchanged.

Be surgical and precise. Do not smooth over frame mismatches. If a number cannot be derived in the exact original frame, say so and explain what would be needed.
