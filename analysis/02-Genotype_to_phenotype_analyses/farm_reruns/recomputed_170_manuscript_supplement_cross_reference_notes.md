# Recomputed 170-Marker Manuscript/Supplement Cross-Reference Notes

These notes are a Codex-derived starting point for independent audit. They should not be treated as final manuscript text until checked against the full PDF extraction and underlying CSV/JSON outputs.

## Current Primary Result Frame

- Primary output root: `farm_outputs/recomputed_170_thresholds/`
- Primary public table: `farm_outputs/recomputed_170_thresholds/manuscript_outputs/Supplementary_File_1.csv`
- Marker-level audit table: `farm_outputs/recomputed_170_thresholds/manuscript_outputs/Supplementary_File_1_corrected_marker_level.csv`
- Validation report: `farm_outputs/recomputed_170_thresholds/validation/recomputed_170_validation_report.md`
- Manuscript replacement draft: `farm_outputs/recomputed_170_thresholds/manuscript_outputs/surgical_result_audit/recomputed_170_manuscript_replacement_text.md`

The public manuscript frame is component-expanded. The marker-level frame contains the 170 binary model predictors. Do not use the marker-level count as a direct replacement for the manuscript's old `88 changes at 81 sites` claim.

## High-Confidence Headline Replacements

| Old manuscript value | Recomputed value | Frame / source | Wording note |
|---|---:|---|---|
| 354 Supplementary File 1 rows | 394 rows | public component-expanded SF1 | Supplement row count changes. |
| 172 amino acid changes with any evidence | 215 unique substitutions with any evidence | public component-expanded SF1 | Prefer `215 unique PBP substitutions`; public any-evidence rows are 236. |
| 147 unique PBP locations | 179 unique PBP positions | public component-expanded SF1 | Wording remains applicable if `locations` means numeric PBP positions. |
| 51.6% of 285 variable positions | 62.8% of 285 variable positions | public component-expanded SF1 | Number update only. |
| 88 changes at 81 sites supported by multiple methods | 117 component-expanded substitutions at 108 PBP positions | public component-expanded SF1 | Wording must distinguish component-expanded substitutions from marker-level predictors. |
| PBP split 18 / 24 / 39 positions | 25 / 42 / 41 positions | public component-expanded SF1 | PBP1A/PBP2B/PBP2X order. |
| Marker-level multi-method count not stated | 86 of 170 markers | marker-level audit table | Use only when explicitly labelled marker-level. |
| PBP2X I371T identified in all four analyses | Still true | public SF1 and marker-level table | Wording remains applicable. |
| Effect range 0.33-3.079 MIC units | 0.374-3.187 MIC units | public multi-method rows | Number update only. |
| Evidence counts VS/S/M/W/None 1/5/82/105/161 | 1/5/111/119/158 | public component-expanded SF1 | Public supplement evidence distribution changes. |
| Strong substitutions | same five PBP2X substitutions | Table 3 top-20 source | Strong set remains E320K, A491V, T490S, A369V, M343T. |

## Main Manuscript Passages Found By PDF Extraction

| Extracted location | Old wording/value | Recomputed replacement | Applicability |
|---|---|---|---|
| Manuscript lines ~515-517 | `172 amino acid changes, at 147 unique PBP locations ... 51.6%` | `215 unique PBP substitutions at 179 unique PBP positions ... 62.8%` | Wording mostly applicable, but `unique substitutions` is clearer than `changes` because public rows can be component-expanded. |
| Manuscript lines ~517-519 | `88 changes at 81 sites ... 18 in PBP1A, 24 in PBP2B, 39 in PBP2X` | `117 component-expanded substitutions at 108 positions ... 25 PBP1A, 42 PBP2B, 41 PBP2X` | Must reframe as component-expanded substitutions; optionally add separate `86 of 170 markers`. |
| Manuscript line ~519 | `PBP2X I371T identified in all four analyses` | unchanged | Still applicable. |
| Manuscript lines ~519-520 | `0.33-3.079 MIC units` | `0.374-3.187 MIC units` | Number update only. |
| Manuscript lines ~520-522 | `Most of the 88 substitutions showed both additive and epistatic effects (89.8%)` | not safely pinned to one recomputed value yet | Do not reuse exact percentage until the old 89.8% recipe is reproduced. Candidate frames include `108/117 = 92.3%` public multi rows with mvLMM display plus epistasis, or `73/86 = 84.9%` marker-level multi-evidence markers with additive plus epistasis. |
| Manuscript lines ~525-526 | `specific effects ... 59.08%, rather than common effects ... 13.63%` | needs exact denominator check | Direction appears still true, but the exact replacement depends on reproducing the old Rmd effect-size categorisation. Marker-level additive evidence gives `67.39%` one-axis specific and `9.78%` same-direction common effects; combined additive-p-threshold plus epistasis-p-threshold universe gives `54.04%` one-axis and `3.41%` same-direction common effects. |
| Manuscript lines ~527-530 | examples `I371T, N444S, S531Y`, `P568S, L609S, G467L/S469N/G483A/A490S` | examples remain present in recomputed top/evidence outputs | Likely wording remains directionally valid, but exact subclass-effect values should be checked. |
| Manuscript lines ~549-553 | `Multivariate mvLMMs detected 94 substitutions ... not identified by univariate LMMs` | needs exact original-frame derivation | Do not update with a guessed value. Marker-level strict display comparison gives mvLMM-only `59` markers; public component-expanded display gives mvLMM-only `153` rows, which is not necessarily the old manuscript frame. |
| Manuscript Table 3, lines ~574-625 | old top-20 values and membership | regenerate from `recomputed_170_table3_top20.csv` | Caption mostly remains applicable, but consider `evidence streams` instead of `methods`. Top row is now I371T beta `3.187`, interactions `40`; top-20 membership and row values changed. |
| Manuscript lines ~657-659 | `identify 88 PBP substitutions` | `117 component-expanded PBP substitutions` or `86 markers`, depending intended framing | Wording must be clarified; old wording alone is ambiguous. |
| Manuscript lines ~700-704 | `0.33-3.079`, `147 of 285 (51.6%)`, `81 positions` | `0.374-3.187`, `179 of 285 (62.8%)`, `108 positions` | Number updates required. |
| Manuscript lines ~729-734 | `most substitutions showed subclass-specific effects (59.08%)` | needs exact denominator check | Direction likely remains true; old percentage should not be retained. |

## Supplementary Information Passages Found By PDF Extraction

| Extracted location | Old wording/value | Recomputed replacement / status | Applicability |
|---|---|---|---|
| SI lines ~1709-1712 | `157 independent tests` | `170 independent tests` | Number update. |
| SI lines ~1712-1713 | `97 ... associated ... exceeding one log2 MIC unit ... p < 0.05` | likely `95` marker-level additive tests using `pv20_adj_galwey <= 0.05` and joint effect `>= 1` | Needs independent verification; wording otherwise applicable. |
| SI lines ~1713-1715 | lowest p-value `5.88 x 10-4` after 100 permutations | additive threshold is raw `0.0001418518635803`, Galwey-adjusted `0.0039718521802484` | Must state raw vs adjusted explicitly. |
| SI lines ~1715-1717 | `89 of the 97 ... below threshold`; PBP split `25/27/37` | likely `92` below recomputed threshold; PBP split `27/26/39` | Needs independent verification from marker-level table. |
| SI line ~1717 | `16 changes` over two MIC units | likely `18` additive threshold-passing markers over two MIC units | Needs independent verification. |
| SI lines ~1724-1726 | `157 independent changes`; `3542 interactions` | `170 independent changes`; `4052 interactions` | Number update. |
| SI lines ~1726-1728 | `2129` significant interactions; `p < 7.62 x 10-4` | `3147` p-threshold-only interactions; adjusted threshold `0.0016750567677524422` | Wording should distinguish p-threshold-only interactions from lower-bound supported interactions. |
| SI lines ~1729-1731 | `1634 interactions remained` after lower-bound filter | `1926` supported interactions | Number update. |
| SI line ~1731 | `138 of 285 variant PBP locations (48.4%)` | likely `159` component-expanded positions (`55.8%`) or `125` marker-level positions | Must resolve exact old frame before replacement. |
| SI lines ~1731-1733 | `30 substitutions ... >40 other positions`; dense split `1A 36.2%, 2X 42%` | marker-level: `39` markers with >40 supported interactions; PBP1A `14` (35.9%), PBP2B `6` (15.4%), PBP2X `19` (48.7%) | Needs exact frame check. |
| SI lines ~1762-1768 | evidence ranking `172`, `147`, `51.6%`, `88`, `81`, `18/24/39`, `89.8%` | same replacements as main Results paragraph; `89.8%` unresolved | Requires update and frame clarification. |
| SI Table S14 lines ~1774-1826 | old top-20 table | replace from `recomputed_170_table3_top20.csv` | Table must be regenerated. |
| SI lines ~1845-1855 | `172 substitutions and 2129 epistatic interactions`; effect-direction counts/percentages | numbers change; exact old denominator must be reproduced before final replacement | Do not retain old percentages. |
| SI lines ~1875-1884 | overlap with previous work, including `97 positions not highlighted` | likely affected by new positions, but external previous-work comparison source must be located | Mark as cannot determine unless overlap inputs/figure source are available. |
| SI lines ~1904-1917 | mvLMM vs uvLMM counts: extracted `82` plus `36/53/4` | extraction/source is internally inconsistent; recomputed values need exact frame | Do not reuse. Marker-level strict display gives `104 either`, `33 both`, `59 mvLMM-only`, `12 uvLMM-only`. |

## Wording Decisions To Make Before Manuscript Revision

- Prefer `evidence streams` over `methods` where additive mvLMM and epistasis LMM are being counted as separate streams within related model families.
- Use `component-expanded substitutions` for manuscript/Supplementary File 1 counts such as `117`.
- Pair manuscript-frame counts with the marker-level count only when useful: `117 component-expanded substitutions at 108 positions; separately, 86 of 170 model markers had two or more evidence streams`.
- Avoid quoting revised `89.8%`, `59.08%`, `13.63%`, `94`, or SI S3.3 values until the original denominator/recipe is reproduced.
- State exact thresholds in methods/provenance, with rounded display only where explicitly labelled.
