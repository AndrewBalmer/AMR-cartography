# Manuscript folder — submission package

This folder holds the final manuscript and supplementary materials
accompanying the AMR Cartography paper.

## Contents

- `Manuscript.docx` / `Manuscript.pdf` — final manuscript.
- `Supplementary_Information.docx` / `Supplementary_Information.pdf` — final
  supplementary information.
- `Supplementary_File_1.csv` — public Supplementary File 1: the
  component-expanded, 170-marker-corrected, recomputed-threshold evidence
  table, including structure-context annotation columns (position, motif
  interval, nearest motif by sequence/3D, distances, prior evidence
  source/category).
- `figures/` — main-text figures extracted at full embedded resolution from
  `Manuscript.docx` (`Figure1.jpg`, `Figure2.jpg`, `Figure3.jpg`; 2048 px wide,
  ~289 dpi at 180 mm column width). Supplied as separate submission assets;
  see the note below.
- `source_data/` — source data and figures underlying the supplementary
  tables/figures (see below).

PDF exports of `Manuscript.docx`/`Supplementary_Information.docx` are
generated directly from Word to preserve table/figure layout.

## `source_data/`

Source data/figures for the supplementary tables and figures, matching the
manuscript's figure numbering (S18 = epistasis, S19 = effect-size axis
categorisation, S20 = prior-work overlap):

- `Table3_and_SupplementaryTableS14_top20.csv`
- `FigureS18_epistatic_interaction_LMM.{pdf,png}`
- `FigureS19_effect_size_axis_categorisation.{pdf,png}`
- `FigureS20_prior_work_overlap_venn.{pdf,png}`
- `StructureContext_distance_histogram.png` + `StructureContext_legend.md`

Supplementary Figure S17 (clustering) is not part of this set — clustering
evidence is unaffected by the 157→170 marker correction, so it is unchanged
from the original analysis and generated elsewhere in `analysis/`.

## Reproducing these files

The recomputed 170-marker mvLMM/uvLMM/epistasis analysis behind
`Supplementary_File_1.csv` and the `source_data/` figures is documented in
`analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/`.

## Main-text figure assets

`figures/` holds Figures 1-3 exactly as embedded in `Manuscript.docx`. They are
2048 px wide JPEGs, which is ~289 dpi at a 180 mm full-page width and ~612 dpi
at an 85 mm single-column width. That clears the 300 dpi bar at single-column
size but sits marginally under it full-width, and they are lossy raster rather
than vector. For final production, regenerate the panels as vector PDF/EPS from
the plotting scripts under `analysis/` (Figure 1 and 2 panels from
`01-Phenotype_and_map_analyses/`, Figure 3 panels from
`02-Genotype_to_phenotype_analyses/`) and re-assemble the composites.
