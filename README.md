
# AMR Cartography

**A Generalisable Framework for Studying Multivariate Drug Resistance**

Andrew J. Balmer\*, Gemma G. R. Murray, Stephanie Lo, Olivier Restif⍅ and Lucy A. Weinert⍅

\*Corresponding author; ⍅ Joint senior author

Contact email: ab69@sanger.ac.uk

---

## Overview

AMR Cartography is a toolkit for analysing and visualising multivariate antibiotic susceptibility profiles (e.g., MICs to several drugs). It can be used to study population level trends in multivariate phenotypes, and link changes in phenotype to underlying genetic variation. It adapts multidimensional scaling (MDS/SMACOF) and multivariate linear mixed models (mvLMM) to:

* build low‑dimensional phenotype maps from log₂‑MIC profiles,
* assess goodness‑of‑fit (cross‑validation, dimensionality tests),
* integrate censored/missing values and even categorical susceptibility,
* associate PBP substitutions with multivariate phenotypes (including epistasis).

This repo contains the analysis scripts to reproduce the figures/tables in the associated manuscript and instructions on to run the core AMR Cartography workflow on similar datasets.

---

## Repository structure

**Brief overview:**

```
AMR-cartography/
├─ analysis/
│  ├─ 01-Phenotype_and_map_analyses/
│  └─ 02-Genotype_to_phenotype_analyses/
├─ manuscript/
│  ├─ Manuscript.pdf
│  ├─ Supplementary_Information.pdf
│  └─ Supplementary_File_1.csv
├─ .gitignore
└─ README.md (this file)
```

---

## Installation & requirements

AMR Cartography analyses are primarily in R, although the mvLMM steps use Python.

### R

* R ≥ 4.2 (4.3+ recommended)
* Primary packages:

  * Core: `tidyverse`, `data.table`, `readr`, `stringr`, `purrr`, `ggplot2`, `patchwork`, `cowplot`
  * MDS: `smacof`
  * R Markdown: `rmarkdown`, `knitr`

---

## Data

* The analysis uses the Active Bacterial Core surveillance dataset of *S. pneumoniae* isolates with MICs for six β‑lactams and PBP TPD sequences.
* **Inputs are referenced/ingested by the scripts**; large raw data and intermediate data files are not included in the repo. Instead, the scripts document how to acquire these using DOIs from relevant papers, or in the case of the intermediate files, how to derive them.

---

## Citation

If you use **AMR Cartography** (methods or code), **please cite**:

> Balmer AJ, Murray GGR, Lo S, Restif O, Weinert LA. *Antimicrobial Resistance Cartography: A Generalisable Framework for Studying Multivariate Drug Resistance*.
>
> **Preprint**: https://doi.org/10.1101/2025.09.12.675231
> **Journal article**: *link coming soon*

A BibTeX entry will be provided once the DOI is available.

---

## License

This project is released under the **MIT License** (permissive).
**Attribution**: in academic or public outputs that use the methods or scripts here, please include the citation above.

See [`LICENSE`](#) (to be added) for full terms.

---

## Funding & acknowledgements

This work was supported by the Biotechnology and Biological Sciences Research Council (BBSRC) (Student Project ID 2113638), BBSRC Doctoral Training Partnership.
We also acknowledge the CDC/ABC programme for making data available, and the colleagues and groups listed in the manuscript’s Acknowledgements.

---

## Contact & contributions

* Lead/contact: **Andrew J. Balmer** (ab69@sanger.ac.uk)
* Issues and improvements via GitHub Issues*/PRs are welcome.
* Please open an issue for reproducibility questions or environment pinning (we can provide an `renv.lock`).

---

## Reproducibility notes

* Rmds are designed to be run independently once inputs are prepared.
* Note, some steps (e.g., bootstrap MDS or mvLMM grids) can be compute‑intensive; scripts cache intermediates (`*.rds`, `*.RData`).
* For this repo, outputs are git‑ignored; re‑running the Rmds will regenerate the figures/tables.
