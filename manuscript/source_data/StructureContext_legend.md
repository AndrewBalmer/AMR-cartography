# Structure Context Annotation Legend

These annotations are post hoc biological context only. They were not used to
assign evidence categories, p-values, or effect sizes.

## Reference Structures

- PBP1A: RCSB 2C6W, chain B.
- PBP2B: RCSB 2WAF, chain A.
- PBP2X: RCSB 5OIZ, chain A.

The selected structures use direct PBP residue numbering for the residues
reported in Supplementary File 1.

## Column Definitions

- `Motif interval`: the candidate residue's linear position relative to the
  conserved active-site motif blocks in that PBP.
- `Nearest motif by sequence`: the conserved active-site motif with the smallest amino-acid
  coordinate distance to the candidate residue.
- `Sequence distance to motif (aa)`: the number of amino acids between the
  candidate residue and the nearest conserved active-site motif residue in the primary
  sequence. A value of 0 means the residue is inside a conserved motif block.
- `Nearest motif by 3D`: the conserved active-site motif containing the closest motif
  residue in the folded reference structure.
- `3D C-alpha distance to motif (Angstrom)`: the minimum Euclidean distance
  between the candidate residue C-alpha atom and any C-alpha atom in a conserved active-site
  motif. `not resolved` means the candidate residue or relevant structure
  coordinate is absent from the selected PDB chain.
- `Previous evidence category`: position-level overlap with the existing
  `CDC_GWAS_overlap_TPD.csv` source.
- `Previous evidence source`: expanded source label for the existing overlap
  file. `Chewapreecha et al. GWAS` corresponds to the GWAS column, `Li et al.
  CDC study` corresponds to the CDC column, and `Laboratory (see Supplementary
  S3.2)` corresponds to the Laboratory column.

## Interpretation Notes

Each PBP transpeptidase domain has a single active-site region organised around
three conserved motif blocks: SXXK, SXN, and KTG/KSG. The annotations identify
which motif block is closest to a candidate substitution. The nearest motif can
differ between sequence and 3D columns because residues distant in the linear
sequence can be close together after protein folding.

## Histogram Definitions

- `structure_distance_histogram_all_evaluated_positions.png`: compares Moderate+
  positions against all unique PBP positions represented in Supplementary File 1.
  This is the fairer analysis-background comparison.
- `structure_distance_histogram_all_resolved_tpd_positions.png`: compares
  Moderate+ positions against all residues in the PBP transpeptidase-domain
  ranges that are resolved in the selected reference structures, faceted into
  separate PBP1A, PBP2B, and PBP2X panels. This is a broader structural-background
  comparison. Dark-to-light red shows Very Strong, Strong, and Moderate
  positions; light pink shows Weak/No evidence positions; grey shows residues
  that were not present in the Supplementary File 1 candidate frame because
  they were invariant, filtered, or otherwise not evaluated/reported.

Histograms use counts, continuous 1-Angstrom bins, and stacked bars.
