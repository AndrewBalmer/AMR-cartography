# The recomputed 170-marker workflow, explained plainly

A plain-language companion to the technical `README.md`. It answers three
questions: **why does this workflow exist, what does each part do (in order),
and how can I trust it without reading every line of code?**

---

## The one-paragraph version

The published analysis tested **157** PBP substitution markers. We later found
**13 valid markers had been left out**, so the corrected panel has **170**.
Instead of editing the original scripts (risky, hard to check), this workflow is
a **separate, self-contained re-run** of the same analysis on the corrected 170
markers. Before it touches anything new, it **re-proves it can reproduce the old
published numbers exactly** — so when the new numbers differ, you know the *only*
thing that changed is the 157→170 correction, not a coding mistake. It runs the
heavy model-fitting on the compute farm, re-derives the significance thresholds
the same way the original did, and rebuilds the evidence table. Every stage
**checks itself** and stops if anything is wrong.

---

## How to trust it *without* reading the code

You do **not** need to read every script. The workflow was built to prove itself
with two commands. Run these from the repo root:

```bash
# 1. Does it reproduce the OLD published analysis exactly?
python3 analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/validate_original_logic.py

# 2. Are the NEW 170-marker outputs complete and consistent?
python3 analysis/02-Genotype_to_phenotype_analyses/recomputed_170_workflow/validate_recomputed_outputs.py \
  --results-dir analysis_outputs/recomputed_170_thresholds
```

Each prints a list of checks that either say **`PASS`** or the script **stops and
tells you exactly what failed**. If both end with "All … checks passed," the
pipeline is internally consistent end-to-end.

There are **three independent self-checks** built in, and this is really the
heart of why it's trustworthy:

1. **Golden test** — it re-derives the old paper's exact numbers (the 354-row
   table, the 1 / 5 / 82 / 105 / 161 evidence counts, 157 markers, 3,542
   epistasis candidates). If it can reproduce the old results, its machinery is
   faithful.
2. **Engine equivalence** — the farm uses a faster maths engine than the
   original. A check confirms the fast engine gives the *same answer* as the
   original slow one (agreement to ~8 decimal places).
3. **Completeness** — it confirms every model finished with **no failed fits**
   and the exact expected number of rows.

On top of that: the **original scripts' logic was never changed** (only a note
was added), and the headline manuscript numbers were independently re-derived
from the raw outputs and matched.

---

## The flow, step by step

```
  INPUTS:  170 markers  ·  relatedness matrix  ·  phenotype map (MIC position)
                │
   ┌────────────┴─────────────────────────────────────────────┐
   │  0. common.py .................... shared stats toolbox    │
   │  1. validate_original_logic.py ... ✅ reproduce OLD paper  │  ← trust gate
   └────────────┬─────────────────────────────────────────────┘
                │
   BUILD THE TEST LISTS + STRICTNESS
   2. generate_corrected_epistasis_interactions.py .. 4,052 marker-pairs
   3. compute_recomputed_meff.R ..................... how strict (28 / 38)
                │
   RUN THE MODELS  (heavy → farm, in chunks; each also run on 100 shuffles)
   4. run_additive_chunk → merge_additive_chunks ..... 170 single markers
   5. run_exact_unilmm_chunk → merge ................. 170 × 6 drugs (compare only)
   6. run_epistasis_chunk (+perms) → merge_epistasis . 4,052 pairs
                │
   SCORE + ASSEMBLE
   7. build_corrected_evidence_table.py .. score each substitution 0–4
                │
   OUTPUTS:  Supplementary File 1 (394 rows)  +  marker-level table (170)
                │
   CHECK + DRAW
   8. validate_recomputed_outputs.py, verify_epistasis_engine_equivalence.py
      generate_recomputed_supplement_figures_original_style.R
```

### The inputs (what goes in)
- **Marker table** — 170 columns, one per PBP substitution; each isolate is a 1
  (has it) or 0 (doesn't).
- **Relatedness matrix** — how genetically related the isolates are. This lets
  the model avoid mistaking "related strains happen to share a mutation" for
  "the mutation causes resistance."
- **Phenotype map** — each isolate's position on the 2D resistance "cartography"
  map (this *is* the MIC phenotype, in map coordinates).

### 0. `common.py` — the shared toolbox
Not run on its own; it's the parts bin every other script imports. It defines the
multiple-testing correction, keeps the **raw** and **adjusted** p-values as two
separate columns so they can never be mixed up, defines the rule for which
marker-pairs are eligible epistasis candidates, and holds the fast model engine.

### 1. `validate_original_logic.py` — the trust gate
Runs **first**, before any new data. It rebuilds the old published results and
confirms they match, and checks the 13 new markers are legitimate (truly binary,
common enough to test, not duplicates of existing markers). If anything is off,
it stops here.

### 2. `generate_corrected_epistasis_interactions.py` — list the pairs to test
Builds the list of marker-pairs to test for epistasis, using the **exact same
rule** the original used (a pair is eligible only if all four combinations of
present/absent for the two markers occur in ≥1% of isolates). Result: **4,052**
candidate pairs (was 3,542; the extra 510 all involve at least one new marker).

### 3. `compute_recomputed_meff.R` — set the strictness
Works out the "effective number of independent tests" (**28** for the 170
markers, **38** for the 4,052 pairs). Because many markers are correlated, this
is smaller than the raw count, and it controls how strict the significance cutoff
must be. Uses the standard `poolr` R library.

### 4. Additive scan — `run_additive_chunk.py` → `merge_additive_chunks.py`
For each of the 170 markers, fit the model: **"does this single substitution
shift the resistance phenotype?"** It's run once for real, then **100 more times
on shuffled data** (permutations) to learn what "significant by pure chance"
looks like. The significance threshold is the strongest chance result across
those 100 shuffles — the same rule the original used. The heavy fitting happens
in chunks on the farm; `merge` stitches them together and, as a check,
reproduces the old threshold to 12 digits.

### 5. Univariate scan — `run_exact_unilmm_chunk.py` → merge
The simpler "one drug at a time" version: 170 markers × 6 drugs = **1,020**
tests. This exists **only for comparison/display** and is deliberately **not**
one of the evidence categories that count toward a substitution's score.

### 6. Epistasis scan — `run_epistasis_chunk.py` (+ permutations) → `merge_epistasis_chunks.py`
For each of the 4,052 pairs, ask: **"do these two substitutions together do
something beyond what each does alone?"** Again once for real + 100 permutations
to set the threshold. A pair only "counts" if it passes the significance
threshold **and** its combined effect is clearly bigger than its error bar. When
a pair counts, credit is given back to **both** markers in it.

### 7. `build_corrected_evidence_table.py` — score and assemble
Each substitution gets an evidence score from **0 to 4** by counting how many of
**four independent methods** flagged it:

> single-substitution effect · clustering · additive mvLMM · epistasis

- **4 = Very Strong** (all four agree) · 3 = Strong · 2 = Moderate · 1 = Weak ·
  0 = Weak/No evidence.

It writes **two tables that are kept carefully separate** (see below).

### 8. The checks and figures
- `validate_recomputed_outputs.py` — the completeness self-check (right counts,
  no failures).
- `verify_epistasis_engine_equivalence.py` — confirms the fast engine ≡ the
  original engine.
- `generate_recomputed_supplement_figures_original_style.R` — redraws the
  supplement figures in the original style from the new numbers.

---

## The two output tables — don't mix them up

This is the single most important thing to keep straight when writing the
manuscript:

| Table | Rows | Use it for |
|---|---|---|
| **Supplementary File 1** (component-expanded, public) | **394** | Every manuscript substitution/position count |
| **Marker-level audit table** | **170** | Only when you explicitly say "marker-level" |

They count the same results two different ways (one row per substitution-component
vs. one row per model marker), so their numbers **should** differ. Quoting a
number from one while implying the other is the easiest mistake to make.

Headline public-frame numbers: 394 rows; evidence **1 / 5 / 111 / 119 / 158**
(Very Strong→none); 215 substitutions / 179 positions with any evidence; **117**
substitutions / **108** positions supported by two or more methods.

---

## If you want to sanity-check it yourself (5 minutes, no coding)

1. Run the two validator commands above — look for "All … checks passed."
2. Open the threshold file and eyeball it:
   `analysis_outputs/recomputed_170_thresholds/thresholds/recomputed_thresholds.json`
3. Open the final table and count the rows (394):
   `analysis_outputs/recomputed_170_thresholds/manuscript_outputs/Supplementary_File_1.csv`

If those three agree with this document, the pipeline is doing what it says.

---

## Mini-glossary

- **Marker / substitution** — a specific amino-acid change in a PBP; the thing
  being tested.
- **mvLMM (additive)** — the model that asks whether a marker shifts the
  phenotype, while accounting for strain relatedness.
- **Epistasis** — when two markers together do something different from the sum
  of their individual effects.
- **Permutation** — re-running on deliberately shuffled data to measure what
  "significant by chance" looks like; used to set the threshold honestly.
- **Galwey meff** — the effective number of independent tests; smaller than the
  raw count because markers are correlated. Sets how strict "significant" is.
- **Raw vs adjusted p-value** — raw is the model's output; adjusted multiplies it
  to account for many tests. The workflow always keeps both, separately.
