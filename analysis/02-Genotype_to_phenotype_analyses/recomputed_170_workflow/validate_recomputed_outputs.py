#!/usr/bin/env python3
"""Validate primary recomputed 170-marker outputs before manuscript use.

Run this after merging and rebuilding outputs. It checks that all observed and
permutation families are complete, that primary fits have no non-ok rows, and
that final public/marker-level tables have the accepted shapes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import EXPECTED_CORRECTED_ADDITIVE_MARKERS, EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--out-report", type=Path)
    parser.add_argument("--expected-permutations", default=100, type=int)
    parser.add_argument("--expected-drugs", default=6, type=int)
    return parser.parse_args()


def check_equal(name: str, actual: object, expected: object, report: list[str]) -> None:
    if actual != expected:
        raise AssertionError(f"{name}: expected {expected!r}, got {actual!r}")
    report.append(f"- PASS `{name}` = `{actual}`.")


def check_true(name: str, condition: bool, report: list[str], detail: str = "") -> None:
    if not condition:
        raise AssertionError(f"{name} failed. {detail}")
    suffix = f" {detail}" if detail else ""
    report.append(f"- PASS `{name}`.{suffix}")


def non_ok_count(df: pd.DataFrame) -> int:
    if df.empty or "fit_status" not in df.columns:
        return 0
    return int((df["fit_status"].notna() & (df["fit_status"] != "ok")).sum())


def load_json(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def main() -> None:
    args = parse_args()
    report = [
        "# Recomputed 170-marker validation report",
        "",
        "Primary outputs must pass these checks before manuscript-facing use.",
        "",
        "## Checks",
        "",
    ]

    additive_dir = args.results_dir / "additive" / "merged"
    uvlmm_dir = args.results_dir / "uvlmm" / "merged"
    epi_dir = args.results_dir / "epistasis" / "merged"
    threshold_dir = args.results_dir / "thresholds"
    manuscript_dir = args.results_dir / "manuscript_outputs"

    additive = read_csv(additive_dir / "mvLMM_p_values_normal_pneumo_low_freq_vars.csv")
    additive_effects = read_csv(additive_dir / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv")
    additive_perm = read_csv(additive_dir / "mvLMM_p_values_normal_pneumo_random_phenotype_FWAS.csv")
    check_equal("additive observed rows", len(additive), EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)
    check_equal("additive observed effect rows", len(additive_effects), EXPECTED_CORRECTED_ADDITIVE_MARKERS * 4, report)
    check_equal(
        "additive permutation rows",
        len(additive_perm),
        EXPECTED_CORRECTED_ADDITIVE_MARKERS * args.expected_permutations,
        report,
    )
    check_equal("additive permutation count", additive_perm["repeat_p_index"].nunique(), args.expected_permutations, report)
    check_equal("additive non-ok rows", non_ok_count(additive) + non_ok_count(additive_perm), 0, report)

    uv = read_csv(uvlmm_dir / "uniLMM_exact_170_marker_drug_results.csv")
    check_equal("uvLMM rows", len(uv), EXPECTED_CORRECTED_ADDITIVE_MARKERS * args.expected_drugs, report)
    check_equal("uvLMM markers", uv["marker"].nunique(), EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)
    check_equal("uvLMM non-ok rows", non_ok_count(uv), 0, report)

    epi = read_csv(epi_dir / "corrected_epistasis_p_values.csv")
    epi_effects = read_csv(epi_dir / "corrected_epistasis_effect_sizes.csv")
    epi_perm = read_csv(epi_dir / "corrected_epistasis_permutation_p_values.csv")
    check_equal("epistasis observed rows", len(epi), EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, report)
    check_equal("epistasis observed effect rows", len(epi_effects), EXPECTED_CORRECTED_EPISTASIS_CANDIDATES * 2, report)
    check_equal(
        "epistasis permutation rows",
        len(epi_perm),
        EXPECTED_CORRECTED_EPISTASIS_CANDIDATES * args.expected_permutations,
        report,
    )
    check_equal("epistasis permutation count", epi_perm["permutation_seed"].nunique(), args.expected_permutations, report)
    check_equal("epistasis non-ok rows", non_ok_count(epi) + non_ok_count(epi_perm), 0, report)

    additive_summary = load_json(threshold_dir / "additive_threshold_summary.json")
    epistasis_summary = load_json(epi_dir / "corrected_epistasis_summary.json")
    meff = load_json(threshold_dir / "recomputed_meff.json")
    check_equal("additive threshold policy", additive_summary["threshold_policy"], "lowest permutation minimum", report)
    check_equal("epistasis threshold policy", epistasis_summary["threshold_policy"], "lowest permutation minimum", report)
    check_true("additive old calibration reproduced", bool(additive_summary["old_additive_reproduces_expected"]), report)
    check_equal("meff additive marker count", meff["additive_marker_count"], EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)
    check_equal("meff epistasis interaction count", meff["epistasis_interaction_count"], EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, report)

    public = read_csv(manuscript_dir / "Supplementary_File_1.csv")
    marker = read_csv(manuscript_dir / "Supplementary_File_1_corrected_marker_level.csv")
    check_true("public table is component-expanded", len(public) >= len(marker), report)
    check_equal("marker-level table rows", len(marker), EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)

    report.extend(["", "## Result", "", "All recomputed 170-marker validation checks passed."])
    out_report = args.out_report or (args.results_dir / "validation" / "recomputed_170_validation_report.md")
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
