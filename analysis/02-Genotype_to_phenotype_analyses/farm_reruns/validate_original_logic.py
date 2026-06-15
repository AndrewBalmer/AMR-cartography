#!/usr/bin/env python3
"""Validate exact original-logic replication before trusting corrected outputs."""

from __future__ import annotations

import argparse
import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    EXPECTED_CORRECTED_ADDITIVE_MARKERS,
    EXPECTED_CORRECTED_EPISTASIS_CANDIDATES,
    EXPECTED_OLD_ADDITIVE_MARKERS,
    EXPECTED_OLD_EPISTASIS_CANDIDATES,
    HISTORICAL_ADDITIVE_GALWEY_MEFF,
    HISTORICAL_ADDITIVE_THRESHOLD,
    HISTORICAL_UV_THRESHOLD,
    add_galwey_adjusted_pvalues,
    four_cell_interaction_metadata,
    marker_order_from_effects,
    read_csv,
)

EXPECTED_OLD_SUPPLEMENT_COUNTS = {
    "Very Strong": 1,
    "Strong": 5,
    "Moderate": 82,
    "Weak": 105,
    "Weak/No Evidence": 161,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-dir", default=Path("AMRC-repo-files/pythonProject1"), type=Path)
    parser.add_argument(
        "--new-dir",
        default=Path("AMRC-repo-files/pythonProject1-additive-production-20260507-150112"),
        type=Path,
    )
    parser.add_argument(
        "--support-dir",
        default=Path("AMRC-repo-files/Streptococcus pneumoniae analysis"),
        type=Path,
    )
    parser.add_argument(
        "--added-marker-file",
        default=Path("farm_outputs/original_logic_rebuild/additive/added_markers.csv"),
        type=Path,
    )
    parser.add_argument("--rebuilt-output-dir", default=None, type=Path)
    parser.add_argument("--out-dir", default=Path("farm_outputs/original_logic_rebuild/validation"), type=Path)
    return parser.parse_args()


def assert_equal(label: str, observed: object, expected: object, report: list[str]) -> None:
    if observed != expected:
        raise AssertionError(f"{label}: expected {expected!r}, observed {observed!r}")
    report.append(f"- PASS `{label}` = `{observed}`.")


def assert_true(label: str, condition: bool, report: list[str], detail: str = "") -> None:
    if not condition:
        raise AssertionError(f"{label} failed. {detail}")
    suffix = f" {detail}" if detail else ""
    report.append(f"- PASS `{label}`.{suffix}")


def old_supplement_from_main() -> pd.DataFrame:
    result = subprocess.run(
        ["git", "show", "main:manuscript/Supplementary_File_1.csv"],
        check=True,
        text=True,
        capture_output=True,
    )
    return pd.read_csv(io.StringIO(result.stdout))


def additive_marker_order(directory: Path) -> list[str]:
    return marker_order_from_effects(directory / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv")


def unordered_pairs(interactions: pd.Series) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for value in interactions.astype(str):
        parent_a, parent_b = value.split(":", 1)
        pairs.add(tuple(sorted((parent_a, parent_b))))
    return pairs


def validate_old_supplement(report: list[str]) -> None:
    df = old_supplement_from_main()
    assert_equal("old Supplementary_File_1 rows", len(df), 354, report)
    counts = df["Evidence"].value_counts().reindex(list(EXPECTED_OLD_SUPPLEMENT_COUNTS), fill_value=0).to_dict()
    assert_equal("old Supplementary_File_1 evidence counts", counts, EXPECTED_OLD_SUPPLEMENT_COUNTS, report)
    multi = df[df["Evidence"].isin(["Very Strong", "Strong", "Moderate"])].copy()
    assert_equal("old multi-method component rows", len(multi), 88, report)
    positions = (
        multi.assign(Loci=multi["Substitution"].astype(str).str[1:4])[["PBP", "Loci"]]
        .drop_duplicates()
        .shape[0]
    )
    assert_equal("old multi-method positions", positions, 81, report)
    assert_equal("old slash-containing public substitution rows", int(df["Substitution"].astype(str).str.contains("/").sum()), 0, report)


def validate_marker_panels(args: argparse.Namespace, report: list[str]) -> tuple[list[str], list[str], list[str]]:
    old_order = additive_marker_order(args.old_dir)
    new_order = additive_marker_order(args.new_dir)
    assert_equal("old fitted additive markers", len(old_order), EXPECTED_OLD_ADDITIVE_MARKERS, report)
    assert_equal("corrected fitted additive markers", len(new_order), EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)
    old_set = set(old_order)
    new_set = set(new_order)
    added = sorted(new_set - old_set)
    assert_true("old fitted additive set is strict subset of corrected", old_set < new_set, report)
    assert_equal("corrected-minus-old fitted markers", len(added), 13, report)

    if args.added_marker_file.exists():
        added_file = sorted(read_csv(args.added_marker_file)["marker"].astype(str).tolist())
        assert_equal("corrected-minus-old equals added_markers.csv", added, added_file, report)
    else:
        report.append(f"- WARN added marker file not found: `{args.added_marker_file}`.")

    markers = read_csv(args.new_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    assert_equal("corrected marker matrix columns", markers.shape[1], EXPECTED_CORRECTED_ADDITIVE_MARKERS, report)
    one_percent = markers.shape[0] * 0.01
    for marker in added:
        values = markers[marker]
        assert_true(f"added marker binary/non-missing: {marker}", values.notna().all() and set(values.unique()).issubset({0, 1}), report)
        present = int(values.sum())
        absent = int((1 - values).sum())
        assert_true(
            f"added marker frequency filter: {marker}",
            present > one_percent and absent > one_percent,
            report,
            detail=f"`Present={present}`, `Absent={absent}`.",
        )
        for other in markers.columns:
            if other == marker:
                continue
            other_values = markers[other]
            duplicate = bool(values.equals(other_values))
            complement = bool((values + other_values).eq(1).all())
            if duplicate or complement:
                raise AssertionError(f"Added marker {marker} duplicates/complements {other}")
    report.append("- PASS all 13 added markers are non-duplicate and non-complement within the corrected panel.")
    return old_order, new_order, added


def validate_epistasis_candidates(args: argparse.Namespace, added: list[str], report: list[str]) -> None:
    old_counts = read_csv(args.old_dir / "Counts_of_each_AA_combination_epi.csv")
    old_matrix = pd.read_csv(args.old_dir / "S.pneumo_map_test_markers_incl_2nd_order_epistatic.csv", nrows=0)
    old_p = read_csv(args.old_dir / "mvLMM_p_values_epi_pneumo.csv")
    assert_equal("old epistasis count table rows", len(old_counts), EXPECTED_OLD_EPISTASIS_CANDIDATES, report)
    assert_equal("old epistasis matrix columns", old_matrix.shape[1], EXPECTED_OLD_EPISTASIS_CANDIDATES, report)
    assert_equal("old epistasis p-value rows", len(old_p), EXPECTED_OLD_EPISTASIS_CANDIDATES, report)
    assert_true(
        "old epistasis count names match matrix columns",
        old_counts["PBP_AA_location"].astype(str).tolist() == old_matrix.columns.astype(str).tolist(),
        report,
    )

    new_markers = read_csv(args.new_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    corrected = four_cell_interaction_metadata(new_markers)
    assert_equal("corrected original-rule epistasis candidates", len(corrected), EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, report)

    old_pairs = unordered_pairs(old_counts["PBP_AA_location"])
    corrected_pairs = unordered_pairs(corrected["interaction"])
    assert_true("old unordered epistasis pairs subset corrected", old_pairs <= corrected_pairs, report)
    corrected_only = corrected_pairs - old_pairs
    assert_equal("corrected-only original-rule epistasis pairs", len(corrected_only), 510, report)
    added_set = set(added)
    bad_new = [pair for pair in corrected_only if pair[0] not in added_set and pair[1] not in added_set]
    assert_equal("corrected-only pairs not involving added marker", len(bad_new), 0, report)


def validate_scoring_frames(args: argparse.Namespace, report: list[str]) -> None:
    old_effects = read_csv(args.support_dir / "Sub_effect_sizes_mv_pneumo.csv")
    adjusted_hits = int(
        (
            (pd.to_numeric(old_effects["pv20_adj_galwey"], errors="coerce") <= HISTORICAL_ADDITIVE_THRESHOLD)
            & (pd.to_numeric(old_effects["Joint_effsize"], errors="coerce") >= 1)
        ).sum()
    )
    raw_hits = int(
        (
            (pd.to_numeric(old_effects["pv20"], errors="coerce") <= HISTORICAL_ADDITIVE_THRESHOLD)
            & (pd.to_numeric(old_effects["Joint_effsize"], errors="coerce") >= 1)
        ).sum()
    )
    assert_equal("old additive adjusted evidence hits", adjusted_hits, 78, report)
    assert_equal("old additive raw-threshold hits", raw_hits, 84, report)
    assert_true("raw additive threshold would change old evidence", raw_hits != adjusted_hits, report)

    old_p = read_csv(args.old_dir / "uniLMM_p_val_normal_MIC_pneumo.csv")
    old_e = read_csv(args.old_dir / "uniLMM_effect_normal_MIC_pneumo.csv")
    old_cand = old_e[old_e["effect_type"].eq("candidate")][["trait", "effect_name"]].reset_index(drop=True)
    uv = old_cand.rename(columns={"trait": "drug", "effect_name": "marker"})
    uv["pv20"] = old_p["pv20"].to_numpy()
    uv = add_galwey_adjusted_pvalues(uv, meff=HISTORICAL_ADDITIVE_GALWEY_MEFF)
    raw_uv_tests = int((uv["pv20"] <= HISTORICAL_UV_THRESHOLD).sum())
    adjusted_uv_tests = int((uv["pv20_adj_galwey"] <= HISTORICAL_UV_THRESHOLD).sum())
    raw_uv_markers = int(uv.loc[uv["pv20"] <= HISTORICAL_UV_THRESHOLD, "marker"].nunique())
    adjusted_uv_markers = int(uv.loc[uv["pv20_adj_galwey"] <= HISTORICAL_UV_THRESHOLD, "marker"].nunique())
    assert_equal("old uvLMM raw significant marker-drug tests", raw_uv_tests, 127, report)
    assert_equal("old uvLMM adjusted significant marker-drug tests", adjusted_uv_tests, 66, report)
    assert_equal("old uvLMM raw significant markers", raw_uv_markers, 61, report)
    assert_equal("old uvLMM adjusted significant markers", adjusted_uv_markers, 40, report)


def validate_rebuilt_outputs(output_dir: Path, report: list[str]) -> None:
    public_path = output_dir / "Supplementary_File_1.csv"
    marker_path = output_dir / "Supplementary_File_1_corrected_marker_level.csv"
    if not public_path.exists() or not marker_path.exists():
        report.append(f"- SKIP rebuilt output checks; expected files are not both present in `{output_dir}`.")
        return
    public = read_csv(public_path)
    marker = read_csv(marker_path)
    assert_true("rebuilt public table is component-expanded", len(public) >= len(marker), report)
    multi = public[public["Evidence"].isin(["Very Strong", "Strong", "Moderate"])].copy()
    positions = (
        multi.assign(Loci=multi["Substitution"].astype(str).str[1:4])[["PBP", "Loci"]]
        .drop_duplicates()
        .shape[0]
    )
    counts = public["Evidence"].value_counts().to_dict()
    marker_counts = marker["Evidence"].value_counts().to_dict()
    report.extend(
        [
            f"- INFO rebuilt component-expanded rows: `{len(public)}`.",
            f"- INFO rebuilt marker-level rows: `{len(marker)}`.",
            f"- INFO rebuilt multi-method component rows: `{len(multi)}`.",
            f"- INFO rebuilt multi-method component positions: `{positions}`.",
            f"- INFO rebuilt public evidence counts: `{json.dumps(counts, sort_keys=True)}`.",
            f"- INFO rebuilt marker-level evidence counts: `{json.dumps(marker_counts, sort_keys=True)}`.",
        ]
    )


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = [
        "# Original Logic Validation Report",
        "",
        "This report is generated by `validate_original_logic.py` and must pass before corrected manuscript-facing outputs are trusted.",
        "",
        "## Checks",
        "",
    ]

    validate_old_supplement(report)
    _, _, added = validate_marker_panels(args, report)
    validate_epistasis_candidates(args, added, report)
    validate_scoring_frames(args, report)
    if args.rebuilt_output_dir is not None:
        validate_rebuilt_outputs(args.rebuilt_output_dir, report)

    report.extend(
        [
            "",
            "## Result",
            "",
            "All required original-logic validation checks passed.",
        ]
    )
    text = "\n".join(report) + "\n"
    (args.out_dir / "original_logic_validation_report.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
