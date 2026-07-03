#!/usr/bin/env python3
"""Merge recomputed 170-marker additive mvLMM chunks and thresholds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import EXPECTED_CORRECTED_ADDITIVE_MARKERS, add_galwey_adjusted_pvalues, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--threshold-dir", required=True, type=Path)
    parser.add_argument("--galwey-meff", required=True, type=float)
    parser.add_argument("--expected-markers", default=EXPECTED_CORRECTED_ADDITIVE_MARKERS, type=int)
    parser.add_argument("--expected-permutations", default=100, type=int)
    parser.add_argument("--old-null-file", type=Path)
    return parser.parse_args()


def read_many(files: list[Path]) -> pd.DataFrame:
    return pd.concat([pd.read_csv(path) for path in files], ignore_index=True) if files else pd.DataFrame()


def non_ok_count(df: pd.DataFrame) -> int:
    if df.empty or "fit_status" not in df.columns:
        return 0
    return int((df["fit_status"].notna() & (df["fit_status"] != "ok")).sum())


def lowest_min_p_threshold(permutations: pd.DataFrame, *, meff: float) -> tuple[pd.DataFrame, dict[str, object]]:
    required = {"repeat_p_index", "pv20"}
    missing = required - set(permutations.columns)
    if missing:
        raise ValueError(f"Additive permutation chunks are missing columns: {sorted(missing)}")
    adjusted = add_galwey_adjusted_pvalues(permutations, meff=meff)
    minima = (
        adjusted.groupby("repeat_p_index", dropna=False)
        .agg(
            raw_min_pv20=("pv20", "min"),
            adjusted_min_pv20=("pv20_adj_galwey", "min"),
            permutation_seed=("permutation_seed", "first"),
            rows=("pv20", "size"),
        )
        .reset_index()
        .sort_values("raw_min_pv20", kind="mergesort")
    )
    threshold_row = minima.iloc[0].to_dict()
    summary = {
        "threshold_policy": "lowest permutation minimum",
        "threshold_source": "additive random-phenotype null",
        "raw_threshold": float(threshold_row["raw_min_pv20"]),
        "adjusted_threshold": float(threshold_row["adjusted_min_pv20"]),
        "threshold_repeat_p_index": int(threshold_row["repeat_p_index"]),
        "threshold_permutation_seed": None
        if pd.isna(threshold_row.get("permutation_seed"))
        else int(threshold_row["permutation_seed"]),
    }
    return minima, summary


def old_additive_calibration(path: Path | None) -> dict[str, object]:
    if path is None or not path.exists():
        return {
            "old_additive_calibration_available": False,
            "old_additive_calibration_note": "old null file was not supplied or does not exist",
        }
    old = read_csv(path)
    minima = old.groupby("repeat_p_index", dropna=False)["pv20"].min().sort_values(kind="mergesort")
    return {
        "old_additive_calibration_available": True,
        "old_additive_lowest_raw_min_p": float(minima.iloc[0]),
        "old_additive_expected_lowest_raw_min_p": 0.0005875153233810342,
        "old_additive_reproduces_expected": abs(float(minima.iloc[0]) - 0.0005875153233810342) < 1e-15,
        "old_additive_permutations": int(old["repeat_p_index"].nunique()),
        "old_additive_rows": int(len(old)),
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.threshold_dir.mkdir(parents=True, exist_ok=True)

    observed_files = sorted(
        path for path in args.chunk_dir.glob("additive_p_values.chunk_*.csv") if "smoke" not in path.name
    )
    effect_files = sorted(
        path for path in args.chunk_dir.glob("additive_effect_sizes.chunk_*.csv") if "smoke" not in path.name
    )
    perm_files = sorted(
        path for path in args.chunk_dir.glob("additive_perm_p_values.chunk_*.csv") if "smoke" not in path.name
    )
    if not observed_files:
        raise FileNotFoundError(f"No observed additive chunks found in {args.chunk_dir}")

    observed = read_many(observed_files).sort_values("marker_index", kind="mergesort")
    effects = read_many(effect_files).sort_values(["marker_index", "test"], kind="mergesort")
    permutations = read_many(perm_files).sort_values(["repeat_p_index", "marker_index"], kind="mergesort")

    expected_effect_rows = args.expected_markers * 4
    expected_permutation_rows = args.expected_markers * args.expected_permutations
    if len(observed) != args.expected_markers:
        raise AssertionError(f"Expected {args.expected_markers} observed additive rows, found {len(observed)}")
    if len(effects) != expected_effect_rows:
        raise AssertionError(f"Expected {expected_effect_rows} observed additive effect rows, found {len(effects)}")
    if len(permutations) != expected_permutation_rows:
        raise AssertionError(f"Expected {expected_permutation_rows} additive permutation rows, found {len(permutations)}")
    if permutations["repeat_p_index"].nunique() != args.expected_permutations:
        raise AssertionError("Additive permutation repeat count is incomplete")
    if non_ok_count(observed) or non_ok_count(permutations):
        raise AssertionError("Primary additive outputs contain non-ok fit rows")

    minima, threshold_summary = lowest_min_p_threshold(permutations, meff=args.galwey_meff)
    observed_adj = add_galwey_adjusted_pvalues(observed.rename(columns={"pv20": "pv20_raw"}), raw_column="pv20_raw", meff=args.galwey_meff)
    permutations_adj = add_galwey_adjusted_pvalues(
        permutations.rename(columns={"pv20": "pv20_raw"}), raw_column="pv20_raw", meff=args.galwey_meff
    )

    observed.to_csv(args.out_dir / "mvLMM_p_values_normal_pneumo_low_freq_vars.csv", index=False)
    effects.to_csv(args.out_dir / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv", index=False)
    permutations.to_csv(args.out_dir / "mvLMM_p_values_normal_pneumo_random_phenotype_FWAS.csv", index=False)
    observed_adj.to_csv(args.out_dir / "mvLMM_p_values_normal_pneumo_low_freq_vars_adjusted.csv", index=False)
    permutations_adj.to_csv(args.out_dir / "mvLMM_p_values_normal_pneumo_random_phenotype_FWAS_adjusted.csv", index=False)
    minima.to_csv(args.threshold_dir / "additive_permutation_minima.csv", index=False)

    summary = {
        "galwey_meff": float(args.galwey_meff),
        "observed_rows": int(len(observed)),
        "effect_rows": int(len(effects)),
        "permutation_rows": int(len(permutations)),
        "permutations": int(permutations["repeat_p_index"].nunique()),
        "observed_non_ok_fit_rows": non_ok_count(observed),
        "permutation_non_ok_fit_rows": non_ok_count(permutations),
        **threshold_summary,
        **old_additive_calibration(args.old_null_file),
    }
    (args.threshold_dir / "additive_threshold_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
