#!/usr/bin/env python3
"""Merge full corrected-panel exact uvLMM chunks.

Input: all uvLMM marker-drug chunks. Output: 1,020 marker-drug tests plus a
marker-level support table used only for the public `Sig. mvLMM/uvLMM` display
column. This script must not add uvLMM as a fifth evidence stream.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import EXPECTED_CORRECTED_ADDITIVE_MARKERS, HISTORICAL_UV_THRESHOLD, add_galwey_adjusted_pvalues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--threshold-dir", required=True, type=Path)
    parser.add_argument("--galwey-meff", required=True, type=float)
    parser.add_argument("--threshold", default=HISTORICAL_UV_THRESHOLD, type=float)
    parser.add_argument("--expected-markers", default=EXPECTED_CORRECTED_ADDITIVE_MARKERS, type=int)
    parser.add_argument("--expected-drugs", default=6, type=int)
    return parser.parse_args()


def read_many(files: list[Path]) -> pd.DataFrame:
    return pd.concat([pd.read_csv(path) for path in files], ignore_index=True) if files else pd.DataFrame()


def non_ok_count(df: pd.DataFrame) -> int:
    if df.empty or "fit_status" not in df.columns:
        return 0
    return int((df["fit_status"].notna() & (df["fit_status"] != "ok")).sum())


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.threshold_dir.mkdir(parents=True, exist_ok=True)

    pvalue_files = sorted(
        path for path in args.chunk_dir.glob("uniLMM_exact_p_values.chunk_*.csv") if "smoke" not in path.name
    )
    effect_files = sorted(
        path for path in args.chunk_dir.glob("uniLMM_exact_effect_sizes.chunk_*.csv") if "smoke" not in path.name
    )
    if not pvalue_files:
        raise FileNotFoundError(f"No uvLMM chunks found in {args.chunk_dir}")

    combined = read_many(pvalue_files).sort_values(["marker_index", "drug"], kind="mergesort")
    effects = read_many(effect_files).sort_values(["marker_index", "drug", "effect_type"], kind="mergesort")
    expected_tests = args.expected_markers * args.expected_drugs
    if len(combined) != expected_tests:
        raise AssertionError(f"Expected {expected_tests} uvLMM rows, found {len(combined)}")
    if combined["marker"].nunique() != args.expected_markers:
        raise AssertionError("uvLMM marker count is incomplete")
    if non_ok_count(combined):
        raise AssertionError("Primary uvLMM output contains non-ok fit rows")

    combined = add_galwey_adjusted_pvalues(combined, meff=args.galwey_meff)
    combined["uvLMM_exact_significant"] = combined["pv20_adj_galwey"] <= args.threshold

    support_rows = []
    for marker, group in combined.groupby("marker", dropna=False):
        sig = sorted(group.loc[group["uvLMM_exact_significant"], "drug"].astype(str).unique())
        support_rows.append(
            {
                "marker": marker,
                "uvLMM_exact_available": True,
                "uvLMM_exact_n_drugs": len(sig),
                "uvLMM_exact_drugs": ";".join(sig),
                "uvLMM_exact_source": "recomputed_full_170",
                "uvLMM_exact_threshold": args.threshold,
                "uvLMM_exact_galwey_meff": args.galwey_meff,
            }
        )
    support = pd.DataFrame(support_rows)

    combined.to_csv(args.out_dir / "uniLMM_exact_170_marker_drug_results.csv", index=False)
    effects.to_csv(args.out_dir / "uniLMM_exact_170_marker_effect_sizes.csv", index=False)
    support.to_csv(args.out_dir / "uniLMM_exact_170_marker_support.csv", index=False)
    summary = {
        "galwey_meff": float(args.galwey_meff),
        "uv_threshold": float(args.threshold),
        "rows": int(len(combined)),
        "markers": int(support["marker"].nunique()),
        "drugs": int(args.expected_drugs),
        "markers_with_any_uv_support": int((support["uvLMM_exact_n_drugs"] > 0).sum()),
        "non_ok_fit_rows": non_ok_count(combined),
    }
    (args.threshold_dir / "uvlmm_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
