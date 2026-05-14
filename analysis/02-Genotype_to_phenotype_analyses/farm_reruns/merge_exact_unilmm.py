#!/usr/bin/env python3
"""Merge historical exact uvLMM and corrected-panel added-marker uvLMM evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import DEFAULT_UV_THRESHOLD, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--new-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--threshold", default=DEFAULT_UV_THRESHOLD, type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    old_p = read_csv(args.legacy_dir / "uniLMM_p_val_normal_MIC_pneumo.csv")
    old_e = read_csv(args.legacy_dir / "uniLMM_effect_normal_MIC_pneumo.csv")
    old_cand = old_e[old_e["effect_type"].eq("candidate")][["trait", "effect_name", "effsize"]].reset_index(drop=True)
    if len(old_cand) != len(old_p):
        raise ValueError("Historical uvLMM p-value/effect rows do not align")
    old = old_cand.rename(columns={"trait": "drug", "effect_name": "marker"})
    old["pv20"] = old_p["pv20"].to_numpy()
    old["source"] = "historical_157"

    added_path = args.new_dir / "uniLMM_exact_added_markers_p_values.csv"
    added = read_csv(added_path)[["marker", "drug", "pv20"]].copy()
    added["source"] = "rerun_added_13"

    combined = pd.concat([old[["marker", "drug", "pv20", "source"]], added], ignore_index=True)
    support_rows = []
    for marker, group in combined.groupby("marker", dropna=False):
        sig = sorted(group.loc[group["pv20"] < args.threshold, "drug"].astype(str).unique())
        support_rows.append(
            {
                "marker": marker,
                "uvLMM_exact_available": True,
                "uvLMM_exact_n_drugs": len(sig),
                "uvLMM_exact_drugs": ";".join(sig),
                "uvLMM_exact_source": ";".join(sorted(group["source"].unique())),
            }
        )
    support = pd.DataFrame(support_rows)
    combined.to_csv(args.out_dir / "uniLMM_exact_170_marker_drug_results.csv", index=False)
    support.to_csv(args.out_dir / "uniLMM_exact_170_marker_support.csv", index=False)
    print("combined", combined.shape)
    print("support", support.shape)
    print("markers", support["marker"].nunique())


if __name__ == "__main__":
    main()
