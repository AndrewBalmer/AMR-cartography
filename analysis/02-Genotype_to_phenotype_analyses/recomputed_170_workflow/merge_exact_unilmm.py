#!/usr/bin/env python3
"""Merge historical exact uvLMM and corrected-panel added-marker uvLMM evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import HISTORICAL_ADDITIVE_GALWEY_MEFF, HISTORICAL_UV_THRESHOLD, add_galwey_adjusted_pvalues, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--new-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--threshold", default=HISTORICAL_UV_THRESHOLD, type=float)
    parser.add_argument("--galwey-meff", default=HISTORICAL_ADDITIVE_GALWEY_MEFF, type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    old_p = read_csv(args.legacy_dir / "uniLMM_p_val_normal_MIC_pneumo.csv")
    old_e = read_csv(args.legacy_dir / "uniLMM_effect_normal_MIC_pneumo.csv")
    old_cand = old_e[old_e["effect_type"].eq("candidate")][["trait", "effect_name", "effsize", "effsize_se"]].reset_index(
        drop=True
    )
    if len(old_cand) != len(old_p):
        raise ValueError("Historical uvLMM p-value/effect rows do not align")
    old = old_cand.rename(columns={"trait": "drug", "effect_name": "marker"})
    old["pv20"] = old_p["pv20"].to_numpy()
    old["source"] = "historical_157"

    added_path = args.new_dir / "uniLMM_exact_added_markers_p_values.csv"
    added = read_csv(added_path)[["marker", "drug", "pv20"]].copy()
    added_effect_path = args.new_dir / "uniLMM_exact_added_markers_effect_sizes.csv"
    if added_effect_path.exists():
        added_effect = read_csv(added_effect_path)
        added_effect = added_effect[added_effect["effect_type"].eq("candidate")][
            ["marker", "drug", "effsize", "effsize_se"]
        ].copy()
        added = added.merge(added_effect, on=["marker", "drug"], how="left", validate="one_to_one")
    else:
        added["effsize"] = pd.NA
        added["effsize_se"] = pd.NA
    added["source"] = "rerun_added_13"

    combined = pd.concat(
        [old[["marker", "drug", "pv20", "effsize", "effsize_se", "source"]], added],
        ignore_index=True,
    )
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
                "uvLMM_exact_source": ";".join(sorted(group["source"].unique())),
                "uvLMM_exact_threshold": args.threshold,
                "uvLMM_exact_galwey_meff": args.galwey_meff,
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
