#!/usr/bin/env python3
"""Compare old manuscript and corrected additive mvLMM outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import DEFAULT_NEW_THRESHOLD, marker_order_from_effects, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-dir", required=True, type=Path)
    parser.add_argument("--new-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--old-threshold", default=0.000588, type=float)
    parser.add_argument("--new-threshold", default=DEFAULT_NEW_THRESHOLD, type=float)
    return parser.parse_args()


def load_additive(directory: Path) -> pd.DataFrame:
    pval_path = directory / "mvLMM_p_values_normal_pneumo_low_freq_vars.csv"
    effect_path = directory / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv"
    pvals = read_csv(pval_path).copy()
    unnamed = [col for col in pvals.columns if str(col).startswith("Unnamed")]
    if unnamed:
        pvals = pvals.rename(columns={unnamed[0]: "result_index"})
    order = marker_order_from_effects(effect_path)
    if len(order) != len(pvals):
        raise ValueError(f"{directory}: {len(order)} marker names for {len(pvals)} p-value rows")
    pvals["marker"] = order

    effects = read_csv(effect_path)
    cand = effects[effects["effect_type"].eq("candidate")][["effect_name", "env", "effsize"]]
    wide = cand.pivot(index="effect_name", columns="env", values="effsize").reset_index()
    wide = wide.rename(columns={"effect_name": "marker", "env1_D1": "effect_axis1", "env1_D2": "effect_axis2"})
    out = pvals.merge(wide, on="marker", how="left")
    out["joint_effect_size"] = np.sqrt(out["effect_axis1"] ** 2 + out["effect_axis2"] ** 2)
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    old = load_additive(args.old_dir).add_prefix("old_").rename(columns={"old_marker": "marker"})
    new = load_additive(args.new_dir).add_prefix("new_").rename(columns={"new_marker": "marker"})
    merged = old.merge(new, on="marker", how="outer", indicator=True)

    merged["old_sig_old_threshold"] = merged["old_pv20"] < args.old_threshold
    merged["old_sig_new_threshold"] = merged["old_pv20"] < args.new_threshold
    merged["new_sig_new_threshold"] = merged["new_pv20"] < args.new_threshold
    merged["panel_status"] = merged["_merge"].map(
        {"both": "shared", "left_only": "dropped_from_corrected", "right_only": "added_in_corrected"}
    )

    shared = merged[merged["panel_status"].eq("shared")].copy()
    shared["significance_change_new_threshold"] = np.select(
        [
            shared["old_sig_new_threshold"] & shared["new_sig_new_threshold"],
            ~shared["old_sig_new_threshold"] & shared["new_sig_new_threshold"],
            shared["old_sig_new_threshold"] & ~shared["new_sig_new_threshold"],
        ],
        ["shared_significant", "gained", "lost"],
        default="shared_not_significant",
    )

    summary = {
        "old_markers": int(merged["old_pv20"].notna().sum()),
        "new_markers": int(merged["new_pv20"].notna().sum()),
        "shared_markers": int(merged["panel_status"].eq("shared").sum()),
        "added_markers": int(merged["panel_status"].eq("added_in_corrected").sum()),
        "dropped_markers": int(merged["panel_status"].eq("dropped_from_corrected").sum()),
        "old_significant_old_threshold": int((merged["old_pv20"] < args.old_threshold).sum()),
        "old_significant_new_threshold": int((merged["old_pv20"] < args.new_threshold).sum()),
        "new_significant_new_threshold": int((merged["new_pv20"] < args.new_threshold).sum()),
        "added_significant_new_threshold": int(
            (merged["panel_status"].eq("added_in_corrected") & merged["new_sig_new_threshold"]).sum()
        ),
        "shared_status_counts": shared["significance_change_new_threshold"].value_counts().to_dict(),
        "shared_log10p_pearson": float(
            np.corrcoef(-np.log10(shared["old_pv20"]), -np.log10(shared["new_pv20"]))[0, 1]
        ),
    }

    merged.to_csv(args.out_dir / "additive_old_vs_new_marker_comparison.csv", index=False)
    shared.to_csv(args.out_dir / "additive_shared_marker_comparison.csv", index=False)
    merged.loc[merged["panel_status"].eq("added_in_corrected"), ["marker"]].to_csv(
        args.out_dir / "added_markers.csv", index=False
    )
    (args.out_dir / "additive_old_vs_new_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
