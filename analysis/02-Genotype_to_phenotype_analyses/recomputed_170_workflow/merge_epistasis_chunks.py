#!/usr/bin/env python3
"""Merge epistasis chunks and derive marker support."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    EXPECTED_CORRECTED_EPISTASIS_CANDIDATES,
    HISTORICAL_EPISTASIS_EFFECT_LOWER_BOUND,
    HISTORICAL_EPISTASIS_GALWEY_MEFF,
    HISTORICAL_EPISTASIS_THRESHOLD,
    add_galwey_adjusted_pvalues,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--alpha", default=0.05, type=float)
    parser.add_argument("--galwey-meff", default=HISTORICAL_EPISTASIS_GALWEY_MEFF, type=float)
    parser.add_argument("--threshold", default=HISTORICAL_EPISTASIS_THRESHOLD, type=float)
    parser.add_argument(
        "--threshold-mode",
        choices=["historical", "lowest-min-p", "explicit"],
        default="historical",
        help="How to choose the adjusted-p support threshold.",
    )
    parser.add_argument("--support-effect-lower-bound", default=HISTORICAL_EPISTASIS_EFFECT_LOWER_BOUND, type=float)
    parser.add_argument("--expected-observed-interactions", default=EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, type=int)
    parser.add_argument("--expected-permutation-rows", type=int)
    parser.add_argument("--expected-permutations", type=int)
    parser.add_argument("--fail-on-non-ok", action="store_true")
    return parser.parse_args()


def read_many(files: list[Path]) -> pd.DataFrame:
    return pd.concat([pd.read_csv(path) for path in files], ignore_index=True) if files else pd.DataFrame()


def candidate_effect_wide(effects: pd.DataFrame) -> pd.DataFrame:
    if effects.empty:
        raise ValueError("Observed effect chunks are required for historical epistasis support")
    candidate = effects[effects["effect_type"].eq("candidate")].copy()
    required = {"test", "effect_name", "env", "effsize", "effsize_se"}
    missing = required - set(candidate.columns)
    if missing:
        raise ValueError(f"Effect chunks are missing columns: {sorted(missing)}")
    wide = candidate.pivot(index=["test", "effect_name"], columns="env", values=["effsize", "effsize_se"])
    wide.columns = [f"{value}_{env}" for value, env in wide.columns]
    wide = wide.reset_index().rename(columns={"effect_name": "interaction"})
    expected = ["effsize_env1_D1", "effsize_env1_D2", "effsize_se_env1_D1", "effsize_se_env1_D2"]
    missing_effects = [column for column in expected if column not in wide.columns]
    if missing_effects:
        raise ValueError(f"Effect chunks are missing candidate env columns: {missing_effects}")
    return wide


def enrich_observed(observed: pd.DataFrame, effects: pd.DataFrame, *, galwey_meff: float, threshold: float, lower_bound: float) -> pd.DataFrame:
    required = {"test", "interaction", "parent_a", "parent_b", "pv20"}
    missing = required - set(observed.columns)
    if missing:
        raise ValueError(f"Observed chunks are missing columns: {sorted(missing)}")

    effects_wide = candidate_effect_wide(effects)
    out = observed.merge(effects_wide, on=["test", "interaction"], how="left", validate="one_to_one")
    if out[["effsize_env1_D1", "effsize_env1_D2", "effsize_se_env1_D1", "effsize_se_env1_D2"]].isna().any().any():
        raise ValueError("Some observed p-value rows did not receive candidate effect sizes")

    out = add_galwey_adjusted_pvalues(out, meff=galwey_meff)
    out["joint_effect_size"] = np.sqrt(out["effsize_env1_D1"] ** 2 + out["effsize_env1_D2"] ** 2)
    out["joint_effect_size_se"] = np.sqrt(out["effsize_se_env1_D1"] ** 2 + out["effsize_se_env1_D2"] ** 2)
    out["joint_effect_lower_bound"] = out["joint_effect_size"] - out["joint_effect_size_se"]
    out["passes_epistasis_threshold"] = out["pv20_adj_galwey"] <= threshold
    out["passes_epistasis_effect_filter"] = out["joint_effect_lower_bound"] >= lower_bound
    out["epistasis_support"] = out["passes_epistasis_threshold"] & out["passes_epistasis_effect_filter"]
    return out


def permutation_minima(permutations: pd.DataFrame, *, galwey_meff: float, alpha: float) -> tuple[pd.DataFrame, dict[str, object]]:
    if permutations.empty:
        return pd.DataFrame(), {
            "permutations": 0,
            "raw_permutation_threshold": np.nan,
            "adjusted_permutation_threshold": np.nan,
        }
    adjusted = add_galwey_adjusted_pvalues(permutations, meff=galwey_meff)
    minima = (
        adjusted.groupby("permutation_seed", dropna=False)
        .agg(raw_min_pv20=("pv20", "min"), adjusted_min_pv20=("pv20_adj_galwey", "min"))
        .reset_index()
        .sort_values("raw_min_pv20", kind="mergesort")
    )
    lowest = minima.iloc[0]
    summary = {
        "permutations": int(minima["permutation_seed"].nunique()),
        "raw_lowest_min_p_threshold": float(lowest["raw_min_pv20"]),
        "adjusted_lowest_min_p_threshold": float(lowest["adjusted_min_pv20"]),
        "lowest_min_p_permutation_seed": int(lowest["permutation_seed"]),
        "raw_permutation_threshold": float(minima["raw_min_pv20"].quantile(alpha, interpolation="higher")),
        "adjusted_permutation_threshold": float(
            minima["adjusted_min_pv20"].quantile(alpha, interpolation="higher")
        ),
    }
    return minima, summary


def marker_support(enriched: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    supported = enriched[enriched["epistasis_support"]].copy()
    p_only = enriched[enriched["passes_epistasis_threshold"]].copy()

    support_parents = pd.concat(
        [
            supported.assign(marker=supported["parent_a"], partner=supported["parent_b"]),
            supported.assign(marker=supported["parent_b"], partner=supported["parent_a"]),
        ],
        ignore_index=True,
    )
    p_only_parents = pd.concat(
        [
            p_only.assign(marker=p_only["parent_a"]),
            p_only.assign(marker=p_only["parent_b"]),
        ],
        ignore_index=True,
    )

    if support_parents.empty:
        support = pd.DataFrame(
            columns=[
                "marker",
                "num_sig_interactions",
                "num_perm_cutoff_interactions",
                "min_epistasis_pv20",
                "min_epistasis_pv20_adj_galwey",
                "min_joint_effect_lower_bound",
                "max_joint_effect_size",
            ]
        )
    else:
        support = (
            support_parents.groupby("marker", dropna=False)
            .agg(
                num_sig_interactions=("interaction", "count"),
                min_epistasis_pv20=("pv20", "min"),
                min_epistasis_pv20_adj_galwey=("pv20_adj_galwey", "min"),
                min_joint_effect_lower_bound=("joint_effect_lower_bound", "min"),
                max_joint_effect_size=("joint_effect_size", "max"),
            )
            .reset_index()
        )

    if p_only_parents.empty:
        p_only_counts = pd.DataFrame(columns=["marker", "num_perm_cutoff_interactions"])
    else:
        p_only_counts = (
            p_only_parents.groupby("marker", dropna=False)
            .agg(num_perm_cutoff_interactions=("interaction", "count"))
            .reset_index()
        )
    support = support.drop(columns=["num_perm_cutoff_interactions"], errors="ignore").merge(
        p_only_counts, on="marker", how="left"
    )
    support["num_perm_cutoff_interactions"] = support["num_perm_cutoff_interactions"].fillna(0).astype(int)
    return supported, support.sort_values(["num_sig_interactions", "marker"], ascending=[False, True])


def non_ok_count(df: pd.DataFrame) -> int:
    if df.empty or "fit_status" not in df.columns:
        return 0
    return int((df["fit_status"].notna() & (df["fit_status"] != "ok")).sum())


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    observed_files = sorted(
        path for path in args.chunk_dir.glob("epistasis_p_values.chunk_*.csv") if "smoke" not in path.name
    )
    effect_files = sorted(
        path for path in args.chunk_dir.glob("epistasis_effect_sizes.chunk_*.csv") if "smoke" not in path.name
    )
    perm_files = sorted(
        path
        for path in args.chunk_dir.glob("epistasis_perm_p_values.chunk_*.csv")
        if "perm_9999_" not in path.name and "smoke" not in path.name
    )
    if not observed_files:
        raise FileNotFoundError(f"No observed chunks found in {args.chunk_dir}")

    observed_raw = read_many(observed_files)
    if args.expected_observed_interactions is not None and len(observed_raw) != args.expected_observed_interactions:
        raise AssertionError(
            f"Expected {args.expected_observed_interactions} observed interactions, found {len(observed_raw)}"
        )
    effects = read_many(effect_files)
    permutations = read_many(perm_files)
    if args.expected_permutation_rows is not None and len(permutations) != args.expected_permutation_rows:
        raise AssertionError(f"Expected {args.expected_permutation_rows} permutation rows, found {len(permutations)}")
    if args.expected_permutations is not None:
        actual_permutations = int(permutations["permutation_seed"].nunique()) if not permutations.empty else 0
        if actual_permutations != args.expected_permutations:
            raise AssertionError(f"Expected {args.expected_permutations} permutations, found {actual_permutations}")
    if args.fail_on_non_ok and (non_ok_count(observed_raw) or non_ok_count(permutations)):
        raise AssertionError("Primary epistasis merge contains non-ok fit rows")

    minima, permutation_summary = permutation_minima(permutations, galwey_meff=args.galwey_meff, alpha=args.alpha)
    if args.threshold_mode == "lowest-min-p":
        support_threshold = float(permutation_summary["adjusted_lowest_min_p_threshold"])
        threshold_policy = "lowest permutation minimum"
    elif args.threshold_mode == "explicit":
        support_threshold = float(args.threshold)
        threshold_policy = "explicit adjusted threshold"
    else:
        support_threshold = float(args.threshold)
        threshold_policy = "historical literal adjusted threshold"

    observed = enrich_observed(
        observed_raw,
        effects,
        galwey_meff=args.galwey_meff,
        threshold=support_threshold,
        lower_bound=args.support_effect_lower_bound,
    )
    supported, support = marker_support(observed)

    observed.to_csv(args.out_dir / "corrected_epistasis_p_values.csv", index=False)
    effects.to_csv(args.out_dir / "corrected_epistasis_effect_sizes.csv", index=False)
    permutations.to_csv(args.out_dir / "corrected_epistasis_permutation_p_values.csv", index=False)
    minima.to_csv(args.out_dir / "corrected_epistasis_permutation_minima.csv", index=False)
    supported.to_csv(args.out_dir / "corrected_epistasis_supported_interactions.csv", index=False)
    support.to_csv(args.out_dir / "corrected_epistasis_marker_support.csv", index=False)

    summary = {
        "support_rule": "pv20_adj_galwey <= threshold AND joint_effect_size - joint_effect_size_se >= lower_bound",
        "galwey_meff": float(args.galwey_meff),
        "threshold_mode": args.threshold_mode,
        "threshold_policy": threshold_policy,
        "epistasis_threshold": float(support_threshold),
        "historical_epistasis_threshold": float(args.threshold),
        "support_effect_lower_bound": float(args.support_effect_lower_bound),
        "observed_interactions": int(len(observed)),
        "effect_rows": int(len(effects)),
        "permutation_rows": int(len(permutations)),
        **permutation_summary,
        "observed_non_ok_fit_rows": non_ok_count(observed),
        "permutation_non_ok_fit_rows": non_ok_count(permutations),
        "alpha": float(args.alpha),
        "pvalue_threshold_only_interactions": int(observed["passes_epistasis_threshold"].sum()),
        "support_interactions": int(observed["epistasis_support"].sum()),
        "markers_with_epistasis_support": int(len(support)),
    }
    (args.out_dir / "corrected_epistasis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
