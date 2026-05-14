#!/usr/bin/env python3
"""Merge corrected-panel epistasis chunks and derive permutation threshold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--alpha", default=0.05, type=float)
    return parser.parse_args()


def concat(pattern: str) -> pd.DataFrame:
    files = sorted(Path().glob(pattern))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(path) for path in files], ignore_index=True)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    observed_files = sorted(
        path
        for path in args.chunk_dir.glob("epistasis_p_values.chunk_*.csv")
        if "smoke" not in path.name
    )
    effect_files = sorted(
        path
        for path in args.chunk_dir.glob("epistasis_effect_sizes.chunk_*.csv")
        if "smoke" not in path.name
    )
    perm_files = sorted(
        path
        for path in args.chunk_dir.glob("epistasis_perm_p_values.chunk_*.csv")
        if "perm_9999_" not in path.name and "smoke" not in path.name
    )
    if not observed_files:
        raise FileNotFoundError(f"No observed chunks found in {args.chunk_dir}")

    observed = pd.concat([pd.read_csv(path) for path in observed_files], ignore_index=True)
    effects = pd.concat([pd.read_csv(path) for path in effect_files], ignore_index=True) if effect_files else pd.DataFrame()
    permutations = pd.concat([pd.read_csv(path) for path in perm_files], ignore_index=True) if perm_files else pd.DataFrame()

    threshold = np.nan
    min_by_perm = pd.DataFrame()
    if not permutations.empty:
        min_by_perm = permutations.groupby("permutation_seed", dropna=False)["pv20"].min().reset_index()
        threshold = float(min_by_perm["pv20"].quantile(args.alpha, interpolation="higher"))
        observed["significant_epistasis_perm"] = observed["pv20"] <= threshold
    else:
        observed["significant_epistasis_perm"] = False

    observed.to_csv(args.out_dir / "corrected_epistasis_p_values.csv", index=False)
    effects.to_csv(args.out_dir / "corrected_epistasis_effect_sizes.csv", index=False)
    permutations.to_csv(args.out_dir / "corrected_epistasis_permutation_p_values.csv", index=False)
    min_by_perm.to_csv(args.out_dir / "corrected_epistasis_permutation_minima.csv", index=False)

    sig = observed[observed["significant_epistasis_perm"]].copy()
    rows: list[dict[str, object]] = []
    for marker in sorted(set(sig["parent_a"]).union(set(sig["parent_b"]))):
        subset = sig[(sig["parent_a"] == marker) | (sig["parent_b"] == marker)]
        rows.append(
            {
                "marker": marker,
                "num_sig_interactions": int(len(subset)),
                "min_epistasis_pv20": float(subset["pv20"].min()) if len(subset) else np.nan,
            }
        )
    support = pd.DataFrame(rows)
    support.to_csv(args.out_dir / "corrected_epistasis_marker_support.csv", index=False)

    observed_non_ok = (
        int((observed["fit_status"].notna() & (observed["fit_status"] != "ok")).sum())
        if "fit_status" in observed
        else 0
    )
    permutation_non_ok = (
        int((permutations["fit_status"].notna() & (permutations["fit_status"] != "ok")).sum())
        if "fit_status" in permutations
        else 0
    )

    summary = {
        "observed_interactions": int(len(observed)),
        "effect_rows": int(len(effects)),
        "permutation_rows": int(len(permutations)),
        "permutations": int(permutations["permutation_seed"].nunique()) if not permutations.empty else 0,
        "observed_non_ok_fit_rows": observed_non_ok,
        "permutation_non_ok_fit_rows": permutation_non_ok,
        "alpha": args.alpha,
        "permutation_threshold": threshold,
        "significant_interactions": int(observed["significant_epistasis_perm"].sum()),
        "markers_with_epistasis_support": int(len(support)),
    }
    (args.out_dir / "corrected_epistasis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
