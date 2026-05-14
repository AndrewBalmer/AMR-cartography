#!/usr/bin/env python3
"""Run a chunk of corrected-panel epistatic mvLMM interaction tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import chunk_bounds, lowrank_multitrait_scan, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--interaction-file", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--array-index", type=int)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--chunk-label", default=None)
    parser.add_argument("--permutation-seed", type=int, default=None)
    parser.add_argument("--pvalues-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_bounds(args: argparse.Namespace, total: int) -> tuple[int, int]:
    if args.array_index is not None:
        return chunk_bounds(total, args.chunk_size, args.array_index)
    if args.start is None or args.end is None:
        raise ValueError("Provide either --array-index/--chunk-size or --start/--end")
    return args.start, min(args.end, total)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    interactions = read_csv(args.interaction_file)
    start, end = resolve_bounds(args, len(interactions))
    label = args.chunk_label or f"{start:06d}_{end:06d}"
    prefix = "epistasis_perm" if args.permutation_seed is not None else "epistasis"
    pvalue_path = args.out_dir / f"{prefix}_p_values.chunk_{label}.csv"
    effect_path = args.out_dir / f"{prefix}_effect_sizes.chunk_{label}.csv"
    if pvalue_path.exists() and not args.force:
        print(f"skipping existing {pvalue_path}")
        return

    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")
    map_coords = read_csv(args.data_dir / "S.pneumo_map_mvlmm_map_coords.csv")
    trait_design = read_csv(args.data_dir / "S.pneumo_map_mvlmm_mod_matrix.csv").to_numpy(dtype=float)
    if args.permutation_seed is not None:
        map_coords = map_coords.sample(frac=1, random_state=args.permutation_seed).reset_index(drop=True)

    stats_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    y = map_coords[["D1", "D2"]].to_numpy(dtype=float)

    for local_index, row in interactions.iloc[start:end].reset_index(drop=True).iterrows():
        interaction = str(row["interaction"])
        parent_a = str(row["parent_a"])
        parent_b = str(row["parent_b"])
        print(f"{start + local_index}: {interaction}", flush=True)

        candidate_values = (
            markers[parent_a].to_numpy(dtype=float) * markers[parent_b].to_numpy(dtype=float)
        ).reshape(-1, 1)
        covariates = relatedness[[parent_a, parent_b]].to_numpy(dtype=float)
        dropped = relatedness.drop(columns=[parent_a, parent_b]).to_numpy(dtype=float)
        result = lowrank_multitrait_scan(
            y=y,
            trait_design=trait_design,
            covariates=covariates,
            relatedness_features=dropped,
            candidate=candidate_values,
        )

        test_id = start + local_index
        stats_rows.append(
            {
                "test": test_id,
                "interaction": interaction,
                "parent_a": parent_a,
                "parent_b": parent_b,
                "n_present": int(row["n_present"]),
                "permutation_seed": args.permutation_seed,
                "lml0": result["lml0"],
                "lml2": result["lml2"],
                "dof20": result["dof20"],
                "scale2": result["scale2"],
                "pv20": result["pv20"],
            }
        )
        if not args.pvalues_only:
            for env, eff, se in zip(["env1_D1", "env1_D2"], result["candidate_eff"], result["candidate_se"]):
                effect_rows.append(
                    {
                        "test": test_id,
                        "trait": "map",
                        "effect_type": "candidate",
                        "effect_name": interaction,
                        "env": env,
                        "effsize": eff,
                        "effsize_se": se,
                    }
                )

    pd.DataFrame(stats_rows).to_csv(pvalue_path, index=False)
    if not args.pvalues_only:
        pd.DataFrame(effect_rows).to_csv(effect_path, index=False)
    print("wrote", pvalue_path, len(stats_rows))
    if not args.pvalues_only:
        print("wrote", effect_path, len(effect_rows))


if __name__ == "__main__":
    main()
