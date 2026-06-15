#!/usr/bin/env python3
"""Run a chunk of corrected-panel additive mvLMM marker tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from xarray import DataArray

from common import chunk_bounds, install_pandas_limix_shim, read_csv


install_pandas_limix_shim()

from limix.qtl import scan  # noqa: E402
from limix.stats import linear_kinship  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--array-index", type=int)
    parser.add_argument("--chunk-size", type=int, default=5)
    parser.add_argument("--chunk-label")
    parser.add_argument("--permutation-index", type=int)
    parser.add_argument("--base-seed", default=2, type=int)
    parser.add_argument("--pvalues-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_bounds(args: argparse.Namespace, total: int) -> tuple[int, int]:
    if args.array_index is not None:
        return chunk_bounds(total, args.chunk_size, args.array_index)
    if args.start is None or args.end is None:
        raise ValueError("Provide either --array-index/--chunk-size or --start/--end")
    return args.start, min(args.end, total)


def stats_frame_to_row(stats: pd.DataFrame, *, marker: str, marker_index: int, test_id: int) -> dict[str, object]:
    if len(stats) != 1:
        raise ValueError(f"Expected one stats row for {marker}, got {len(stats)}")
    row = stats.reset_index(drop=True).iloc[0].to_dict()
    row["test"] = test_id
    row["marker"] = marker
    row["marker_index"] = marker_index
    row["fit_status"] = "ok"
    row["fit_error"] = ""
    return row


def effect_rows(effects: pd.DataFrame, *, marker: str, marker_index: int, test_id: int) -> list[dict[str, object]]:
    out = effects.reset_index(drop=True).copy()
    out["test"] = test_id
    out["marker"] = marker
    out["marker_index"] = marker_index
    return out.to_dict(orient="records")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")
    map_coords = read_csv(args.data_dir / "S.pneumo_map_mvlmm_map_coords.csv")
    if not (len(markers) == len(relatedness) == len(map_coords)):
        raise ValueError("Marker, relatedness, and phenotype row counts do not match")

    if args.permutation_index is not None:
        seed = args.base_seed + args.permutation_index
        map_coords = map_coords.sample(frac=1, random_state=seed).reset_index(drop=True)
    else:
        seed = None

    start, end = resolve_bounds(args, markers.shape[1])
    label = args.chunk_label or f"{start:03d}_{end:03d}"
    prefix = "additive_perm" if args.permutation_index is not None else "additive"
    pvalue_path = args.out_dir / f"{prefix}_p_values.chunk_{label}.csv"
    effect_path = args.out_dir / f"{prefix}_effect_sizes.chunk_{label}.csv"
    if pvalue_path.exists() and not args.force:
        print(f"skipping existing {pvalue_path}")
        return

    stats_rows: list[dict[str, object]] = []
    effects_out: list[dict[str, object]] = []
    trait_labels = map_coords.columns.astype(object).to_numpy()
    phenotype_cov = map_coords.cov().to_numpy(dtype=float)
    null_trait_design = DataArray(
        np.empty((len(trait_labels), 0)),
        dims=["sample", "env"],
        coords={"env": np.asarray([], dtype=object)},
    )
    candidate_trait_design = DataArray(
        np.eye(len(trait_labels)),
        dims=["sample", "env"],
        coords={"env": trait_labels},
    )

    for marker_index in range(start, end):
        marker = str(markers.columns[marker_index])
        print(f"{marker_index}: {marker}", flush=True)
        test_marker = markers[[marker]]
        kinship_features = relatedness.drop(columns=[marker]).to_numpy(dtype=float)
        kinship = linear_kinship(kinship_features)
        result = scan(
            G=test_marker,
            Y=map_coords,
            K=kinship,
            A=phenotype_cov,
            A0=null_trait_design,
            A1=candidate_trait_design,
            verbose=False,
        )
        test_id = marker_index
        stats_row = stats_frame_to_row(result.stats, marker=marker, marker_index=marker_index, test_id=test_id)
        if args.permutation_index is not None:
            stats_row["repeat_p_index"] = args.permutation_index
            stats_row["permutation_seed"] = seed
        stats_rows.append(stats_row)

        if not args.pvalues_only:
            rows = effect_rows(result.effsizes["h2"], marker=marker, marker_index=marker_index, test_id=test_id)
            if args.permutation_index is not None:
                for row in rows:
                    row["repeat_eff_index"] = args.permutation_index
                    row["permutation_seed"] = seed
            effects_out.extend(rows)

    pd.DataFrame(stats_rows).to_csv(pvalue_path, index=False)
    print("wrote", pvalue_path, len(stats_rows))
    if not args.pvalues_only:
        pd.DataFrame(effects_out).to_csv(effect_path, index=False)
        print("wrote", effect_path, len(effects_out))


if __name__ == "__main__":
    main()
