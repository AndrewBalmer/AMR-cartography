#!/usr/bin/env python3
"""Run exact leave-one-marker-out univariate LMM chunks for the corrected panel.

The uvLMM stream is a display/comparison analysis, not an evidence stream. This
runner fits 170 markers x 6 drugs in chunks while excluding the tested marker
from the relatedness features for that test.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs_linear
from scipy.stats import chi2

from common import chunk_bounds, normalized_relatedness_features, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--array-index", type=int)
    parser.add_argument("--chunk-size", type=int, default=5)
    parser.add_argument("--chunk-label")
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

    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")
    mic_values = read_csv(args.legacy_dir / "S.pneumo_map_mvlmm_MIC_values.csv")
    if not (len(markers) == len(relatedness) == len(mic_values)):
        raise ValueError("Marker, relatedness, and MIC row counts do not match")

    start, end = resolve_bounds(args, markers.shape[1])
    label = args.chunk_label or f"{start:03d}_{end:03d}"
    pvalue_path = args.out_dir / f"uniLMM_exact_p_values.chunk_{label}.csv"
    effect_path = args.out_dir / f"uniLMM_exact_effect_sizes.chunk_{label}.csv"
    if pvalue_path.exists() and not args.force:
        print(f"skipping existing {pvalue_path}")
        return

    stats_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    test_index = start * len(mic_values.columns)

    for marker_index in range(start, end):
        marker = str(markers.columns[marker_index])
        print("marker", marker_index, marker, flush=True)
        qs = economic_qs_linear(normalized_relatedness_features(relatedness.drop(columns=[marker]).to_numpy()))
        candidate = markers[[marker]].to_numpy(dtype=float)
        for drug in mic_values.columns:
            print("  drug", drug, flush=True)
            y = mic_values[drug].to_numpy(dtype=float)
            lmm = LMM(y, np.ones((len(y), 1)), qs, restricted=False)
            lmm.fit(verbose=False)
            result = lmm.get_fast_scanner().scan(candidate)

            lml0 = float(lmm.lml())
            lml2 = float(result["lml"])
            stats_rows.append(
                {
                    "test": test_index,
                    "lml0": lml0,
                    "lml2": lml2,
                    "dof20": 1,
                    "scale2": float(result["scale"]),
                    "pv20": float(chi2.sf(max(0.0, 2.0 * (lml2 - lml0)), 1)),
                    "marker": marker,
                    "marker_index": marker_index,
                    "drug": drug,
                    "fit_status": "ok",
                    "fit_error": "",
                }
            )
            effect_rows.extend(
                [
                    {
                        "test": test_index,
                        "trait": drug,
                        "effect_type": "covariate",
                        "effect_name": "offset",
                        "effsize": float(result["effsizes0"][0]),
                        "effsize_se": float(result["effsizes0_se"][0]),
                        "marker": marker,
                        "marker_index": marker_index,
                        "drug": drug,
                    },
                    {
                        "test": test_index,
                        "trait": drug,
                        "effect_type": "candidate",
                        "effect_name": marker,
                        "effsize": float(result["effsizes1"][0]),
                        "effsize_se": float(result["effsizes1_se"][0]),
                        "marker": marker,
                        "marker_index": marker_index,
                        "drug": drug,
                    },
                ]
            )
            test_index += 1

    pd.DataFrame(stats_rows).to_csv(pvalue_path, index=False)
    pd.DataFrame(effect_rows).to_csv(effect_path, index=False)
    print("wrote", pvalue_path, len(stats_rows))
    print("wrote", effect_path, len(effect_rows))


if __name__ == "__main__":
    main()
