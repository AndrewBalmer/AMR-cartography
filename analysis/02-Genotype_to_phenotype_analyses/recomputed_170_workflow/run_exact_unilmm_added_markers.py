#!/usr/bin/env python3
"""Exact leave-one-marker-out univariate LMM for corrected-panel added markers."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs_linear
from scipy.stats import chi2

from common import normalized_relatedness_features, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--added-marker-file", required=True, type=Path)
    parser.add_argument("--out-prefix", default="uniLMM_exact_added_markers")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")
    mic_values = read_csv(args.legacy_dir / "S.pneumo_map_mvlmm_MIC_values.csv")
    added = read_csv(args.added_marker_file)["marker"].astype(str).tolist()

    missing = [marker for marker in added if marker not in markers.columns or marker not in relatedness.columns]
    if missing:
        raise ValueError(f"Added markers missing from marker/relatedness matrices: {missing}")
    if not (len(markers) == len(relatedness) == len(mic_values)):
        raise ValueError("Marker, relatedness, and MIC row counts do not match")

    stats_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    test_index = 0

    for marker in added:
        print("marker", marker, flush=True)
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
                    "drug": drug,
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
                        "drug": drug,
                    },
                ]
            )
            test_index += 1

    pvals = pd.DataFrame(stats_rows)
    effects = pd.DataFrame(effect_rows)
    pvals.to_csv(args.data_dir / f"{args.out_prefix}_p_values.csv", index=False)
    effects.to_csv(args.data_dir / f"{args.out_prefix}_effect_sizes.csv", index=False)
    print("wrote", pvals.shape, effects.shape)


if __name__ == "__main__":
    main()
