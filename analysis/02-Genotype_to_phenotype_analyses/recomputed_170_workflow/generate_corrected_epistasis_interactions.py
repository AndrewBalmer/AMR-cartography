#!/usr/bin/env python3
"""Generate epistasis candidates with the original manuscript four-cell rule.

Input: the corrected 170-marker binary matrix.
Output: interaction metadata for marker pairs where all four parent genotype
cells are present in at least 1% of isolates. This script defines the epistasis
testing universe, so review it before reviewing any epistasis model results.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from common import (
    EXPECTED_CORRECTED_ADDITIVE_MARKERS,
    EXPECTED_CORRECTED_EPISTASIS_CANDIDATES,
    binary_marker_matrix,
    four_cell_interaction_metadata,
    read_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-fraction", default=0.01, type=float)
    parser.add_argument("--expected-marker-count", default=EXPECTED_CORRECTED_ADDITIVE_MARKERS, type=int)
    parser.add_argument("--expected-interaction-count", default=EXPECTED_CORRECTED_EPISTASIS_CANDIDATES, type=int)
    parser.add_argument(
        "--write-matrix",
        action="store_true",
        help="Also write the full interaction marker matrix for direct inspection.",
    )
    return parser.parse_args()


def write_interaction_matrix(markers: pd.DataFrame, interactions: pd.DataFrame, path: Path) -> None:
    matrix = binary_marker_matrix(markers)
    out = {
        row.interaction: (
            matrix[str(row.parent_a)].to_numpy(dtype="uint8")
            * matrix[str(row.parent_b)].to_numpy(dtype="uint8")
        )
        for row in interactions.itertuples(index=False)
    }
    pd.DataFrame(out).to_csv(path, index=False)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv",)
    missing = [marker for marker in markers.columns if marker not in relatedness.columns]
    if missing:
        raise ValueError(f"Test markers missing from relatedness matrix: {missing[:10]}")
    if args.expected_marker_count is not None and markers.shape[1] != args.expected_marker_count:
        raise AssertionError(f"Expected {args.expected_marker_count} markers, found {markers.shape[1]}")

    interactions = four_cell_interaction_metadata(markers, min_fraction=args.min_fraction)
    if args.expected_interaction_count is not None and len(interactions) != args.expected_interaction_count:
        raise AssertionError(
            f"Expected {args.expected_interaction_count} four-cell interactions, found {len(interactions)}"
        )

    interactions.to_csv(args.out_dir / "corrected_epistasis_interactions.csv", index=False)
    if args.write_matrix:
        write_interaction_matrix(markers, interactions, args.out_dir / "corrected_epistasis_interaction_matrix.csv")

    possible = math.comb(markers.shape[1], 2)
    summary = {
        "generation_rule": "original_model_matrix_pair_products_four_genotype_cells_ge_1_percent",
        "hash_deduplicated": False,
        "additive_markers": int(markers.shape[1]),
        "isolates": int(markers.shape[0]),
        "possible_pairwise_interactions": int(possible),
        "valid_four_cell_interactions": int(len(interactions)),
        "failed_four_cell_filter": int(possible - len(interactions)),
        "min_fraction": float(args.min_fraction),
        "one_percent_count": float(markers.shape[0] * args.min_fraction),
        "min_cell_count_min": int(interactions["min_cell_count"].min()) if len(interactions) else None,
        "min_cell_count_max": int(interactions["min_cell_count"].max()) if len(interactions) else None,
    }
    (args.out_dir / "corrected_epistasis_interaction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame([{"metric": key, "value": value} for key, value in summary.items()]).to_csv(
        args.out_dir / "corrected_epistasis_interaction_summary.csv", index=False
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
