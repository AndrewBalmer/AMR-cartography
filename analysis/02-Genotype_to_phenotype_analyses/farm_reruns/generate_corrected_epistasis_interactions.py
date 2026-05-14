#!/usr/bin/env python3
"""Generate valid pairwise interaction metadata for the corrected additive panel."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--write-matrix", action="store_true", help="Also write the full interaction marker matrix")
    return parser.parse_args()


def marker_hash(values: np.ndarray) -> str:
    packed = np.packbits(values.astype(np.uint8))
    return hashlib.sha1(packed.tobytes()).hexdigest()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")

    missing = [marker for marker in markers.columns if marker not in relatedness.columns]
    if missing:
        raise ValueError(f"Test markers missing from relatedness matrix: {missing[:10]}")

    seen: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    matrix: dict[str, np.ndarray] = {}
    counts = {"valid": 0, "constant": 0, "same_parent": 0, "duplicate": 0}
    names = list(markers.columns)

    for i, parent_a in enumerate(names):
        a = markers[parent_a].to_numpy(dtype=np.uint8)
        for parent_b in names[i + 1 :]:
            b = markers[parent_b].to_numpy(dtype=np.uint8)
            interaction = (a * b).astype(np.uint8)
            interaction_name = f"{parent_a}:{parent_b}"

            if interaction.min() == interaction.max():
                counts["constant"] += 1
                continue
            if np.array_equal(interaction, a) or np.array_equal(interaction, b):
                counts["same_parent"] += 1
                continue

            digest = marker_hash(interaction)
            if digest in seen:
                counts["duplicate"] += 1
                continue

            seen[digest] = interaction_name
            counts["valid"] += 1
            rows.append(
                {
                    "interaction": interaction_name,
                    "parent_a": parent_a,
                    "parent_b": parent_b,
                    "n_present": int(interaction.sum()),
                    "frequency": float(interaction.mean()),
                    "hash": digest,
                }
            )
            if args.write_matrix:
                matrix[interaction_name] = interaction

    meta = pd.DataFrame(rows)
    meta.to_csv(args.out_dir / "corrected_epistasis_interactions.csv", index=False)
    if args.write_matrix:
        pd.DataFrame(matrix).to_csv(args.out_dir / "corrected_epistasis_interaction_matrix.csv", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "additive_markers", "value": len(names)},
            {"metric": "possible_pairwise_interactions", "value": math.comb(len(names), 2)},
            *[{"metric": key, "value": value} for key, value in counts.items()],
        ]
    )
    summary.to_csv(args.out_dir / "corrected_epistasis_interaction_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
