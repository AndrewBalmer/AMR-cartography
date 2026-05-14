#!/usr/bin/env python3
"""Wrapper for a corrected-panel epistasis permutation chunk."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--interaction-file", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--array-index", type=int)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--permutation-index", required=True, type=int)
    parser.add_argument("--base-seed", default=20260514, type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script = Path(__file__).with_name("run_epistasis_chunk.py")
    seed = args.base_seed + args.permutation_index
    cmd = [
        sys.executable,
        str(script),
        "--data-dir",
        str(args.data_dir),
        "--interaction-file",
        str(args.interaction_file),
        "--out-dir",
        str(args.out_dir),
        "--array-index",
        str(args.array_index),
        "--chunk-size",
        str(args.chunk_size),
        "--chunk-label",
        f"perm_{args.permutation_index:04d}_array_{args.array_index:06d}",
        "--permutation-seed",
        str(seed),
        "--pvalues-only",
    ]
    if args.force:
        cmd.append("--force")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
