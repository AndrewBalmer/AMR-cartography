#!/usr/bin/env python3
"""Build recomputed-threshold manifest and old-vs-new manuscript statistics."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from common import HISTORICAL_PUBLIC_SUPPLEMENT, load_historical_public_supplement


EVIDENCE_ORDER = ["Very Strong", "Strong", "Moderate", "Weak", "Weak/No Evidence"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--original-logic-public", type=Path)
    parser.add_argument(
        "--historical-public-supplement",
        default=HISTORICAL_PUBLIC_SUPPLEMENT,
        type=Path,
        help="Frozen, checksum-verified copy of the historical 354-row public Supplementary File 1.",
    )
    parser.add_argument("--out-dir", type=Path)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {"missing": str(path)}


def load_historical_public(path: Path | None = None) -> pd.DataFrame:
    """Historical (157-marker) public table, from the frozen fixture."""
    return load_historical_public_supplement(path)


def evidence_counts(df: pd.DataFrame) -> str:
    counts = df["Evidence"].value_counts().reindex(EVIDENCE_ORDER, fill_value=0)
    return "/".join(str(int(counts[label])) for label in EVIDENCE_ORDER)


def position_series(df: pd.DataFrame) -> pd.Series:
    return df["Substitution"].astype(str).str.extract(r"([0-9]+)", expand=False)


def table_stats(df: pd.DataFrame) -> dict[str, object]:
    any_evidence = df[df["Evidence"].isin(["Very Strong", "Strong", "Moderate", "Weak"])].copy()
    multi = df[df["Evidence"].isin(["Very Strong", "Strong", "Moderate"])].copy()
    any_positions = any_evidence.assign(position=position_series(any_evidence))[["PBP", "position"]].drop_duplicates()
    multi_positions = multi.assign(position=position_series(multi))[["PBP", "position"]].drop_duplicates()
    beta = pd.to_numeric(multi["β Joint"], errors="coerce")
    return {
        "public_rows": int(len(df)),
        "evidence_counts": evidence_counts(df),
        "any_evidence_rows": int(len(any_evidence)),
        "any_unique_substitutions": int(any_evidence[["PBP", "Substitution"]].drop_duplicates().shape[0]),
        "any_unique_positions": int(any_positions.shape[0]),
        "multi_rows": int(len(multi)),
        "multi_unique_substitutions": int(multi[["PBP", "Substitution"]].drop_duplicates().shape[0]),
        "multi_unique_positions": int(multi_positions.shape[0]),
        "multi_positions_by_pbp": "/".join(
            str(int((multi_positions["PBP"] == pbp).sum())) for pbp in ["1A", "2B", "2X"]
        ),
        "multi_beta_range": "NA" if beta.dropna().empty else f"{beta.min():.3f}-{beta.max():.3f}",
    }


def markdown_table(rows: list[dict[str, object]]) -> list[str]:
    metrics = [
        ("public_rows", "Public rows"),
        ("evidence_counts", "Evidence VS/S/M/W/None"),
        ("any_evidence_rows", "Any-evidence rows"),
        ("any_unique_substitutions", "Any-evidence unique substitutions"),
        ("any_unique_positions", "Any-evidence unique positions"),
        ("multi_rows", "Multi-method rows"),
        ("multi_unique_substitutions", "Multi-method unique substitutions"),
        ("multi_unique_positions", "Multi-method unique positions"),
        ("multi_positions_by_pbp", "Multi positions 1A/2B/2X"),
        ("multi_beta_range", "Multi beta range"),
    ]
    names = [row["name"] for row in rows]
    out = ["| Metric | " + " | ".join(names) + " |", "|---|" + "|".join(["---"] * len(names)) + "|"]
    for key, label in metrics:
        out.append("| " + label + " | " + " | ".join(str(row["stats"].get(key, "NA")) for row in rows) + " |")
    return out


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir or (args.results_dir / "manuscript_outputs")
    out_dir.mkdir(parents=True, exist_ok=True)

    threshold_dir = args.results_dir / "thresholds"
    epistasis_dir = args.results_dir / "epistasis" / "merged"
    manifest = {
        "threshold_policy": "lowest permutation minimum",
        "meff": load_json(threshold_dir / "recomputed_meff.json"),
        "additive": load_json(threshold_dir / "additive_threshold_summary.json"),
        "uvlmm": load_json(threshold_dir / "uvlmm_summary.json"),
        "epistasis": load_json(epistasis_dir / "corrected_epistasis_summary.json"),
    }
    (threshold_dir / "recomputed_thresholds.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    rows: list[dict[str, object]] = []
    rows.append(
        {
            "name": "historical old preprint (frozen fixture)",
            "stats": table_stats(load_historical_public(args.historical_public_supplement)),
        }
    )

    original_public = args.original_logic_public
    if original_public is not None and original_public.exists():
        rows.append({"name": "historical-threshold corrected panel", "stats": table_stats(pd.read_csv(original_public))})

    recomputed_public = out_dir / "Supplementary_File_1.csv"
    if recomputed_public.exists():
        rows.append({"name": "recomputed 170-marker analysis", "stats": table_stats(pd.read_csv(recomputed_public))})

    report = [
        "# Recomputed 170-marker manuscript statistics",
        "",
        "Threshold policy: lowest observed permutation minimum for each null family.",
        "",
        "## Threshold Manifest",
        "",
        "```json",
        json.dumps(manifest, indent=2, sort_keys=True),
        "```",
        "",
        "## Headline Comparison",
        "",
        *markdown_table(rows),
        "",
    ]
    (out_dir / "recomputed_170_manuscript_statistics.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
