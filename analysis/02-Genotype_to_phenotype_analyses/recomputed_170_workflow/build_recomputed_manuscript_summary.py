#!/usr/bin/env python3
"""Build old-vs-recomputed manuscript result summary artifacts.

This is a reporting script, not a model-fitting step. It compares the old
preprint public table with the recomputed public table, writes gained/lost
multi-method rows, and drafts manuscript replacement text for review.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd


EVIDENCE_ORDER = ["Very Strong", "Strong", "Moderate", "Weak", "Weak/No Evidence"]
ANY_EVIDENCE = ["Very Strong", "Strong", "Moderate", "Weak"]
MULTI_METHOD = ["Very Strong", "Strong", "Moderate"]
PBP_ORDER = ["1A", "2B", "2X"]
VARIABLE_PBP_POSITIONS = 285


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument(
        "--historical-public",
        type=Path,
        default=Path("analysis_outputs/original_logic_rebuild/manuscript_outputs/Supplementary_File_1.csv"),
    )
    parser.add_argument("--out-dir", type=Path)
    return parser.parse_args()


def load_main_public() -> pd.DataFrame:
    text = subprocess.check_output(["git", "show", "main:manuscript/Supplementary_File_1.csv"], text=True)
    return pd.read_csv(StringIO(text))


def position_series(df: pd.DataFrame) -> pd.Series:
    return df["Substitution"].astype(str).str.extract(r"([0-9]+)", expand=False)


def evidence_counts(df: pd.DataFrame) -> dict[str, int]:
    counts = df["Evidence"].value_counts().reindex(EVIDENCE_ORDER, fill_value=0)
    return {label: int(counts[label]) for label in EVIDENCE_ORDER}


def summarize(df: pd.DataFrame) -> dict[str, object]:
    any_df = df[df["Evidence"].isin(ANY_EVIDENCE)].copy()
    multi_df = df[df["Evidence"].isin(MULTI_METHOD)].copy()
    any_positions = any_df.assign(position=position_series(any_df))[["PBP", "position"]].drop_duplicates()
    multi_positions = multi_df.assign(position=position_series(multi_df))[["PBP", "position"]].drop_duplicates()
    beta = pd.to_numeric(multi_df["β Joint"], errors="coerce")
    return {
        "public_rows": int(len(df)),
        "evidence_counts": evidence_counts(df),
        "any_rows": int(len(any_df)),
        "any_unique_substitutions": int(any_df[["PBP", "Substitution"]].drop_duplicates().shape[0]),
        "any_unique_positions": int(any_positions.shape[0]),
        "any_positions_by_pbp": {pbp: int((any_positions["PBP"] == pbp).sum()) for pbp in PBP_ORDER},
        "any_position_fraction_of_285": float(any_positions.shape[0] / VARIABLE_PBP_POSITIONS),
        "multi_rows": int(len(multi_df)),
        "multi_unique_substitutions": int(multi_df[["PBP", "Substitution"]].drop_duplicates().shape[0]),
        "multi_unique_positions": int(multi_positions.shape[0]),
        "multi_positions_by_pbp": {pbp: int((multi_positions["PBP"] == pbp).sum()) for pbp in PBP_ORDER},
        "multi_rows_by_pbp": {pbp: int((multi_df["PBP"] == pbp).sum()) for pbp in PBP_ORDER},
        "multi_beta_min": None if beta.dropna().empty else float(beta.min()),
        "multi_beta_max": None if beta.dropna().empty else float(beta.max()),
    }


def key_df(df: pd.DataFrame) -> pd.DataFrame:
    return df[["PBP", "Substitution"]].drop_duplicates().copy()


def multi_df(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["Evidence"].isin(MULTI_METHOD)].copy()


def diff_rows(new: pd.DataFrame, old: pd.DataFrame, *, kind: str) -> pd.DataFrame:
    new_multi = multi_df(new)
    old_multi = multi_df(old)
    if kind == "gained":
        keys = key_df(new_multi).merge(key_df(old_multi), on=["PBP", "Substitution"], how="left", indicator=True)
        keys = keys[keys["_merge"].eq("left_only")].drop(columns=["_merge"])
        return new_multi.merge(keys, on=["PBP", "Substitution"], how="inner")
    if kind == "lost":
        keys = key_df(old_multi).merge(key_df(new_multi), on=["PBP", "Substitution"], how="left", indicator=True)
        keys = keys[keys["_merge"].eq("left_only")].drop(columns=["_merge"])
        return old_multi.merge(keys, on=["PBP", "Substitution"], how="inner")
    raise ValueError(kind)


def changed_shared_rows(new: pd.DataFrame, old: pd.DataFrame) -> pd.DataFrame:
    cols = ["Evidence", "PBP", "Substitution", "β Joint", "Adj. p-value", "No. Sig. interactions", "total"]
    shared = key_df(multi_df(new)).merge(key_df(multi_df(old)), on=["PBP", "Substitution"], how="inner")
    left = old[cols].merge(shared, on=["PBP", "Substitution"], how="inner").add_prefix("old_")
    right = new[cols].merge(shared, on=["PBP", "Substitution"], how="inner").add_prefix("new_")
    out = left.merge(
        right,
        left_on=["old_PBP", "old_Substitution"],
        right_on=["new_PBP", "new_Substitution"],
        how="inner",
    )
    comparable = ["Evidence", "β Joint", "Adj. p-value", "No. Sig. interactions", "total"]
    mask = False
    for col in comparable:
        mask = mask | (out[f"old_{col}"].astype(str) != out[f"new_{col}"].astype(str))
    return out[mask].copy()


def markdown_table(rows: list[list[object]]) -> str:
    header = rows[0]
    body = rows[1:]
    out = ["| " + " | ".join(map(str, header)) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for row in body:
        out.append("| " + " | ".join(map(str, row)) + " |")
    return "\n".join(out)


def frame_to_markdown(df: pd.DataFrame) -> str:
    rows: list[list[object]] = [df.columns.tolist()]
    rows.extend(df.fillna("").astype(str).values.tolist())
    return markdown_table(rows)


def evidence_count_string(summary: dict[str, object]) -> str:
    counts = summary["evidence_counts"]
    return "/".join(str(counts[label]) for label in EVIDENCE_ORDER)


def pbp_string(values: dict[str, int]) -> str:
    return "/".join(str(values[pbp]) for pbp in PBP_ORDER)


def fmt_beta(summary: dict[str, object]) -> str:
    return f"{summary['multi_beta_min']:.3f}-{summary['multi_beta_max']:.3f}"


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {}


def build_replacement_text(
    summaries: dict[str, dict[str, object]],
    recomputed: pd.DataFrame,
    manifest: dict[str, object],
    marker_summary: dict[str, object],
) -> str:
    old = summaries["old_preprint"]
    hist = summaries["historical_threshold"]
    rec = summaries["recomputed_170"]
    top = multi_df(recomputed).head(1).iloc[0]
    strong = multi_df(recomputed)[multi_df(recomputed)["Evidence"].eq("Strong")]
    strong_list = ", ".join(f"PBP{row.PBP} {row.Substitution}" for row in strong.itertuples(index=False))
    marker_multi = int(sum(marker_summary[label] for label in ["Very Strong", "Strong", "Moderate"]))
    epistasis = manifest.get("epistasis", {})
    additive = manifest.get("additive", {})

    return "\n".join(
        [
            "# Manuscript Replacement Text Draft",
            "",
            "Use this as draft wording only; review against the manuscript before pasting.",
            "",
            "## Results: Evidence Overview",
            "",
            (
                "In the recomputed 170-marker analysis, Supplementary File 1 contains "
                f"{rec['public_rows']} component-expanded rows. "
                f"{rec['any_rows']} rows show at least one evidence stream, corresponding to "
                f"{rec['any_unique_substitutions']} unique PBP substitutions at "
                f"{rec['any_unique_positions']} unique PBP positions "
                f"({rec['any_position_fraction_of_285'] * 100:.1f}% of the original "
                f"{VARIABLE_PBP_POSITIONS} variable PBP positions)."
            ),
            "",
            (
                f"{rec['multi_rows']} component-expanded substitutions at "
                f"{rec['multi_unique_positions']} PBP positions were supported by two or more methods "
                f"({rec['multi_positions_by_pbp']['1A']} PBP1A, "
                f"{rec['multi_positions_by_pbp']['2B']} PBP2B, "
                f"{rec['multi_positions_by_pbp']['2X']} PBP2X positions). "
                f"In the separate marker-level table, {marker_multi} of 170 markers had "
                "two or more evidence streams. "
                f"This compares with {old['multi_rows']} substitutions at {old['multi_unique_positions']} positions "
                "in the old preprint and "
                f"{hist['multi_rows']} substitutions at {hist['multi_unique_positions']} positions "
                "in the historical-threshold corrected-panel rebuild."
            ),
            "",
            (
                f"The sole Very Strong substitution remains PBP{top['PBP']} {top['Substitution']} "
                f"(beta Joint {top['β Joint']}, adjusted p-value {top['Adj. p-value']}, "
                f"{top['No. Sig. interactions']} significant interactions). "
                f"The Strong substitutions remain: {strong_list}."
            ),
            "",
            (
                f"Among substitutions supported by two or more methods, joint mvLMM effect sizes ranged from "
                f"{rec['multi_beta_min']:.3f} to {rec['multi_beta_max']:.3f} MIC units."
            ),
            "",
            "## Results: Evidence Counts",
            "",
            (
                "Public component-expanded evidence counts are: "
                f"Very Strong {rec['evidence_counts']['Very Strong']}, "
                f"Strong {rec['evidence_counts']['Strong']}, "
                f"Moderate {rec['evidence_counts']['Moderate']}, "
                f"Weak {rec['evidence_counts']['Weak']}, and "
                f"Weak/No Evidence {rec['evidence_counts']['Weak/No Evidence']}."
            ),
            "",
            (
                "Marker-level evidence counts are: "
                f"Very Strong {marker_summary['Very Strong']}, "
                f"Strong {marker_summary['Strong']}, "
                f"Moderate {marker_summary['Moderate']}, "
                f"Weak {marker_summary['Weak']}, and "
                f"Weak/No Evidence {marker_summary['Weak/No Evidence']}."
            ),
            "",
            "## Methods/Provenance Sentence",
            "",
            (
                "For the revised 170-marker analysis, Galwey effective-test counts and empirical thresholds "
                "were recomputed on the corrected marker and interaction universes. The additive threshold "
                f"was the lowest of 100 random-phenotype minimum p-values "
                f"(raw {additive.get('raw_threshold')}, Galwey-adjusted {additive.get('adjusted_threshold')}); "
                "the epistasis threshold was the lowest of 100 permutation minimum p-values "
                f"(raw {epistasis.get('raw_lowest_min_p_threshold')}, "
                f"Galwey-adjusted {epistasis.get('epistasis_threshold')})."
            ),
            "",
            "## Table 3",
            "",
            "Use `recomputed_170_table3_top20.csv` as the recomputed Table 3 candidate row source.",
            "",
            "## Important Frame Note",
            "",
            (
                "The counts above use the original manuscript's component-expanded public frame. "
                "Marker-level counts should be discussed only where explicitly labelled as marker-level."
            ),
            "",
        ]
    )


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir or (args.results_dir / "manuscript_outputs" / "manuscript_summary")
    out_dir.mkdir(parents=True, exist_ok=True)

    recomputed_public = pd.read_csv(args.results_dir / "manuscript_outputs" / "Supplementary_File_1.csv")
    marker_level = pd.read_csv(args.results_dir / "manuscript_outputs" / "Supplementary_File_1_corrected_marker_level.csv")
    historical_public = pd.read_csv(args.historical_public)
    old_public = load_main_public()
    manifest = load_json(args.results_dir / "thresholds" / "recomputed_thresholds.json")

    tables = {
        "old_preprint": old_public,
        "historical_threshold": historical_public,
        "recomputed_170": recomputed_public,
    }
    summaries = {name: summarize(df) for name, df in tables.items()}

    comparisons = [
        ("old_preprint", old_public, "old_preprint"),
        ("historical_threshold", historical_public, "historical_threshold_corrected"),
    ]
    for _, baseline, label in comparisons:
        diff_rows(recomputed_public, baseline, kind="gained").to_csv(
            out_dir / f"gained_multi_vs_{label}.csv", index=False
        )
        diff_rows(recomputed_public, baseline, kind="lost").to_csv(
            out_dir / f"lost_multi_vs_{label}.csv", index=False
        )
        changed_shared_rows(recomputed_public, baseline).to_csv(
            out_dir / f"shared_multi_value_changes_vs_{label}.csv", index=False
        )

    table3 = multi_df(recomputed_public).head(20)
    table3.to_csv(out_dir / "recomputed_170_table3_top20.csv", index=False)

    marker_counts = marker_level["Evidence"].value_counts().reindex(EVIDENCE_ORDER, fill_value=0).to_dict()
    replacement = build_replacement_text(summaries, recomputed_public, manifest, marker_counts)
    (out_dir / "recomputed_170_manuscript_replacement_text.md").write_text(replacement)

    headline_rows = [
        ["Metric", "Old preprint", "Historical-threshold corrected", "Recomputed 170-marker"],
        [
            "Public rows",
            summaries["old_preprint"]["public_rows"],
            summaries["historical_threshold"]["public_rows"],
            summaries["recomputed_170"]["public_rows"],
        ],
        [
            "Evidence VS/S/M/W/None",
            evidence_count_string(summaries["old_preprint"]),
            evidence_count_string(summaries["historical_threshold"]),
            evidence_count_string(summaries["recomputed_170"]),
        ],
        [
            "Any-evidence rows/subs/positions",
            f"{summaries['old_preprint']['any_rows']}/{summaries['old_preprint']['any_unique_substitutions']}/{summaries['old_preprint']['any_unique_positions']}",
            f"{summaries['historical_threshold']['any_rows']}/{summaries['historical_threshold']['any_unique_substitutions']}/{summaries['historical_threshold']['any_unique_positions']}",
            f"{summaries['recomputed_170']['any_rows']}/{summaries['recomputed_170']['any_unique_substitutions']}/{summaries['recomputed_170']['any_unique_positions']}",
        ],
        [
            "Any-evidence positions 1A/2B/2X",
            pbp_string(summaries["old_preprint"]["any_positions_by_pbp"]),
            pbp_string(summaries["historical_threshold"]["any_positions_by_pbp"]),
            pbp_string(summaries["recomputed_170"]["any_positions_by_pbp"]),
        ],
        [
            "Multi-method rows/subs/positions",
            f"{summaries['old_preprint']['multi_rows']}/{summaries['old_preprint']['multi_unique_substitutions']}/{summaries['old_preprint']['multi_unique_positions']}",
            f"{summaries['historical_threshold']['multi_rows']}/{summaries['historical_threshold']['multi_unique_substitutions']}/{summaries['historical_threshold']['multi_unique_positions']}",
            f"{summaries['recomputed_170']['multi_rows']}/{summaries['recomputed_170']['multi_unique_substitutions']}/{summaries['recomputed_170']['multi_unique_positions']}",
        ],
        [
            "Multi positions 1A/2B/2X",
            pbp_string(summaries["old_preprint"]["multi_positions_by_pbp"]),
            pbp_string(summaries["historical_threshold"]["multi_positions_by_pbp"]),
            pbp_string(summaries["recomputed_170"]["multi_positions_by_pbp"]),
        ],
        [
            "Multi beta range",
            fmt_beta(summaries["old_preprint"]),
            fmt_beta(summaries["historical_threshold"]),
            fmt_beta(summaries["recomputed_170"]),
        ],
    ]

    report = [
        "# Recomputed 170-Marker Manuscript Result Summary",
        "",
        "All counts use the original component-expanded public Supplementary File 1 frame unless labelled otherwise.",
        "",
        "## Headline Counts",
        "",
        markdown_table(headline_rows),
        "",
        "## Gained/Lost Multi-Method Rows",
        "",
        "| Comparison | Gained | Lost | Shared rows with changed displayed values |",
        "|---|---:|---:|---:|",
    ]
    for _, baseline, label in comparisons:
        gained = len(diff_rows(recomputed_public, baseline, kind="gained"))
        lost = len(diff_rows(recomputed_public, baseline, kind="lost"))
        changed = len(changed_shared_rows(recomputed_public, baseline))
        report.append(f"| recomputed vs {label} | {gained} | {lost} | {changed} |")

    report.extend(
        [
            "",
            "## Generated Files",
            "",
            "- `gained_multi_vs_old_preprint.csv`",
            "- `lost_multi_vs_old_preprint.csv`",
            "- `shared_multi_value_changes_vs_old_preprint.csv`",
            "- `gained_multi_vs_historical_threshold_corrected.csv`",
            "- `lost_multi_vs_historical_threshold_corrected.csv`",
            "- `shared_multi_value_changes_vs_historical_threshold_corrected.csv`",
            "- `recomputed_170_table3_top20.csv`",
            "- `recomputed_170_manuscript_replacement_text.md`",
            "",
            "## Table 3 Top 20 Preview",
            "",
            frame_to_markdown(
                table3[["Evidence", "PBP", "Substitution", "β Joint", "Adj. p-value", "No. Sig. interactions", "total"]]
            ),
            "",
            "## Replacement Text",
            "",
            replacement,
        ]
    )
    (out_dir / "recomputed_170_manuscript_result_summary.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
