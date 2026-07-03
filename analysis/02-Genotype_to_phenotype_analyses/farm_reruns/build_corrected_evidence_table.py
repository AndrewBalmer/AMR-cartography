#!/usr/bin/env python3
"""Build manuscript-facing evidence tables after corrected farm reruns."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    HISTORICAL_ADDITIVE_GALWEY_MEFF,
    HISTORICAL_ADDITIVE_THRESHOLD,
    HISTORICAL_UV_THRESHOLD,
    add_galwey_adjusted_pvalues,
    marker_order_from_effects,
    read_csv,
)


LEGACY_PUBLIC_COLUMNS = [
    "Evidence",
    "PBP",
    "Substitution",
    "Phenotypic Distance",
    "SE",
    "No. isolates",
    "Confidence",
    "β Joint",
    "Adj. p-value",
    "Effect Size Axis 1",
    "Effect Size Axis 2",
    "No. Sig. interactions",
    "total",
    "Sig. mvLMM/uvLMM (No. drugs",
]

PBP_DISPLAY = {"PBP1A": "1A", "PBP2B": "2B", "PBP2X": "2X"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--new-dir", required=True, type=Path)
    parser.add_argument("--support-dir", required=True, type=Path)
    parser.add_argument("--analysis-out", required=True, type=Path)
    parser.add_argument("--additive-dir", type=Path)
    parser.add_argument("--uv-dir", type=Path)
    parser.add_argument("--epistasis-dir", type=Path)
    parser.add_argument("--threshold-file", type=Path)
    parser.add_argument("--analysis-label", default="Corrected rerun")
    parser.add_argument(
        "--output-dir",
        default=Path("farm_outputs/original_logic_rebuild/manuscript_outputs"),
        type=Path,
        help="Directory for rebuilt outputs. Defaults away from manuscript/ until validation passes.",
    )
    parser.add_argument(
        "--manuscript-dir",
        type=Path,
        default=None,
        help="Deprecated alias for --output-dir; use only for the final deliberate replacement step.",
    )
    parser.add_argument("--additive-threshold", default=HISTORICAL_ADDITIVE_THRESHOLD, type=float)
    parser.add_argument("--uv-threshold", default=HISTORICAL_UV_THRESHOLD, type=float)
    parser.add_argument("--additive-galwey-meff", default=HISTORICAL_ADDITIVE_GALWEY_MEFF, type=float)
    return parser.parse_args()


def parse_marker(marker: str) -> list[dict[str, object]]:
    return [
        {
            "protein": prot,
            "aa": f"{aa1}{pos}{aa2}",
            "pos": int(pos),
            "pos_key": f"{prot}_{int(pos)}",
            "loc": f"{prot}_{aa1}{pos}",
            "full": f"{prot}_{aa1}{pos}_{aa2}",
        }
        for prot, aa1, pos, aa2 in re.findall(r"(PBP[12][ABX])_([A-Z])([0-9]+)_([A-Z])", str(marker))
    ]


def compact_label(marker: str) -> str:
    parts = parse_marker(marker)
    return "/".join(str(part["aa"]) for part in parts) if parts else str(marker)


def primary_pbp(marker: str) -> str:
    parts = parse_marker(marker)
    proteins = list(dict.fromkeys(str(part["protein"]) for part in parts))
    return proteins[0] if len(proteins) == 1 else ("Mixed" if proteins else "")


def split_location_rows(df: pd.DataFrame, column: str = "Location") -> pd.DataFrame:
    out = df.copy()
    out[column] = [
        str(value).split("/") if pd.notna(value) else [np.nan]
        for value in out[column]
    ]
    return out.explode(column, ignore_index=True)


def as_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def column_or_na(df: pd.DataFrame, column: str) -> pd.Series:
    return df[column] if column in df.columns else pd.Series(np.nan, index=df.index)


def format_count(value: object) -> str:
    if pd.isna(value):
        return "-"
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(number) and float(number).is_integer():
        return str(int(number))
    return str(value)


def format_decimal(value: object, digits: int = 3, missing: str = "-") -> str:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number):
        return missing
    text = f"{float(number):.{digits}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def format_p_value(value: object) -> str:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number):
        return "Not Tested Due to Low Sample Size"
    number = min(float(number), 1.0)
    return "<1e-16" if number < 1e-16 else f"{float(number):.3e}"


def format_interaction_count(value: object) -> str:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number) or int(number) == 0:
        return "-"
    return str(int(number))


def evidence_category(total: int) -> str:
    return {0: "Weak/No Evidence", 1: "Weak", 2: "Moderate", 3: "Strong", 4: "Very Strong"}.get(
        total, "Very Strong"
    )


def load_additive(new_dir: Path, threshold: float, galwey_meff: float, additive_dir: Path | None = None) -> pd.DataFrame:
    source_dir = additive_dir or new_dir
    pvals = read_csv(source_dir / "mvLMM_p_values_normal_pneumo_low_freq_vars.csv")
    unnamed = [col for col in pvals.columns if str(col).startswith("Unnamed")]
    if unnamed:
        pvals = pvals.rename(columns={unnamed[0]: "mvLMM_result_index"})
    effect_file = source_dir / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv"
    if "marker" not in pvals.columns:
        order = marker_order_from_effects(effect_file)
        if len(order) != len(pvals):
            raise ValueError("Additive mvLMM p-values/effects do not align")
        pvals["marker"] = order

    effects = read_csv(effect_file)
    cand = effects[effects["effect_type"].eq("candidate")][["effect_name", "env", "effsize"]]
    wide = cand.pivot(index="effect_name", columns="env", values="effsize").reset_index()
    wide = wide.rename(columns={"effect_name": "marker", "env1_D1": "effect_axis1", "env1_D2": "effect_axis2"})

    out = pvals.merge(wide, on="marker", how="left")
    if "pv20_raw" not in out.columns:
        out = out.rename(columns={"pv20": "pv20_raw"})
    out = add_galwey_adjusted_pvalues(out, raw_column="pv20_raw", meff=galwey_meff)
    out["joint_effect_size"] = np.sqrt(out["effect_axis1"] ** 2 + out["effect_axis2"] ** 2)
    out["mvLMM_significant_historical"] = out["pv20_adj_galwey"] <= threshold
    out["mvLMM_significant_p001"] = out["pv20_adj_galwey"] <= HISTORICAL_UV_THRESHOLD
    out["mvLMM_evidence"] = out["mvLMM_significant_historical"] & (out["joint_effect_size"] >= 1)
    out["PBP"] = out["marker"].map(primary_pbp)
    out["Substitution"] = out["marker"].map(compact_label)
    out["constituents"] = out["marker"].map(
        lambda marker: ";".join(f"{part['protein']}:{part['aa']}" for part in parse_marker(marker))
    )
    return out


def load_uv_support(
    analysis_out: Path,
    legacy_dir: Path,
    new_dir: Path,
    threshold: float,
    galwey_meff: float,
    uv_dir: Path | None = None,
) -> pd.DataFrame:
    support_path = (uv_dir or analysis_out) / "uniLMM_exact_170_marker_support.csv"
    if support_path.exists():
        support = read_csv(support_path)
        return support.rename(
            columns={
                "uvLMM_exact_n_drugs": "uvLMM_exact_n_drugs_p001",
                "uvLMM_exact_drugs": "uvLMM_exact_drugs_p001",
            }
        )

    old_p = read_csv(legacy_dir / "uniLMM_p_val_normal_MIC_pneumo.csv")
    old_e = read_csv(legacy_dir / "uniLMM_effect_normal_MIC_pneumo.csv")
    old = old_e[old_e["effect_type"].eq("candidate")][["trait", "effect_name"]].reset_index(drop=True)
    old = old.rename(columns={"trait": "drug", "effect_name": "marker"})
    old["pv20"] = old_p["pv20"].to_numpy()
    old["source"] = "historical_157"

    added = read_csv(new_dir / "uniLMM_exact_added_markers_p_values.csv")[["marker", "drug", "pv20"]]
    added["source"] = "rerun_added_13"
    uv = pd.concat([old, added], ignore_index=True)
    uv = add_galwey_adjusted_pvalues(uv, meff=galwey_meff)

    rows = []
    for marker, group in uv.groupby("marker"):
        sig = sorted(group.loc[group["pv20_adj_galwey"] <= threshold, "drug"].astype(str).unique())
        rows.append(
            {
                "marker": marker,
                "uvLMM_exact_n_drugs_p001": len(sig),
                "uvLMM_exact_drugs_p001": ";".join(sig),
                "uvLMM_exact_source": ";".join(sorted(group["source"].unique())),
            }
        )
    return pd.DataFrame(rows)


def load_component_support(support_dir: Path) -> pd.Series:
    single = read_csv(support_dir / "Single_subs_all_S.pneumo.csv")
    single["PBP"] = single["PBP"].astype(str).str.strip()
    single["aa"] = (
        single["Substitution_of_interest"]
        .astype(str)
        .str.replace(r"\([0-9]+\)", "", regex=True)
        .str.replace(r"\s+", "", regex=True)
    )
    single["single_key"] = single["PBP"] + ":" + single["aa"]
    single["median_phenotypic_distance"] = as_numeric(single["median_phenotypic_distance"])
    single_hits = set(single.loc[single["median_phenotypic_distance"] >= 1, "single_key"])

    cluster = load_legacy_clustering(support_dir)
    cluster_hits = set(
        cluster.loc[cluster["Confidence"].eq("Strong"), "PBP"].astype(str)
        + "_"
        + cluster.loc[cluster["Confidence"].eq("Strong"), "Loci"].astype(str)
    )

    prior = read_csv(support_dir / "CDC_GWAS_overlap_TPD.csv")
    prior["PBP"] = prior["PBP"].replace({"1a": "PBP1A", "2b": "PBP2B", "2x": "PBP2X"})
    prior["pos_key"] = prior.apply(lambda row: f"{row['PBP']}_{int(row['Position'])}", axis=1)
    for col in ["GWAS", "CDC", "Laboratory"]:
        prior[col] = prior[col].astype(str).str.lower().eq("yes")
    prior_hits = set(prior.loc[prior[["GWAS", "CDC", "Laboratory"]].any(axis=1), "pos_key"])

    return pd.Series({"single_hits": single_hits, "cluster_hits": cluster_hits, "prior_hits": prior_hits})


def load_legacy_single_subs(support_dir: Path) -> pd.DataFrame:
    single = read_csv(support_dir / "Single_subs_all_S.pneumo.csv").copy()
    single["Location"] = single["AA_1"].astype(str) + single["Loci"].astype(str) + single["AA_2"].astype(str)
    single["Single_subs"] = True
    return single


def load_legacy_clustering(support_dir: Path) -> pd.DataFrame:
    cluster_files = sorted(support_dir.glob("Sig_AA_subs_Cluster_gen*_S.pneumo.csv"))
    if not cluster_files:
        raise FileNotFoundError(f"No clustering support files found in {support_dir}")
    cluster = pd.concat([read_csv(path) for path in cluster_files], ignore_index=True)
    cluster = cluster.sort_values("Confidence", kind="stable").drop_duplicates("Location", keep="first")
    cluster[["PBP", "Loci"]] = cluster["Location"].astype(str).str.split("_", n=1, expand=True)
    return cluster.drop(columns=["Location"])


def additive_for_legacy(additive: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PBP": additive["PBP"],
            "Location": additive["Substitution"],
            "marker_mv_LMM": additive["marker"],
            "pv20_raw_mv_LMM": additive["pv20_raw"],
            "pv20_adj_galwey_mv_LMM": additive["pv20_adj_galwey"],
            "Joint_effsize_mv_LMM": additive["joint_effect_size"],
            "effsize_D1_mv_LMM": additive["effect_axis1"],
            "effsize_D2_mv_LMM": additive["effect_axis2"],
            "mvLMM_significant_historical_mv_LMM": additive["mvLMM_significant_historical"],
            "mvLMM_significant_p001_mv_LMM": additive["mvLMM_significant_p001"],
            "mvLMM_evidence_mv_LMM": additive["mvLMM_evidence"],
        }
    )


def epistasis_for_legacy(epi: pd.DataFrame) -> pd.DataFrame:
    if epi.empty:
        return pd.DataFrame(
            columns=[
                "PBP",
                "Location",
            "marker_mv_LMM_epistatic",
            "PBP_count_mv_LMM_epistatic",
            "min_epistasis_pv20_mv_LMM_epistatic",
            "min_epistasis_pv20_adj_galwey_mv_LMM_epistatic",
        ]
        )

    out = epi.copy()
    out["PBP"] = out["marker"].map(primary_pbp)
    out["Location"] = out["marker"].map(compact_label)
    out = out[out["PBP"].isin(PBP_DISPLAY)].copy()
    out = split_location_rows(out, "Location")
    return out.rename(
        columns={
            "marker": "marker_mv_LMM_epistatic",
            "num_sig_interactions": "PBP_count_mv_LMM_epistatic",
            "min_epistasis_pv20": "min_epistasis_pv20_mv_LMM_epistatic",
            "min_epistasis_pv20_adj_galwey": "min_epistasis_pv20_adj_galwey_mv_LMM_epistatic",
        }
    )


def build_marker_level_table(additive: pd.DataFrame, uv: pd.DataFrame, component: pd.Series, epi: pd.DataFrame) -> pd.DataFrame:
    table = additive.merge(uv, on="marker", how="left").merge(epi, on="marker", how="left")
    table["uvLMM_exact_n_drugs_p001"] = as_numeric(table["uvLMM_exact_n_drugs_p001"]).fillna(0).astype(int)
    table["num_sig_interactions"] = as_numeric(table["num_sig_interactions"]).fillna(0).astype(int)
    table["single_sub_evidence"] = table["marker"].map(
        lambda marker: any(f"{part['protein']}:{part['aa']}" in component["single_hits"] for part in parse_marker(marker))
    )
    table["cluster_evidence"] = table["marker"].map(
        lambda marker: any(str(part["loc"]) in component["cluster_hits"] for part in parse_marker(marker))
    )
    table["prior_any_position"] = table["marker"].map(
        lambda marker: any(str(part["pos_key"]) in component["prior_hits"] for part in parse_marker(marker))
    )
    table["epistasis_evidence"] = table["num_sig_interactions"] > 0
    table["evidence_total"] = (
        table["single_sub_evidence"].astype(int)
        + table["cluster_evidence"].astype(int)
        + table["mvLMM_evidence"].astype(int)
        + table["epistasis_evidence"].astype(int)
    )
    table["Evidence"] = table["evidence_total"].map(evidence_category)
    table["Sig. mvLMM/uvLMM (No. drugs)"] = table.apply(
        lambda row: f"{'Yes' if row['mvLMM_significant_p001'] else 'No'}/"
        f"{'Yes' if row['uvLMM_exact_n_drugs_p001'] > 0 else 'No'} ({row['uvLMM_exact_n_drugs_p001']})",
        axis=1,
    )
    return table.sort_values(["evidence_total", "joint_effect_size"], ascending=[False, False])


def build_component_sig_summary(additive: pd.DataFrame, uv: pd.DataFrame) -> pd.DataFrame:
    comparison = additive.merge(uv, on="marker", how="left")
    comparison["uvLMM_exact_n_drugs_p001"] = as_numeric(comparison["uvLMM_exact_n_drugs_p001"]).fillna(0).astype(int)
    comparison = comparison[["PBP", "Substitution", "mvLMM_significant_p001", "uvLMM_exact_n_drugs_p001"]].rename(
        columns={"Substitution": "Location"}
    )
    comparison = split_location_rows(comparison, "Location")
    comparison = (
        comparison.groupby(["PBP", "Location"], as_index=False)
        .agg({"mvLMM_significant_p001": "max", "uvLMM_exact_n_drugs_p001": "max"})
    )
    comparison["sig_summary"] = comparison.apply(
        lambda row: f"{'Yes' if bool(row['mvLMM_significant_p001']) else 'No'}/"
        f"{'Yes' if int(row['uvLMM_exact_n_drugs_p001']) > 0 else 'No'} ({int(row['uvLMM_exact_n_drugs_p001'])})",
        axis=1,
    )
    return comparison[["PBP", "Location", "sig_summary"]]


def build_legacy_public_table(
    additive: pd.DataFrame,
    uv: pd.DataFrame,
    epi: pd.DataFrame,
    support_dir: Path,
    additive_threshold: float,
) -> pd.DataFrame:
    single = load_legacy_single_subs(support_dir)
    mv = additive_for_legacy(additive)
    table = pd.merge(single, mv, on=["Location", "PBP"], how="outer")
    table = split_location_rows(table, "Location")
    table["Loci"] = table["Location"].where(table["Location"].isna(), table["Location"].astype(str).str[:-1])

    clustering = load_legacy_clustering(support_dir)
    table = pd.merge(clustering, table, on=["Loci", "PBP"], how="outer")
    table["Location"] = table["Location"].where(table["Location"].notna(), table["Loci"])

    epi_legacy = epistasis_for_legacy(epi)
    table = pd.merge(epi_legacy, table, on=["Location", "PBP"], how="outer")

    median = as_numeric(column_or_na(table, "median_phenotypic_distance"))
    confidence = column_or_na(table, "Confidence").fillna("-")
    additive_p = as_numeric(column_or_na(table, "pv20_adj_galwey_mv_LMM"))
    joint = as_numeric(column_or_na(table, "Joint_effsize_mv_LMM"))
    epi_count = as_numeric(column_or_na(table, "PBP_count_mv_LMM_epistatic"))

    public = pd.DataFrame(index=table.index)
    public["single_sub_evidence"] = (median >= 1).astype(int)
    public["clustering_evidence"] = confidence.eq("Strong").astype(int)
    public["mvLMM_evidence"] = ((additive_p <= additive_threshold) & (joint >= 1)).astype(int)
    public["epistasis_evidence"] = (epi_count >= 1).astype(int)
    public["total"] = public[["single_sub_evidence", "clustering_evidence", "mvLMM_evidence", "epistasis_evidence"]].sum(
        axis=1
    )
    public["Evidence"] = public["total"].map(evidence_category)
    public["PBP"] = table["PBP"].replace(PBP_DISPLAY)
    public["Substitution"] = table["Location"]
    public["Phenotypic Distance"] = median.map(format_decimal)
    public["SE"] = column_or_na(table, "se").map(format_decimal)
    public["No. isolates"] = (
        column_or_na(table, "number_of_isolates_of_reference").map(format_count)
        + "/"
        + column_or_na(table, "number_of_isolates_of_comparison_group").map(format_count)
    )
    public["Confidence"] = confidence
    public["β Joint"] = joint.round(3)
    public["Adj. p-value"] = additive_p.map(format_p_value)
    public["Effect Size Axis 1"] = column_or_na(table, "effsize_D1_mv_LMM").map(format_decimal)
    public["Effect Size Axis 2"] = column_or_na(table, "effsize_D2_mv_LMM").map(format_decimal)
    public["No. Sig. interactions"] = epi_count.map(format_interaction_count)

    comparison = build_component_sig_summary(additive, uv)
    public["_join_pbp"] = table["PBP"].to_numpy()
    public["_join_location"] = table["Location"].to_numpy()
    public = public.merge(
        comparison,
        left_on=["_join_pbp", "_join_location"],
        right_on=["PBP", "Location"],
        how="left",
        suffixes=("", "_comparison"),
    )
    public["Sig. mvLMM/uvLMM (No. drugs"] = public["sig_summary"].fillna("No/No (0)")

    sort_joint = public["β Joint"].fillna(-np.inf)
    public = public.assign(_sort_joint=sort_joint).sort_values(
        ["total", "_sort_joint"], ascending=[False, False], kind="mergesort"
    )
    public = public.drop(
        columns=[
            "single_sub_evidence",
            "clustering_evidence",
            "mvLMM_evidence",
            "epistasis_evidence",
            "_join_pbp",
            "_join_location",
            "PBP_comparison",
            "Location",
            "sig_summary",
            "_sort_joint",
        ],
        errors="ignore",
    )
    return public[LEGACY_PUBLIC_COLUMNS]


def write_marker_level_display(table: pd.DataFrame, manuscript_dir: Path) -> None:
    table.to_csv(manuscript_dir / "Supplementary_File_1_corrected_marker_level.csv", index=False)


def write_legacy_public(table: pd.DataFrame, manuscript_dir: Path) -> None:
    table.to_csv(manuscript_dir / "Supplementary_File_1.csv", index=False, na_rep="NA")


def build_summary(
    marker_table: pd.DataFrame,
    legacy_table: pd.DataFrame,
    analysis_out: Path,
    additive_threshold: float,
    *,
    analysis_label: str,
    threshold_file: Path | None,
) -> list[str]:
    summary_path = analysis_out / "corrected_epistasis_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    thresholds = json.loads(threshold_file.read_text()) if threshold_file is not None and threshold_file.exists() else {}
    marker_counts = marker_table["Evidence"].value_counts().reindex(
        ["Very Strong", "Strong", "Moderate", "Weak", "Weak/No Evidence"], fill_value=0
    )
    public_counts = legacy_table["Evidence"].value_counts().reindex(
        ["Very Strong", "Strong", "Moderate", "Weak", "Weak/No Evidence"], fill_value=0
    )

    summary_lines = [
        f"# {analysis_label} manuscript summary",
        "",
        f"- Additive markers tested: `{len(marker_table)}`.",
        f"- Public legacy-format Supplementary File 1 rows: `{len(legacy_table)}`.",
        f"- Public unique PBP/substitution rows: `{legacy_table[['PBP', 'Substitution']].drop_duplicates().shape[0]}`.",
        f"- Additive markers significant at analysis adjusted threshold `{additive_threshold}`: `{int(marker_table['mvLMM_significant_historical'].sum())}`.",
        f"- Additive markers with manuscript mvLMM evidence (`threshold + joint effect >= 1`): `{int(marker_table['mvLMM_evidence'].sum())}`.",
        f"- Additive markers significant at adjusted `0.001` for mvLMM/uvLMM display: `{int(marker_table['mvLMM_significant_p001'].sum())}`.",
        f"- Exact uvLMM marker-drug tests expected: `{len(marker_table) * 6}`.",
        f"- Exact uvLMM markers with >=1 significant drug: `{int((marker_table['uvLMM_exact_n_drugs_p001'] > 0).sum())}`.",
        f"- Markers with corrected epistasis evidence: `{int(marker_table['epistasis_evidence'].sum())}`.",
        f"- Markers with >=2 evidence streams: `{int((marker_table['evidence_total'] >= 2).sum())}`.",
        "",
        "## Evidence counts",
        "",
        "- Marker-level corrected table:",
        *[f"  - {label}: `{int(count)}`." for label, count in marker_counts.items()],
        "- Public legacy-format table:",
        *[f"  - {label}: `{int(count)}`." for label, count in public_counts.items()],
    ]

    if summary:
        summary_lines.extend(
            [
                "",
                "## Corrected epistasis merge summary",
                "",
                f"- Observed interactions: `{summary.get('observed_interactions', 'NA')}`.",
                f"- Permutation rows: `{summary.get('permutation_rows', 'NA')}`.",
                f"- Permutations: `{summary.get('permutations', 'NA')}`.",
                f"- Epistasis threshold: `{summary.get('epistasis_threshold', summary.get('historical_epistasis_threshold', 'NA'))}`.",
                f"- Epistasis threshold policy: `{summary.get('threshold_policy', 'NA')}`.",
                f"- P-value-threshold-only interactions: `{summary.get('pvalue_threshold_only_interactions', 'NA')}`.",
                f"- Supported interactions after lower-bound effect filter: `{summary.get('support_interactions', 'NA')}`.",
                f"- Observed non-ok fit rows: `{summary.get('observed_non_ok_fit_rows', 'NA')}`.",
                f"- Permutation non-ok fit rows: `{summary.get('permutation_non_ok_fit_rows', 'NA')}`.",
            ]
        )
        if summary.get("permutation_non_ok_fit_rows", 0):
            summary_lines.append(
                "- Non-ok permutation rows were retained as conservative null rows in the merged output; "
                "inspect `fit_status`/`fit_error` before changing manuscript claims."
            )

    if thresholds:
        summary_lines.extend(
            [
                "",
                "## Recomputed threshold summary",
                "",
                "```json",
                json.dumps(thresholds, indent=2, sort_keys=True),
                "```",
            ]
        )

    summary_lines.extend(
        [
            "",
            "## Output note",
            "",
            "- `Supplementary_File_1.csv` preserves the original manuscript/thesis component-level aggregation and column names.",
            "- `Supplementary_File_1_corrected_marker_level.csv` is the 170-marker marker-level summary table.",
            "- The public legacy column `Adj. p-value` is retained for compatibility and contains Galwey-adjusted additive mvLMM p-values, never raw `pv20` values.",
        ]
    )
    return summary_lines


def main() -> None:
    args = parse_args()
    output_dir = args.manuscript_dir if args.manuscript_dir is not None else args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    additive = load_additive(args.new_dir, args.additive_threshold, args.additive_galwey_meff, args.additive_dir)
    uv = load_uv_support(
        args.analysis_out,
        args.legacy_dir,
        args.new_dir,
        args.uv_threshold,
        args.additive_galwey_meff,
        args.uv_dir,
    )
    component = load_component_support(args.support_dir)
    epistasis_dir = args.epistasis_dir or args.analysis_out
    epistasis_path = epistasis_dir / "corrected_epistasis_marker_support.csv"
    epi = read_csv(epistasis_path) if epistasis_path.exists() else pd.DataFrame(columns=["marker", "num_sig_interactions"])

    marker_table = build_marker_level_table(additive, uv, component, epi)
    legacy_table = build_legacy_public_table(additive, uv, epi, args.support_dir, args.additive_threshold)

    write_marker_level_display(marker_table, output_dir)
    write_legacy_public(legacy_table, output_dir)

    summary_lines = build_summary(
        marker_table,
        legacy_table,
        epistasis_dir,
        args.additive_threshold,
        analysis_label=args.analysis_label,
        threshold_file=args.threshold_file,
    )
    (output_dir / "corrected_rerun_manuscript_summary.md").write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
