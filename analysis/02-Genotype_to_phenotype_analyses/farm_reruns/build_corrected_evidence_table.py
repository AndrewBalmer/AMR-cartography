#!/usr/bin/env python3
"""Build manuscript-facing evidence tables after corrected farm reruns."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import DEFAULT_NEW_THRESHOLD, DEFAULT_UV_THRESHOLD, marker_order_from_effects, read_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--new-dir", required=True, type=Path)
    parser.add_argument("--support-dir", required=True, type=Path)
    parser.add_argument("--analysis-out", required=True, type=Path)
    parser.add_argument("--manuscript-dir", default=Path("manuscript"), type=Path)
    parser.add_argument("--additive-threshold", default=DEFAULT_NEW_THRESHOLD, type=float)
    parser.add_argument("--uv-threshold", default=DEFAULT_UV_THRESHOLD, type=float)
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
        for prot, aa1, pos, aa2 in re.findall(r"(PBP[12][ABX])_([A-Z])([0-9]+)_([A-Z])", marker)
    ]


def compact_label(marker: str) -> str:
    parts = parse_marker(marker)
    return "/".join(str(part["aa"]) for part in parts) if parts else marker


def primary_pbp(marker: str) -> str:
    parts = parse_marker(marker)
    proteins = list(dict.fromkeys(str(part["protein"]) for part in parts))
    return proteins[0] if len(proteins) == 1 else ("Mixed" if proteins else "")


def load_additive(new_dir: Path, threshold: float) -> pd.DataFrame:
    pvals = read_csv(new_dir / "mvLMM_p_values_normal_pneumo_low_freq_vars.csv")
    unnamed = [col for col in pvals.columns if str(col).startswith("Unnamed")]
    if unnamed:
        pvals = pvals.rename(columns={unnamed[0]: "mvLMM_result_index"})
    order = marker_order_from_effects(new_dir / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv")
    if len(order) != len(pvals):
        raise ValueError("Additive mvLMM p-values/effects do not align")
    pvals["marker"] = order
    effects = read_csv(new_dir / "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv")
    cand = effects[effects["effect_type"].eq("candidate")][["effect_name", "env", "effsize"]]
    wide = cand.pivot(index="effect_name", columns="env", values="effsize").reset_index()
    wide = wide.rename(columns={"effect_name": "marker", "env1_D1": "effect_axis1", "env1_D2": "effect_axis2"})
    out = pvals.merge(wide, on="marker", how="left")
    out["joint_effect_size"] = np.sqrt(out["effect_axis1"] ** 2 + out["effect_axis2"] ** 2)
    out["mvLMM_significant_perm"] = out["pv20"] < threshold
    out["mvLMM_evidence"] = out["mvLMM_significant_perm"] & (out["joint_effect_size"] >= 1)
    out["PBP"] = out["marker"].map(primary_pbp)
    out["Substitution"] = out["marker"].map(compact_label)
    out["constituents"] = out["marker"].map(
        lambda marker: ";".join(f"{part['protein']}:{part['aa']}" for part in parse_marker(marker))
    )
    return out


def load_uv_support(analysis_out: Path, legacy_dir: Path, new_dir: Path, threshold: float) -> pd.DataFrame:
    support_path = analysis_out / "uniLMM_exact_170_marker_support.csv"
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
    rows = []
    for marker, group in uv.groupby("marker"):
        sig = sorted(group.loc[group["pv20"] < threshold, "drug"].astype(str).unique())
        rows.append(
            {
                "marker": marker,
                "uvLMM_exact_n_drugs_p001": len(sig),
                "uvLMM_exact_drugs_p001": ";".join(sig),
                "uvLMM_exact_source": ";".join(sorted(group["source"].unique())),
            }
        )
    return pd.DataFrame(rows)


def load_component_support(support_dir: Path) -> pd.DataFrame:
    single = read_csv(support_dir / "Single_subs_all_S.pneumo.csv")
    single["PBP"] = single["PBP"].astype(str).str.strip()
    single["aa"] = (
        single["Substitution_of_interest"]
        .astype(str)
        .str.replace(r"\([0-9]+\)", "", regex=True)
        .str.replace(r"\s+", "", regex=True)
    )
    single["single_key"] = single["PBP"] + ":" + single["aa"]
    single["median_phenotypic_distance"] = pd.to_numeric(single["median_phenotypic_distance"], errors="coerce")
    single_hits = set(single.loc[single["median_phenotypic_distance"] >= 1, "single_key"])

    cluster_files = sorted(support_dir.glob("Sig_AA_subs_Cluster_gen*_S.pneumo.csv"))
    cluster = pd.concat([pd.read_csv(path) for path in cluster_files], ignore_index=True)
    cluster_hits = set(cluster.loc[cluster["Confidence"].eq("Strong"), "Location"].astype(str))

    prior = read_csv(support_dir / "CDC_GWAS_overlap_TPD.csv")
    prior["PBP"] = prior["PBP"].replace({"1a": "PBP1A", "2b": "PBP2B", "2x": "PBP2X"})
    prior["pos_key"] = prior.apply(lambda row: f"{row['PBP']}_{int(row['Position'])}", axis=1)
    for col in ["GWAS", "CDC", "Laboratory"]:
        prior[col] = prior[col].astype(str).str.lower().eq("yes")
    prior_hits = set(prior.loc[prior[["GWAS", "CDC", "Laboratory"]].any(axis=1), "pos_key"])

    return pd.DataFrame(
        [
            {
                "support_type": "sets",
                "single_hits": single_hits,
                "cluster_hits": cluster_hits,
                "prior_hits": prior_hits,
            }
        ]
    )


def evidence_category(total: int) -> str:
    return {0: "Weak/No Evidence", 1: "Weak", 2: "Moderate", 3: "Strong", 4: "Very Strong"}.get(total, "Very Strong")


def main() -> None:
    args = parse_args()
    args.manuscript_dir.mkdir(parents=True, exist_ok=True)
    additive = load_additive(args.new_dir, args.additive_threshold)
    uv = load_uv_support(args.analysis_out, args.legacy_dir, args.new_dir, args.uv_threshold)
    component = load_component_support(args.support_dir).iloc[0]
    epistasis_path = args.analysis_out / "corrected_epistasis_marker_support.csv"
    epi = read_csv(epistasis_path) if epistasis_path.exists() else pd.DataFrame(columns=["marker", "num_sig_interactions"])

    table = additive.merge(uv, on="marker", how="left").merge(epi, on="marker", how="left")
    table["uvLMM_exact_n_drugs_p001"] = pd.to_numeric(table["uvLMM_exact_n_drugs_p001"], errors="coerce").fillna(0).astype(int)
    table["num_sig_interactions"] = pd.to_numeric(table["num_sig_interactions"], errors="coerce").fillna(0).astype(int)
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
        lambda row: f"{'Yes' if row['mvLMM_significant_perm'] else 'No'}/"
        f"{'Yes' if row['uvLMM_exact_n_drugs_p001'] > 0 else 'No'} ({row['uvLMM_exact_n_drugs_p001']})",
        axis=1,
    )
    table = table.sort_values(["evidence_total", "joint_effect_size"], ascending=[False, False])
    table.to_csv(args.manuscript_dir / "Supplementary_File_1_corrected_marker_level.csv", index=False)

    display = pd.DataFrame(
        {
            "Evidence": table["Evidence"],
            "PBP": table["PBP"].replace({"PBP1A": "1A", "PBP2B": "2B", "PBP2X": "2X"}),
            "Substitution": table["Substitution"],
            "Marker": table["marker"],
            "Constituent substitutions": table["constituents"],
            "mvLMM p-value": table["pv20"].map(lambda value: "<1e-16" if value < 1e-16 else f"{value:.3e}"),
            "beta Joint": table["joint_effect_size"].map(lambda value: f"{value:.3f}"),
            "Effect Size Axis 1": table["effect_axis1"].map(lambda value: f"{value:.3f}"),
            "Effect Size Axis 2": table["effect_axis2"].map(lambda value: f"{value:.3f}"),
            "No. Sig. interactions": table["num_sig_interactions"].replace(0, "-"),
            "total": table["evidence_total"],
            "Sig. mvLMM/uvLMM (No. drugs)": table["Sig. mvLMM/uvLMM (No. drugs)"],
            "Prior GWAS/CDC/lab position": table["prior_any_position"],
        }
    )
    display.to_csv(args.manuscript_dir / "Supplementary_File_1.csv", index=False)

    audit = [
        "# Corrected rerun manuscript audit",
        "",
        f"- Additive markers tested: `{len(table)}`.",
        f"- Additive markers significant at `{args.additive_threshold}`: `{int(table['mvLMM_significant_perm'].sum())}`.",
        f"- Exact uvLMM marker-drug tests expected: `1020`.",
        f"- Exact uvLMM markers with >=1 significant drug: `{int((table['uvLMM_exact_n_drugs_p001'] > 0).sum())}`.",
        f"- Markers with corrected epistasis evidence: `{int(table['epistasis_evidence'].sum())}`.",
        f"- Markers with >=2 evidence streams: `{int((table['evidence_total'] >= 2).sum())}`.",
    ]
    (args.manuscript_dir / "corrected_rerun_manuscript_audit.md").write_text("\n".join(audit) + "\n")
    print("\n".join(audit))


if __name__ == "__main__":
    main()
