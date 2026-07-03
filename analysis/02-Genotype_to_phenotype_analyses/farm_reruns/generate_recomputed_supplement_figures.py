#!/usr/bin/env python3
"""Deprecated review plotter for recomputed Supplementary Figure S17/S18/S19.

This script generated exploratory review plots and should not be used for
publication-style Supplementary Figure S17/S18/S19 replacement panels.
Use generate_recomputed_supplement_figures_original_style.R instead, which
mirrors the original R/Rmd figure styles while reading recomputed data.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


PROJECT_ROOT = Path(__file__).resolve().parents[3]
FARM_ROOT = PROJECT_ROOT / "farm_outputs" / "recomputed_170_thresholds"
OUTDIR = FARM_ROOT / "manuscript_outputs" / "supplement_figures"

PBP_COLORS = {"1A": "#4DAF4A", "2B": "#377EB8", "2X": "#984EA3"}
EVIDENCE_RANK = {
    "Weak/No Evidence": 0,
    "Weak": 1,
    "Moderate": 2,
    "Strong": 3,
    "Very Strong": 4,
}


def read_thresholds() -> dict:
    return json.loads((FARM_ROOT / "thresholds" / "recomputed_thresholds.json").read_text())


def parse_pbp(marker: str) -> str:
    match = re.match(r"^PBP(1A|2B|2X)_", str(marker))
    if not match:
        raise ValueError(f"Cannot parse PBP from marker: {marker}")
    return match.group(1)


def parse_marker_label(marker: str) -> str:
    """Convert marker names to compact labels, preserving compound markers."""
    parts = re.split(r"_(?=PBP)", str(marker))
    labels: list[str] = []
    pbp = parse_pbp(parts[0])
    for part in parts:
        match = re.match(r"^PBP(1A|2B|2X)_([A-Z])(\d+)_([A-Z])$", part)
        if match:
            labels.append(f"{match.group(2)}{match.group(3)}{match.group(4)}")
        else:
            labels.append(part.replace(f"PBP{pbp}_", "").replace("_", ""))
    return f"{pbp} " + "/".join(labels)


def parse_position(substitution: str) -> int | None:
    match = re.match(r"^[A-Z](\d+)", str(substitution))
    return int(match.group(1)) if match else None


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def neg_log10(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").clip(lower=np.nextafter(0, 1))
    return -np.log10(numeric)


def load_epistasis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    epi = pd.read_csv(FARM_ROOT / "epistasis" / "merged" / "corrected_epistasis_p_values.csv")
    supported = pd.read_csv(
        FARM_ROOT / "epistasis" / "merged" / "corrected_epistasis_supported_interactions.csv"
    )
    marker_support = pd.read_csv(
        FARM_ROOT / "epistasis" / "merged" / "corrected_epistasis_marker_support.csv"
    )
    for frame in [epi, supported]:
        for col in [
            "passes_epistasis_threshold",
            "passes_epistasis_effect_filter",
            "epistasis_support",
        ]:
            frame[col] = bool_series(frame[col])
    marker_support["PBP"] = marker_support["marker"].map(parse_pbp)
    marker_support["label"] = marker_support["marker"].map(parse_marker_label)
    return epi, supported, marker_support


def marker_layout(marker_support: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Deterministic grouped circular layout for the supported-interaction network."""
    centers = {"1A": (-1.15, 0.55), "2B": (1.15, 0.55), "2X": (0.0, -0.95)}
    layout: dict[str, tuple[float, float]] = {}
    for pbp, group in marker_support.groupby("PBP", sort=False):
        group = group.sort_values(["num_sig_interactions", "marker"], ascending=[False, True])
        n = len(group)
        radius = 0.42 + 0.006 * n
        cx, cy = centers[pbp]
        for i, marker in enumerate(group["marker"]):
            angle = 2 * math.pi * i / max(n, 1)
            layout[marker] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))
    return layout


def generate_s17(thresholds: dict) -> dict:
    epi, supported, marker_support = load_epistasis()
    epi_threshold = thresholds["epistasis"]["epistasis_threshold"]

    epi = epi.copy()
    epi["plot_category"] = "Not significant"
    epi.loc[epi["passes_epistasis_threshold"], "plot_category"] = "P-threshold only"
    epi.loc[epi["epistasis_support"], "plot_category"] = "Supported"
    epi["neg_log10_adj_p"] = neg_log10(epi["pv20_adj_galwey"])
    epi["edge_weight"] = neg_log10(epi["pv20_adj_galwey"])

    OUTDIR.mkdir(parents=True, exist_ok=True)
    epi.to_csv(OUTDIR / "recomputed_S17_epistasis_scatter_source.csv", index=False)
    marker_support.to_csv(OUTDIR / "recomputed_S17_marker_support_source.csv", index=False)

    fig = plt.figure(figsize=(18, 5.8), constrained_layout=True)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.2, 1.05, 1.0])
    ax_scatter = fig.add_subplot(gs[0, 0])
    ax_network = fig.add_subplot(gs[0, 1])
    ax_hist = fig.add_subplot(gs[0, 2])

    scatter_styles = {
        "Not significant": {"color": "#D9D9D9", "alpha": 0.45, "zorder": 1},
        "P-threshold only": {"color": "#FCAE91", "alpha": 0.55, "zorder": 2},
        "Supported": {"color": "#CB181D", "alpha": 0.75, "zorder": 3},
    }
    for category, style in scatter_styles.items():
        data = epi[epi["plot_category"].eq(category)]
        sizes = 8 + 1.3 * np.sqrt(pd.to_numeric(data["n_present"], errors="coerce").fillna(0))
        ax_scatter.scatter(
            data["joint_effect_size"],
            data["neg_log10_adj_p"],
            s=sizes,
            linewidth=0,
            label=f"{category} (n={len(data):,})",
            **style,
        )
    ax_scatter.axvline(1, color="black", linewidth=1)
    ax_scatter.axhline(-math.log10(epi_threshold), color="black", linewidth=1)
    ax_scatter.set_xlabel("Joint effect size (MIC units)")
    ax_scatter.set_ylabel("-log10 adjusted p-value")
    ax_scatter.set_title("A. Epistasis LMM")
    ax_scatter.set_xlim(left=0)
    ax_scatter.legend(frameon=False, fontsize=8, loc="upper right")
    ax_scatter.grid(True, color="#DDDDDD", linewidth=0.8)

    layout = marker_layout(marker_support)
    support_lookup = marker_support.set_index("marker")
    max_edge = max(1.0, supported["pv20_adj_galwey"].map(lambda x: -math.log10(max(float(x), 1e-320))).max())
    for row in supported.itertuples(index=False):
        a = row.parent_a
        b = row.parent_b
        if a not in layout or b not in layout:
            continue
        x1, y1 = layout[a]
        x2, y2 = layout[b]
        edge_score = -math.log10(max(float(row.pv20_adj_galwey), 1e-320))
        ax_network.plot(
            [x1, x2],
            [y1, y2],
            color="#666666",
            alpha=0.035,
            linewidth=0.25 + 1.0 * edge_score / max_edge,
            zorder=1,
        )
    for pbp, group in marker_support.groupby("PBP"):
        xs = [layout[m][0] for m in group["marker"]]
        ys = [layout[m][1] for m in group["marker"]]
        sizes = 20 + 4.2 * group["num_sig_interactions"]
        ax_network.scatter(
            xs,
            ys,
            s=sizes,
            color=PBP_COLORS[pbp],
            edgecolor="black",
            linewidth=0.35,
            label=f"PBP{pbp}",
            zorder=2,
            alpha=0.95,
        )
    label_rows = (
        marker_support.sort_values(["PBP", "num_sig_interactions"], ascending=[True, False])
        .groupby("PBP", as_index=False)
        .head(1)
    )
    for row in label_rows.itertuples(index=False):
        x, y = layout[row.marker]
        ax_network.text(
            x,
            y + 0.08,
            row.label,
            fontsize=7,
            ha="center",
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.2},
        )
    ax_network.set_title("B. Supported interaction network")
    ax_network.set_aspect("equal")
    ax_network.axis("off")
    ax_network.legend(frameon=False, fontsize=8, loc="upper right")

    bins = np.arange(0, marker_support["num_sig_interactions"].max() + 5, 5)
    bottom = np.zeros(len(bins) - 1)
    for pbp in ["1A", "2B", "2X"]:
        values = marker_support.loc[marker_support["PBP"].eq(pbp), "num_sig_interactions"]
        counts, edges = np.histogram(values, bins=bins)
        ax_hist.bar(
            edges[:-1],
            counts,
            width=np.diff(edges),
            align="edge",
            bottom=bottom,
            color=PBP_COLORS[pbp],
            edgecolor="black",
            linewidth=0.5,
            label=f"PBP{pbp}",
        )
        bottom += counts
    ax_hist.set_xlabel("Supported interactions per marker")
    ax_hist.set_ylabel("Marker count")
    ax_hist.set_title("C. Interaction-degree histogram")
    ax_hist.legend(frameon=False, fontsize=8)
    ax_hist.grid(True, axis="y", color="#DDDDDD", linewidth=0.8)

    fig.suptitle("Recomputed Supplementary Figure S17: epistatic interaction LMM", fontsize=14)
    out = OUTDIR / "recomputed_S17_epistasis.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)

    return {
        "figure": str(out),
        "tested_interactions": int(len(epi)),
        "p_threshold_only_interactions": int(epi["passes_epistasis_threshold"].sum()),
        "supported_interactions": int(epi["epistasis_support"].sum()),
        "markers_with_support": int(len(marker_support)),
        "markers_gt_40_supported_interactions": int((marker_support["num_sig_interactions"] > 40).sum()),
        "epistasis_threshold": epi_threshold,
    }


def pvalue_to_float(value: object) -> float:
    text = str(value).strip()
    if text.startswith("<"):
        text = text[1:]
    try:
        return float(text)
    except ValueError:
        return float("nan")


def category(value: float) -> str:
    if value < -1:
        return "-"
    if value > 1:
        return "+"
    return "/"


def additive_public_effects(thresholds: dict) -> pd.DataFrame:
    sf1 = pd.read_csv(FARM_ROOT / "manuscript_outputs" / "Supplementary_File_1.csv")
    additive_threshold = thresholds["additive"]["adjusted_threshold"]
    data = sf1.copy()
    data["P-value numeric"] = data["Adj. p-value"].map(pvalue_to_float)
    data["Effect Size Axis 1"] = pd.to_numeric(data["Effect Size Axis 1"], errors="coerce")
    data["Effect Size Axis 2"] = pd.to_numeric(data["Effect Size Axis 2"], errors="coerce")
    data = data[
        data["P-value numeric"].le(additive_threshold)
        & data["Effect Size Axis 1"].notna()
        & data["Effect Size Axis 2"].notna()
    ].copy()
    data["source"] = "Additive mvLMM"
    data["effect_x"] = data["Effect Size Axis 1"]
    data["effect_y"] = data["Effect Size Axis 2"]
    data["label"] = data["PBP"].astype(str) + " " + data["Substitution"].astype(str)
    return data


def epistasis_p_threshold_effects(thresholds: dict) -> pd.DataFrame:
    epi, _, _ = load_epistasis()
    epi_threshold = thresholds["epistasis"]["epistasis_threshold"]
    data = epi[pd.to_numeric(epi["pv20_adj_galwey"], errors="coerce").le(epi_threshold)].copy()
    data["source"] = "Epistatic LMM"
    data["effect_x"] = pd.to_numeric(data["effsize_env1_D1"], errors="coerce")
    data["effect_y"] = pd.to_numeric(data["effsize_env1_D2"], errors="coerce")
    data["label"] = data["interaction"]
    return data


def generate_s18(thresholds: dict) -> dict:
    additive = additive_public_effects(thresholds)
    epistasis = epistasis_p_threshold_effects(thresholds)
    additive_source = additive[
        [
            "Evidence",
            "PBP",
            "Substitution",
            "β Joint",
            "Adj. p-value",
            "Effect Size Axis 1",
            "Effect Size Axis 2",
            "source",
            "effect_x",
            "effect_y",
            "label",
        ]
    ].copy()
    epistasis_source = epistasis[
        [
            "test",
            "interaction",
            "parent_a",
            "parent_b",
            "pv20_adj_galwey",
            "joint_effect_size",
            "source",
            "effect_x",
            "effect_y",
            "label",
        ]
    ].copy()
    additive_source.to_csv(OUTDIR / "recomputed_S18_additive_source.csv", index=False)
    epistasis_source.to_csv(OUTDIR / "recomputed_S18_epistasis_source.csv", index=False)

    combined = pd.concat(
        [
            additive[["source", "effect_x", "effect_y", "label"]],
            epistasis[["source", "effect_x", "effect_y", "label"]],
        ],
        ignore_index=True,
    )
    combined["axis1_category"] = combined["effect_x"].map(category)
    combined["axis2_category"] = combined["effect_y"].map(category)
    counts = pd.crosstab(
        pd.Categorical(combined["axis2_category"], categories=["-", "/", "+"]),
        pd.Categorical(combined["axis1_category"], categories=["-", "/", "+"]),
        dropna=False,
    )
    counts.index.name = "axis2_category"
    counts.columns.name = "axis1_category"
    counts.to_csv(OUTDIR / "recomputed_S18_axis_category_counts.csv")

    summary_rows = []
    total = len(combined)
    for axis2 in ["-", "/", "+"]:
        for axis1 in ["-", "/", "+"]:
            n = int(counts.loc[axis2, axis1])
            summary_rows.append(
                {
                    "axis2_category": axis2,
                    "axis1_category": axis1,
                    "count": n,
                    "percent": round(100 * n / total, 2) if total else 0,
                }
            )
    pd.DataFrame(summary_rows).to_csv(OUTDIR / "recomputed_S18_axis_category_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    ax.scatter(
        epistasis["effect_x"],
        epistasis["effect_y"],
        s=10,
        color="#CB181D",
        alpha=0.18,
        linewidth=0,
        label=f"Epistatic LMM (n={len(epistasis):,})",
    )
    ax.scatter(
        additive["effect_x"],
        additive["effect_y"],
        s=34,
        facecolor="#FB6A4A",
        edgecolor="black",
        linewidth=0.35,
        alpha=0.85,
        label=f"Additive mvLMM (n={len(additive):,})",
    )
    for value in [-1, 1]:
        ax.axvline(value, color="black", linewidth=0.8)
        ax.axhline(value, color="black", linewidth=0.8)
    ax.axvline(0, color="#888888", linewidth=0.5, linestyle=":")
    ax.axhline(0, color="#888888", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Effect size - map axis 1")
    ax.set_ylabel("Effect size - map axis 2")
    ax.set_title("Recomputed Supplementary Figure S18")
    ax.set_xlim(-4.5, 5.0)
    ax.set_ylim(-4.0, 3.4)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#DDDDDD", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left", fontsize=8)

    x_positions = {"-": -3.65, "/": 0.0, "+": 3.55}
    y_positions = {"-": -3.25, "/": 0.0, "+": 2.65}
    for row in summary_rows:
        label = f"{row['count']:,}\n({row['percent']:.2f}%)"
        ax.text(
            x_positions[row["axis1_category"]],
            y_positions[row["axis2_category"]],
            label,
            ha="center",
            va="center",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.72, "pad": 2},
        )

    out = OUTDIR / "recomputed_S18_effect_size_axes.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)

    return {
        "figure": str(out),
        "additive_component_expanded_rows": int(len(additive)),
        "epistasis_p_threshold_rows": int(len(epistasis)),
        "total_plotted_effects": int(total),
        "axis_category_counts_csv": str(OUTDIR / "recomputed_S18_axis_category_counts.csv"),
    }


def load_overlap_sets() -> tuple[dict[str, set[str]], pd.DataFrame]:
    sf1 = pd.read_csv(FARM_ROOT / "manuscript_outputs" / "Supplementary_File_1.csv")
    sf1["Position"] = sf1["Substitution"].map(parse_position)
    sf1 = sf1[sf1["Position"].notna()].copy()
    sf1["PBP_position"] = sf1["PBP"].astype(str) + "_" + sf1["Position"].astype(int).astype(str)
    sf1["rank"] = sf1["Evidence"].map(EVIDENCE_RANK)
    strongest = sf1.sort_values("rank").groupby("PBP_position", as_index=False).tail(1)

    prior = pd.read_csv(
        PROJECT_ROOT
        / "AMRC-repo-files"
        / "Streptococcus pneumoniae analysis"
        / "CDC_GWAS_overlap_TPD.csv",
        encoding="utf-8-sig",
    )
    prior["PBP"] = prior["PBP"].replace({"1a": "1A", "2b": "2B", "2x": "2X"})
    prior["PBP_position"] = (
        prior["PBP"].astype(str) + "_" + prior["Position"].astype(int).astype(str)
    )
    for col in ["GWAS", "CDC", "Laboratory"]:
        prior[col] = prior[col].astype(str).str.lower().eq("yes")

    sets = {
        "Cartography Weak+": set(strongest.loc[strongest["rank"].ge(1), "PBP_position"]),
        "Cartography Moderate+": set(strongest.loc[strongest["rank"].ge(2), "PBP_position"]),
        "GWAS": set(prior.loc[prior["GWAS"], "PBP_position"]),
        "CDC": set(prior.loc[prior["CDC"], "PBP_position"]),
        "Laboratory": set(prior.loc[prior["Laboratory"], "PBP_position"]),
    }
    return sets, prior


def intersection_frame(sets: dict[str, set[str]]) -> pd.DataFrame:
    set_names = list(sets)
    positions = sorted(set().union(*sets.values()))
    rows = []
    for position in positions:
        row = {"PBP_position": position}
        for name in set_names:
            row[name] = position in sets[name]
        rows.append(row)
    membership = pd.DataFrame(rows)
    group = membership.groupby(set_names, dropna=False).size().reset_index(name="count")
    group = group[group["count"].gt(0)].sort_values(["count"] + set_names, ascending=[False] + [False] * len(set_names))
    group["combination"] = [
        " + ".join(name for name in set_names if row[name]) or "none"
        for row in group.to_dict("records")
    ]
    return group


def generate_s19() -> dict:
    sets, prior = load_overlap_sets()
    intersections = intersection_frame(sets)
    intersections.to_csv(OUTDIR / "recomputed_S19_overlap_intersections.csv", index=False)

    set_sizes = pd.DataFrame(
        [{"set": name, "size": len(values)} for name, values in sets.items()]
    )
    set_sizes.to_csv(OUTDIR / "recomputed_S19_overlap_set_sizes.csv", index=False)

    any_previous = sets["GWAS"] | sets["CDC"] | sets["Laboratory"]
    lab_only = sets["Laboratory"] - sets["GWAS"] - sets["CDC"]
    gwas_only = sets["GWAS"] - sets["Laboratory"] - sets["CDC"]
    summary = {
        "cartography_weak_plus_positions": len(sets["Cartography Weak+"]),
        "cartography_moderate_plus_positions": len(sets["Cartography Moderate+"]),
        "previous_any_source_positions": len(any_previous),
        "weak_plus_overlap_any_previous": len(sets["Cartography Weak+"] & any_previous),
        "weak_plus_not_previous": len(sets["Cartography Weak+"] - any_previous),
        "weak_plus_lab_only_recovered": len(sets["Cartography Weak+"] & lab_only),
        "previous_not_recovered_by_weak_plus": len(any_previous - sets["Cartography Weak+"]),
        "previous_not_recovered_gwas_only": len(gwas_only - sets["Cartography Weak+"]),
        "previous_not_recovered_lab_only": len(lab_only - sets["Cartography Weak+"]),
        "moderate_plus_overlap_any_previous": len(sets["Cartography Moderate+"] & any_previous),
        "moderate_plus_not_previous": len(sets["Cartography Moderate+"] - any_previous),
        "moderate_plus_lab_only_recovered": len(sets["Cartography Moderate+"] & lab_only),
        "previous_not_recovered_by_moderate_plus": len(any_previous - sets["Cartography Moderate+"]),
    }
    pd.DataFrame([summary]).to_csv(OUTDIR / "recomputed_S19_overlap_summary.csv", index=False)

    set_names = list(sets)
    top = intersections.head(28).reset_index(drop=True)
    fig = plt.figure(figsize=(13.5, 8.2), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.2, 4.3], height_ratios=[3.2, 1.25])
    ax_sizes = fig.add_subplot(gs[:, 0])
    ax_bars = fig.add_subplot(gs[0, 1])
    ax_matrix = fig.add_subplot(gs[1, 1])

    size_colors = {
        "Cartography Weak+": "#FDE0DD",
        "Cartography Moderate+": "#CB181D",
        "GWAS": "#4DAF4A",
        "CDC": "#FFD92F",
        "Laboratory": "#377EB8",
    }
    size_df = set_sizes.set_index("set").loc[set_names].reset_index()
    y = np.arange(len(size_df))
    size_bars = ax_sizes.barh(
        y,
        size_df["size"],
        color=[size_colors[name] for name in size_df["set"]],
        edgecolor="black",
        linewidth=0.6,
    )
    ax_sizes.set_yticks(y)
    ax_sizes.set_yticklabels(size_df["set"], fontsize=8)
    ax_sizes.invert_yaxis()
    ax_sizes.set_xlabel("Set size")
    ax_sizes.set_title("Sets")
    ax_sizes.grid(True, axis="x", color="#DDDDDD", linewidth=0.8)
    ax_sizes.bar_label(size_bars, labels=[str(v) for v in size_df["size"]], padding=3, fontsize=8)

    bar_colors = []
    for row in top.to_dict("records"):
        if row["Cartography Moderate+"]:
            bar_colors.append("#CB181D")
        elif row["Cartography Weak+"]:
            bar_colors.append("#FCAEAE")
        else:
            bar_colors.append("#999999")
    bars = ax_bars.bar(
        np.arange(len(top)), top["count"], color=bar_colors, edgecolor="black", linewidth=0.5
    )
    ax_bars.bar_label(bars, labels=[str(v) for v in top["count"]], padding=2, fontsize=7)
    ax_bars.set_ylabel("Position count")
    ax_bars.set_title("Recomputed Supplementary Figure S19: overlap with previous work")
    ax_bars.set_xticks([])
    ax_bars.grid(True, axis="y", color="#DDDDDD", linewidth=0.8)

    for i, row in top.iterrows():
        active = [j for j, name in enumerate(set_names) if row[name]]
        if active:
            ax_matrix.plot([i, i], [min(active), max(active)], color="#555555", linewidth=1.0, zorder=1)
        for j, name in enumerate(set_names):
            active_here = bool(row[name])
            ax_matrix.scatter(
                i,
                j,
                s=42 if active_here else 18,
                color="#222222" if active_here else "#DDDDDD",
                zorder=2,
            )
    ax_matrix.set_yticks(np.arange(len(set_names)))
    ax_matrix.set_yticklabels(set_names, fontsize=8)
    ax_matrix.set_ylim(len(set_names) - 0.5, -0.5)
    ax_matrix.set_xlim(-0.7, len(top) - 0.3)
    ax_matrix.set_xlabel("Intersection combinations, sorted by size")
    ax_matrix.grid(False)
    for spine in ["top", "right", "left"]:
        ax_matrix.spines[spine].set_visible(False)

    out = OUTDIR / "recomputed_S19_overlap_upset.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)

    return {"figure": str(out), **summary}


def write_summary(s17: dict, s18: dict, s19: dict) -> None:
    lines = [
        "# Recomputed Supplementary Figure Review Outputs",
        "",
        "These files were generated from `farm_outputs/recomputed_170_thresholds/` and do not modify manuscript or supplement sources.",
        "",
        "## S17 Epistasis",
        "",
        f"- Figure: `{Path(s17['figure']).relative_to(PROJECT_ROOT)}`",
        f"- Tested interactions: `{s17['tested_interactions']}`",
        f"- P-threshold-only interactions: `{s17['p_threshold_only_interactions']}`",
        f"- Supported interactions after lower-bound effect filter: `{s17['supported_interactions']}`",
        f"- Markers with epistasis support: `{s17['markers_with_support']}`",
        f"- Markers with >40 supported interactions: `{s17['markers_gt_40_supported_interactions']}`",
        f"- Adjusted epistasis threshold: `{s17['epistasis_threshold']}`",
        "",
        "## S18 Effect Directions",
        "",
        f"- Figure: `{Path(s18['figure']).relative_to(PROJECT_ROOT)}`",
        f"- Additive component-expanded public rows plotted: `{s18['additive_component_expanded_rows']}`",
        f"- Epistasis p-threshold rows plotted: `{s18['epistasis_p_threshold_rows']}`",
        f"- Total plotted effects: `{s18['total_plotted_effects']}`",
        "- Additive points use the component-expanded public Supplementary File 1 frame.",
        "- Epistasis points use the p-threshold-only interaction frame, matching the original S18 logic before the lower-bound support filter.",
        "",
        "## S19 Previous-Work Overlap",
        "",
        f"- Figure: `{Path(s19['figure']).relative_to(PROJECT_ROOT)}`",
        f"- Cartography Weak+ positions: `{s19['cartography_weak_plus_positions']}`",
        f"- Cartography Moderate+ positions: `{s19['cartography_moderate_plus_positions']}`",
        f"- Previous-source positions: `{s19['previous_any_source_positions']}`",
        f"- Weak+ overlap with any previous source: `{s19['weak_plus_overlap_any_previous']}`",
        f"- Weak+ positions not highlighted by previous sources: `{s19['weak_plus_not_previous']}`",
        f"- Laboratory-only positions recovered by Weak+ cartography: `{s19['weak_plus_lab_only_recovered']}`",
        f"- Previous-source positions not recovered by Weak+ cartography: `{s19['previous_not_recovered_by_weak_plus']}`",
        f"- Previous-source not recovered split: `{s19['previous_not_recovered_gwas_only']}` GWAS-only, `{s19['previous_not_recovered_lab_only']}` laboratory-only",
        f"- Moderate+ overlap with any previous source: `{s19['moderate_plus_overlap_any_previous']}`",
        f"- Moderate+ positions not highlighted by previous sources: `{s19['moderate_plus_not_previous']}`",
        "",
        "## Remaining Literature-Curation Work",
        "",
        "- The overlap file is position-level. It does not distinguish exact-substitution matches from same-position matches.",
        "- Manual curation is still needed for exact-substitution support, same-motif/nearby-region support, and in vitro mechanism notes.",
        "- The prior-overlap input is currently read through the `AMRC-repo-files` symlink and should be provenance-stamped before release.",
        "",
    ]
    (OUTDIR / "recomputed_supplement_figures_summary.md").write_text("\n".join(lines))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    thresholds = read_thresholds()
    s17 = generate_s17(thresholds)
    s18 = generate_s18(thresholds)
    s19 = generate_s19()
    write_summary(s17, s18, s19)
    print(f"Wrote recomputed supplement figure outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
