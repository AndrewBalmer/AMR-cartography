#!/usr/bin/env python3
"""Annotate recomputed Supplementary File 1 with PBP motif/structure context.

This is a post hoc biological annotation step. It does not alter evidence
scores, p-values, or manuscript-facing source files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SF1 = (
    PROJECT_ROOT
    / "analysis_outputs"
    / "recomputed_170_thresholds"
    / "manuscript_outputs"
    / "Supplementary_File_1.csv"
)
DEFAULT_PRIOR = (
    PROJECT_ROOT
    / "AMRC-repo-files"
    / "Streptococcus pneumoniae analysis"
    / "CDC_GWAS_overlap_TPD.csv"
)
DEFAULT_OUTDIR = (
    PROJECT_ROOT
    / "analysis_outputs"
    / "recomputed_170_thresholds"
    / "manuscript_outputs"
    / "structure_context"
)


@dataclass(frozen=True)
class StructureConfig:
    pdb_id: str
    chain_id: str
    motifs: dict[str, range]
    transpeptidase_range: range
    source_url: str
    note: str


STRUCTURES: dict[str, StructureConfig] = {
    "1A": StructureConfig(
        pdb_id="2C6W",
        chain_id="B",
        motifs={
            "S370TMK": range(370, 374),
            "S466SN": range(466, 469),
            "K557TG": range(557, 560),
        },
        transpeptidase_range=range(371, 648),
        source_url="https://www.rcsb.org/structure/2C6W",
        note="PBP1A reference chain B; direct residue numbering used.",
    ),
    "2B": StructureConfig(
        pdb_id="2WAF",
        chain_id="A",
        motifs={
            "S386VVK": range(386, 390),
            "S443SN": range(443, 446),
            "K615TG": range(615, 618),
        },
        transpeptidase_range=range(379, 657),
        source_url="https://www.rcsb.org/structure/2WAF",
        note="PBP2B reference chain A; direct residue numbering used.",
    ),
    "2X": StructureConfig(
        pdb_id="5OIZ",
        chain_id="A",
        motifs={
            "S337TMK": range(337, 341),
            "S395SN": range(395, 398),
            "K547SG": range(547, 550),
        },
        transpeptidase_range=range(229, 588),
        source_url="https://www.rcsb.org/structure/5OIZ",
        note="PBP2X reference chain A; direct residue numbering used.",
    ),
}


EVIDENCE_RANK = {
    "Weak/No Evidence": 0,
    "Weak": 1,
    "Moderate": 2,
    "Strong": 3,
    "Very Strong": 4,
}
PLOT_CATEGORIES = [
    ("Very Strong", "very_strong_count", "#7F0000"),
    ("Strong", "strong_count", "#CB181D"),
    ("Moderate", "moderate_count", "#FB6A4A"),
    ("Weak/No evidence", "weak_no_evidence_count", "#FDE0DD"),
]
INVARIANT_COLOR = "#EEEEEE"

SUB_RE = re.compile(r"^[A-Z](?P<position>[0-9]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sf1", type=Path, default=DEFAULT_SF1)
    parser.add_argument("--prior-overlap", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Optional recomputed-analysis output root. When supplied, the "
        "default --sf1 and --outdir paths are resolved under this directory.",
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=None,
        help="Directory of committed/cached PDB files for offline, reproducible runs. "
        "Defaults to <outdir>/pdb. A file named <accession>.pdb here is used verbatim "
        "instead of downloading from RCSB.",
    )
    parser.add_argument("--bin-width", type=int, default=1)
    return parser.parse_args()


def _looks_like_pdb(data: bytes) -> bool:
    head = data[:4000].decode("latin-1", errors="ignore")
    return ("ATOM  " in head) or ("HEADER" in head) or ("CRYST1" in head)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_pdb(pdb_id: str, pdb_dir: Path) -> Path:
    """Return a local PDB file, preferring a committed/cached copy for reproducibility.

    A copy already present in ``pdb_dir`` (which can be committed for offline
    reproduction) is used verbatim. Otherwise the structure is downloaded from
    RCSB by accession and validated, so an HTML error page or truncated response
    is never silently accepted.
    """
    pdb_dir.mkdir(parents=True, exist_ok=True)
    path = pdb_dir / f"{pdb_id}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
    except Exception as exc:  # offline / RCSB unreachable
        raise RuntimeError(
            f"Could not download {pdb_id} from {url} ({exc}). For offline or "
            f"reproducible runs, place {pdb_id}.pdb in a directory and pass "
            f"--pdb-dir."
        ) from exc
    if not _looks_like_pdb(data):
        raise RuntimeError(
            f"Downloaded {pdb_id} from {url} but it does not look like a PDB file "
            f"({len(data)} bytes); refusing to use it."
        )
    path.write_bytes(data)
    return path


def parse_ca_coordinates(pdb_path: Path, chain_id: str) -> dict[int, tuple[float, float, float]]:
    coords: dict[int, tuple[float, float, float]] = {}
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[21].strip() != chain_id:
                continue
            try:
                residue = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords[residue] = (x, y, z)
    return coords


def euclidean(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def parse_position(substitution: str) -> int | None:
    match = SUB_RE.match(str(substitution))
    if not match:
        return None
    return int(match.group("position"))


def motif_residues(config: StructureConfig) -> list[tuple[str, int]]:
    residues: list[tuple[str, int]] = []
    for motif, positions in config.motifs.items():
        for position in positions:
            residues.append((motif, position))
    return residues


def nearest_by_sequence(pbp: str, position: int | None) -> tuple[str, int | None]:
    if position is None:
        return ("not parsed", None)
    best: tuple[int, str] | None = None
    for motif, motif_position in motif_residues(STRUCTURES[pbp]):
        distance = abs(position - motif_position)
        candidate = (distance, motif)
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return (best[1], best[0])


def motif_interval(pbp: str, position: int | None) -> str:
    if position is None:
        return "not parsed"
    config = STRUCTURES[pbp]
    motifs = list(config.motifs.items())
    first_name, first_range = motifs[0]
    second_name, second_range = motifs[1]
    third_name, third_range = motifs[2]

    if any(position in residues for _, residues in motifs):
        return "within active-site motif"
    if position < first_range.start:
        return f"before {first_name}"
    if first_range.stop <= position < second_range.start:
        return "SXXK-SXN interval"
    if second_range.stop <= position < third_range.start:
        return "SXN-KTG/KSG interval"
    if position >= third_range.stop:
        return "after KTG/KSG"
    return "between motif residues"


def nearest_by_structure(
    pbp: str,
    position: int | None,
    coords_by_pbp: dict[str, dict[int, tuple[float, float, float]]],
) -> tuple[str, float | str]:
    if position is None:
        return ("not parsed", "not resolved")
    coords = coords_by_pbp[pbp]
    if position not in coords:
        return ("not resolved", "not resolved")
    motif_candidates = [
        (motif, motif_position)
        for motif, motif_position in motif_residues(STRUCTURES[pbp])
        if motif_position in coords
    ]
    if not motif_candidates:
        return ("not resolved", "not resolved")
    best: tuple[float, str, int] | None = None
    for motif, motif_position in motif_candidates:
        distance = euclidean(coords[position], coords[motif_position])
        candidate = (distance, motif, motif_position)
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return (best[1], round(best[0], 3))


def load_prior_overlap(path: Path) -> dict[tuple[str, int], tuple[str, str]]:
    prior = pd.read_csv(path, encoding="utf-8-sig")
    prior["PBP"] = prior["PBP"].replace({"1a": "1A", "2b": "2B", "2x": "2X"})
    source_labels = {
        "GWAS": "Chewapreecha et al. GWAS",
        "CDC": "Li et al. CDC study",
        "Laboratory": "Laboratory (see Supplementary S3.2)",
    }
    mapping: dict[tuple[str, int], tuple[str, str]] = {}
    for row in prior.itertuples(index=False):
        sources = []
        for source in ["GWAS", "CDC", "Laboratory"]:
            if str(getattr(row, source)).strip().lower() == "yes":
                sources.append(source_labels[source])
        category = "same-position overlap" if sources else "not previously highlighted"
        source_text = "; ".join(sources) if sources else "none in existing overlap file"
        mapping[(row.PBP, int(row.Position))] = (category, source_text)
    return mapping


def annotate_sf1(
    sf1: Path,
    prior_overlap: Path,
    coords_by_pbp: dict[str, dict[int, tuple[float, float, float]]],
) -> pd.DataFrame:
    df = pd.read_csv(sf1)
    prior = load_prior_overlap(prior_overlap)

    parsed_positions: list[int | None] = []
    motif_intervals: list[str] = []
    nearest_sequence_motifs: list[str] = []
    sequence_distances: list[int | None] = []
    nearest_3d_motifs: list[str] = []
    distances_3d: list[float | str] = []
    prior_categories: list[str] = []
    prior_sources: list[str] = []

    for row in df.itertuples(index=False):
        pbp = str(row.PBP)
        substitution = str(row.Substitution)
        position = parse_position(substitution)
        parsed_positions.append(position)
        motif_intervals.append(motif_interval(pbp, position))
        nearest_seq, seq_distance = nearest_by_sequence(pbp, position)
        nearest_sequence_motifs.append(nearest_seq)
        sequence_distances.append(seq_distance)
        nearest_3d, distance_3d = nearest_by_structure(pbp, position, coords_by_pbp)
        nearest_3d_motifs.append(nearest_3d)
        distances_3d.append(distance_3d)
        prior_category, prior_source = prior.get(
            (pbp, position),
            ("not previously highlighted", "none in existing overlap file"),
        )
        prior_categories.append(prior_category)
        prior_sources.append(prior_source)

    annotated = df.copy()
    annotated["Position"] = parsed_positions
    annotated["Motif interval"] = motif_intervals
    annotated["Nearest motif by sequence"] = nearest_sequence_motifs
    annotated["Sequence distance to motif (aa)"] = sequence_distances
    annotated["Nearest motif by 3D"] = nearest_3d_motifs
    annotated["3D C-alpha distance to motif (Angstrom)"] = distances_3d
    annotated["Previous evidence category"] = prior_categories
    annotated["Previous evidence source"] = prior_sources
    return annotated


def unique_position_frame(annotated: pd.DataFrame) -> pd.DataFrame:
    position_df = annotated[
        [
            "PBP",
            "Position",
            "Motif interval",
            "Nearest motif by sequence",
            "Sequence distance to motif (aa)",
            "Nearest motif by 3D",
            "3D C-alpha distance to motif (Angstrom)",
            "Evidence",
        ]
    ].copy()
    position_df["_rank"] = position_df["Evidence"].map(EVIDENCE_RANK).fillna(-1)
    position_df = position_df.sort_values(["PBP", "Position", "_rank"], ascending=[True, True, False])
    position_df = position_df.drop_duplicates(["PBP", "Position"], keep="first")
    position_df["Structure plot evidence category"] = position_df["Evidence"].where(
        position_df["_rank"] >= EVIDENCE_RANK["Moderate"],
        "Weak/No evidence",
    )
    return position_df.drop(columns=["_rank"])


def resolved_tpd_background(
    coords_by_pbp: dict[str, dict[int, tuple[float, float, float]]]
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for pbp, config in STRUCTURES.items():
        for position in config.transpeptidase_range:
            if position not in coords_by_pbp[pbp]:
                continue
            nearest_seq, seq_distance = nearest_by_sequence(pbp, position)
            nearest_3d, distance_3d = nearest_by_structure(pbp, position, coords_by_pbp)
            records.append(
                {
                    "PBP": pbp,
                    "Position": position,
                    "Motif interval": motif_interval(pbp, position),
                    "Nearest motif by sequence": nearest_seq,
                    "Sequence distance to motif (aa)": seq_distance,
                    "Nearest motif by 3D": nearest_3d,
                    "3D C-alpha distance to motif (Angstrom)": distance_3d,
                }
            )
    return pd.DataFrame.from_records(records)


def histogram_counts(
    background: pd.DataFrame,
    bin_width: int,
) -> pd.DataFrame:
    frame = background.copy()
    frame["_distance"] = pd.to_numeric(
        background["3D C-alpha distance to motif (Angstrom)"], errors="coerce"
    )
    frame = frame.dropna(subset=["_distance"])
    max_distance = float(frame["_distance"].max())
    upper = int(math.ceil(max_distance / bin_width) * bin_width)
    bins = list(range(0, upper + bin_width, bin_width))
    rows = []
    for start, end in zip(bins[:-1], bins[1:]):
        if end == bins[-1]:
            in_bin = frame[(frame["_distance"] >= start) & (frame["_distance"] <= end)]
        else:
            in_bin = frame[(frame["_distance"] >= start) & (frame["_distance"] < end)]
        row = {
            "bin_start_A": start,
            "bin_end_A": end,
            "bin_label": f"{start}-{end}",
            "total_evaluated_count": int(len(in_bin)),
        }
        for label, column, _ in PLOT_CATEGORIES:
            row[column] = int((in_bin["Structure plot evidence category"] == label).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def resolved_tpd_histogram_counts(
    resolved_background: pd.DataFrame,
    evaluated_positions: pd.DataFrame,
    bin_width: int,
) -> pd.DataFrame:
    evaluated_category = {
        (row["PBP"], row["Position"]): row["Structure plot evidence category"]
        for _, row in evaluated_positions.iterrows()
    }

    frame = resolved_background.copy()
    frame["_distance"] = pd.to_numeric(
        frame["3D C-alpha distance to motif (Angstrom)"], errors="coerce"
    )
    frame = frame.dropna(subset=["_distance"])

    def category(row: pd.Series) -> str:
        key = (row["PBP"], row["Position"])
        return evaluated_category.get(key, "Invariant or removed")

    frame["_category"] = frame.apply(category, axis=1)
    max_distance = float(frame["_distance"].max())
    upper = int(math.ceil(max_distance / bin_width) * bin_width)
    bins = list(range(0, upper + bin_width, bin_width))
    rows = []
    for pbp in ["1A", "2B", "2X"]:
        pbp_frame = frame[frame["PBP"] == pbp]
        for start, end in zip(bins[:-1], bins[1:]):
            if end == bins[-1]:
                in_bin = pbp_frame[(pbp_frame["_distance"] >= start) & (pbp_frame["_distance"] <= end)]
            else:
                in_bin = pbp_frame[(pbp_frame["_distance"] >= start) & (pbp_frame["_distance"] < end)]
            invariant_or_removed_count = int((in_bin["_category"] == "Invariant or removed").sum())
            row = {
                "PBP": pbp,
                "bin_start_A": start,
                "bin_end_A": end,
                "bin_label": f"{start}-{end}",
                "invariant_or_removed_count": invariant_or_removed_count,
                "total_resolved_tpd_count": int(len(in_bin)),
            }
            for label, column, _ in PLOT_CATEGORIES:
                row[column] = int((in_bin["_category"] == label).sum())
            rows.append(
                row
            )
    return pd.DataFrame(rows)


def plot_histogram(counts: pd.DataFrame, title: str, background_label: str, path: Path) -> None:
    starts = counts["bin_start_A"].to_numpy()
    widths = (counts["bin_end_A"] - counts["bin_start_A"]).to_numpy()
    bottom = pd.Series(0, index=counts.index).to_numpy()
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for label, column, color in PLOT_CATEGORIES:
        heights = counts[column].to_numpy()
        ax.bar(
            starts,
            heights,
            width=widths,
            align="edge",
            bottom=bottom,
            color=color,
            edgecolor="black",
            linewidth=0.35,
            label=label,
        )
        bottom = bottom + heights
    max_end = counts["bin_end_A"].max()
    ax.set_xlim(0, max_end)
    ax.set_xticks(list(range(0, int(max_end) + 5, 5)))
    ax.set_xlabel("C-alpha distance to nearest active-site motif (Angstrom)")
    ax.set_ylabel("Position count")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(title)
    ax.grid(True, axis="both", color="#D9D9D9", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_resolved_tpd_three_category_histogram(counts: pd.DataFrame, title: str, path: Path) -> None:
    pbp_order = ["1A", "2B", "2X"]
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.5), sharex=True, sharey=False)
    fig.patch.set_facecolor("white")
    handles = labels = None

    max_end = counts["bin_end_A"].max()
    x_ticks = list(range(0, int(max_end) + 5, 5))

    for ax, pbp in zip(axes, pbp_order):
        subset = counts[counts["PBP"] == pbp].copy()
        starts = subset["bin_start_A"].to_numpy()
        widths = (subset["bin_end_A"] - subset["bin_start_A"]).to_numpy()
        bottom = pd.Series(0, index=subset.index).to_numpy()

        ax.set_facecolor("white")
        for label, column, color in PLOT_CATEGORIES:
            heights = subset[column].to_numpy()
            ax.bar(
                starts,
                heights,
                width=widths,
                align="edge",
                bottom=bottom,
                color=color,
                edgecolor="black",
                linewidth=0.35,
                label=label,
            )
            bottom = bottom + heights
        ax.bar(
            starts,
            subset["invariant_or_removed_count"].to_numpy(),
            width=widths,
            align="edge",
            bottom=bottom,
            color=INVARIANT_COLOR,
            edgecolor="black",
            linewidth=0.35,
            label="Invariant or removed",
        )
        ax.set_xlim(0, max_end)
        ax.set_xticks(x_ticks)
        ax.set_title(f"PBP{pbp}", loc="left")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, axis="both", color="#D9D9D9", linewidth=0.7)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.0)
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()

    for ax in axes:
        ax.set_ylabel("Residue count")
    fig.supxlabel("C-alpha distance to nearest active-site motif (Angstrom)")
    fig.suptitle(title, y=0.98)
    if handles is not None:
        fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.945))
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(path, dpi=300)
    plt.close(fig)


def write_outputs(
    annotated: pd.DataFrame,
    outdir: Path,
    coords_by_pbp: dict[str, dict[int, tuple[float, float, float]]],
    bin_width: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    annotated_path = outdir / "Supplementary_File_1_with_structure_context.csv"
    annotated.to_csv(annotated_path, index=False)
    annotated.head(10).to_csv(outdir / "structure_context_top10_preview.csv", index=False)

    positions = unique_position_frame(annotated)
    positions.to_csv(outdir / "structure_context_unique_sf1_positions.csv", index=False)

    evaluated_counts = histogram_counts(positions, bin_width)
    evaluated_counts.to_csv(outdir / "structure_distance_histogram_all_evaluated_positions.csv", index=False)
    plot_histogram(
        evaluated_counts,
        "Distance to nearest active-site motif",
        "All evaluated positions",
        outdir / "structure_distance_histogram_all_evaluated_positions.png",
    )

    resolved_background = resolved_tpd_background(coords_by_pbp)
    resolved_background.to_csv(outdir / "structure_context_resolved_tpd_background.csv", index=False)
    resolved_counts = resolved_tpd_histogram_counts(resolved_background, positions, bin_width)
    resolved_counts.to_csv(outdir / "structure_distance_histogram_all_resolved_tpd_positions.csv", index=False)
    plot_resolved_tpd_three_category_histogram(
        resolved_counts,
        "Distance to nearest active-site motif",
        outdir / "structure_distance_histogram_all_resolved_tpd_positions.png",
    )

    metadata = {
        "purpose": (
            "Post hoc structural and prior-literature annotation. These columns are "
            "not used in evidence scoring."
        ),
        "input_supplementary_file_1": str(DEFAULT_SF1),
        "moderate_plus_definition": (
            "Evidence in Moderate, Strong, or Very Strong, represented as separate "
            "stacked histogram categories after collapsing to unique PBP-position by "
            "the strongest evidence category at that position."
        ),
        "weak_no_evidence_definition": (
            "Evaluated positions whose strongest Supplementary File 1 evidence category is "
            "Weak or Weak/No Evidence."
        ),
        "invariant_or_removed_definition": (
            "Resolved reference-structure residues absent from the Supplementary File 1 "
            "candidate-position frame; this includes invariant, filtered, or otherwise "
            "not evaluated/reported residues."
        ),
        "distance_method": "Minimum C-alpha Euclidean distance from residue to any residue in the conserved active-site motifs.",
        "sequence_numbering": "Direct manuscript/PBP residue numbering was used for the selected PDB chains.",
        "reference_structures": {
            pbp: {
                "pdb_id": config.pdb_id,
                "chain_id": config.chain_id,
                "source_url": config.source_url,
                "motifs": {
                    motif: [positions.start, positions.stop - 1]
                    for motif, positions in config.motifs.items()
                },
                "transpeptidase_range": [
                    config.transpeptidase_range.start,
                    config.transpeptidase_range.stop - 1,
                ],
                "note": config.note,
            }
            for pbp, config in STRUCTURES.items()
        },
        "counts": {
            "supplementary_file_1_rows": int(len(annotated)),
            "unique_evaluated_positions": int(len(positions)),
            "unique_very_strong_positions": int((positions["Structure plot evidence category"] == "Very Strong").sum()),
            "unique_strong_positions": int((positions["Structure plot evidence category"] == "Strong").sum()),
            "unique_moderate_positions": int((positions["Structure plot evidence category"] == "Moderate").sum()),
            "unique_weak_no_evidence_positions": int(
                (positions["Structure plot evidence category"] == "Weak/No evidence").sum()
            ),
            "resolved_tpd_background_positions": int(len(resolved_background)),
            "unresolved_evaluated_positions": int(
                pd.to_numeric(
                    positions["3D C-alpha distance to motif (Angstrom)"], errors="coerce"
                )
                .isna()
                .sum()
            ),
        },
        "outputs": {
            "annotated_supplementary_file_1": str(annotated_path),
            "top10_preview": str(outdir / "structure_context_top10_preview.csv"),
            "evaluated_histogram_png": str(outdir / "structure_distance_histogram_all_evaluated_positions.png"),
            "resolved_tpd_histogram_png": str(outdir / "structure_distance_histogram_all_resolved_tpd_positions.png"),
        },
        "histogram_definitions": {
            "all_evaluated_positions": (
                "Unique PBP positions represented in Supplementary File 1. This is the fairer "
                "comparison for association enrichment because it only includes positions that "
                "entered the candidate/evidence table."
            ),
            "all_resolved_tpd_positions": (
                "All residues in the PBP transpeptidase-domain ranges that are resolved in the "
                "selected reference structures. This is a broader structural background."
            ),
        },
    }
    (outdir / "structure_context_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (outdir / "structure_context_legend.md").write_text(
        """# Structure Context Annotation Legend

These annotations are post hoc biological context only. They were not used to
assign evidence categories, p-values, or effect sizes.

## Reference Structures

- PBP1A: RCSB 2C6W, chain B.
- PBP2B: RCSB 2WAF, chain A.
- PBP2X: RCSB 5OIZ, chain A.

The selected structures use direct PBP residue numbering for the residues
reported in Supplementary File 1.

## Column Definitions

- `Motif interval`: the candidate residue's linear position relative to the
  conserved active-site motif blocks in that PBP.
- `Nearest motif by sequence`: the conserved active-site motif with the smallest amino-acid
  coordinate distance to the candidate residue.
- `Sequence distance to motif (aa)`: the number of amino acids between the
  candidate residue and the nearest conserved active-site motif residue in the primary
  sequence. A value of 0 means the residue is inside a conserved motif block.
- `Nearest motif by 3D`: the conserved active-site motif containing the closest motif
  residue in the folded reference structure.
- `3D C-alpha distance to motif (Angstrom)`: the minimum Euclidean distance
  between the candidate residue C-alpha atom and any C-alpha atom in a conserved active-site
  motif. `not resolved` means the candidate residue or relevant structure
  coordinate is absent from the selected PDB chain.
- `Previous evidence category`: position-level overlap with the existing
  `CDC_GWAS_overlap_TPD.csv` source.
- `Previous evidence source`: expanded source label for the existing overlap
  file. `Chewapreecha et al. GWAS` corresponds to the GWAS column, `Li et al.
  CDC study` corresponds to the CDC column, and `Laboratory (see Supplementary
  S3.2)` corresponds to the Laboratory column.

## Interpretation Notes

Each PBP transpeptidase domain has a single active-site region organised around
three conserved motif blocks: SXXK, SXN, and KTG/KSG. The annotations identify
which motif block is closest to a candidate substitution. The nearest motif can
differ between sequence and 3D columns because residues distant in the linear
sequence can be close together after protein folding.

## Histogram Definitions

- `structure_distance_histogram_all_evaluated_positions.png`: compares Moderate+
  positions against all unique PBP positions represented in Supplementary File 1.
  This is the fairer analysis-background comparison.
- `structure_distance_histogram_all_resolved_tpd_positions.png`: compares
  Moderate+ positions against all residues in the PBP transpeptidase-domain
  ranges that are resolved in the selected reference structures, faceted into
  separate PBP1A, PBP2B, and PBP2X panels. This is a broader structural-background
  comparison. Dark-to-light red shows Very Strong, Strong, and Moderate
  positions; light pink shows Weak/No evidence positions; grey shows residues
  that were not present in the Supplementary File 1 candidate frame because
  they were invariant, filtered, or otherwise not evaluated/reported.

Histograms use counts, continuous 1-Angstrom bins, and stacked bars.
""",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    if args.results_dir is not None:
        if args.sf1 == DEFAULT_SF1:
            args.sf1 = args.results_dir / "manuscript_outputs" / "Supplementary_File_1.csv"
        if args.outdir == DEFAULT_OUTDIR:
            args.outdir = args.results_dir / "manuscript_outputs" / "structure_context"

    pdb_dir = args.pdb_dir or (args.outdir / "pdb")
    coords_by_pbp: dict[str, dict[int, tuple[float, float, float]]] = {}
    provenance: dict[str, dict[str, str]] = {}
    for pbp, config in STRUCTURES.items():
        pdb_path = fetch_pdb(config.pdb_id, pdb_dir)
        coords = parse_ca_coordinates(pdb_path, config.chain_id)
        if not coords:
            raise RuntimeError(f"No C-alpha coordinates parsed for {pbp} {config.pdb_id} chain {config.chain_id}")
        coords_by_pbp[pbp] = coords
        provenance[pbp] = {
            "pdb_id": config.pdb_id,
            "chain_id": config.chain_id,
            "source_url": config.source_url,
            "local_path": str(pdb_path),
            "sha256": _sha256(pdb_path),
        }
    pdb_dir.mkdir(parents=True, exist_ok=True)
    (pdb_dir / "pdb_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    annotated = annotate_sf1(args.sf1, args.prior_overlap, coords_by_pbp)
    write_outputs(annotated, args.outdir, coords_by_pbp, args.bin_width)
    print(f"Wrote structure-context outputs to {args.outdir}")
    print(f"PDB provenance (accession + sha256) written to {pdb_dir / 'pdb_provenance.json'}")


if __name__ == "__main__":
    main()
