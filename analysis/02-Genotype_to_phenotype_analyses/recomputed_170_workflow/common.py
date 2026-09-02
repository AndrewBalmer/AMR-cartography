#!/usr/bin/env python3
"""Shared utilities for corrected-panel mvLMM recomputation.

This file contains the small pieces of logic reused across the rerun scripts:
historical constants, expected shape checks, marker parsing, Galwey-adjusted
p-values, interaction generation, and model-fitting helpers. Start here when
checking whether raw p-values, adjusted p-values, marker names, and interaction
rules are handled consistently.
"""

from __future__ import annotations

import hashlib
import math
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


def install_pandas_limix_shim() -> None:
    """Provide the pandas.core.index symbol expected by older Limix stacks."""
    if "pandas.core.index" not in sys.modules:
        module = types.ModuleType("pandas.core.index")
        module.InvalidIndexError = pd.errors.InvalidIndexError
        sys.modules["pandas.core.index"] = module
    else:
        setattr(sys.modules["pandas.core.index"], "InvalidIndexError", pd.errors.InvalidIndexError)


install_pandas_limix_shim()

HISTORICAL_ADDITIVE_THRESHOLD = 0.000588
HISTORICAL_UV_THRESHOLD = 0.001
HISTORICAL_EPISTASIS_THRESHOLD = 0.0007620121
HISTORICAL_ADDITIVE_GALWEY_MEFF = 28.0
HISTORICAL_EPISTASIS_GALWEY_MEFF = 39.0
HISTORICAL_EPISTASIS_EFFECT_LOWER_BOUND = 1.0

EXPECTED_OLD_ADDITIVE_MARKERS = 157
EXPECTED_CORRECTED_ADDITIVE_MARKERS = 170
EXPECTED_OLD_EPISTASIS_CANDIDATES = 3542
EXPECTED_CORRECTED_EPISTASIS_CANDIDATES = 4052

# Frozen copy of the historical (157-marker) public Supplementary File 1, as
# published in the bioRxiv release commit c7890d7a76ebe611f0ad6e0001d0dcf8a03bb572.
# The historical comparisons must never read this table from a mutable branch
# pointer such as `main`, which now carries the corrected 170-marker table.
HISTORICAL_PUBLIC_SUPPLEMENT = Path(__file__).resolve().parent / "fixtures" / "Supplementary_File_1_historical_354rows.csv"
HISTORICAL_PUBLIC_SUPPLEMENT_SHA256 = "0fdd239d438dcbab610e3bf05f2611def83f2bd6672704753914bca0a4ec06f5"
HISTORICAL_PUBLIC_SUPPLEMENT_SOURCE_COMMIT = "c7890d7a76ebe611f0ad6e0001d0dcf8a03bb572"
EXPECTED_HISTORICAL_PUBLIC_SUPPLEMENT_ROWS = 354


def load_historical_public_supplement(path: Path | None = None) -> pd.DataFrame:
    """Load the frozen historical public Supplementary File 1, checksum-verified.

    This replaces the previous `git show main:manuscript/Supplementary_File_1.csv`
    lookup, which silently changed meaning once the corrected table was merged
    to main.
    """
    fixture = Path(path) if path is not None else HISTORICAL_PUBLIC_SUPPLEMENT
    if not fixture.exists():
        raise FileNotFoundError(
            f"Historical public supplement fixture not found: {fixture}\n"
            f"Restore it with: git show {HISTORICAL_PUBLIC_SUPPLEMENT_SOURCE_COMMIT}"
            ":manuscript/Supplementary_File_1.csv > " + str(fixture)
        )
    digest = hashlib.sha256(fixture.read_bytes()).hexdigest()
    if digest != HISTORICAL_PUBLIC_SUPPLEMENT_SHA256:
        raise AssertionError(
            f"Historical public supplement fixture checksum mismatch for {fixture}: "
            f"expected {HISTORICAL_PUBLIC_SUPPLEMENT_SHA256}, observed {digest}"
        )
    return pd.read_csv(fixture)


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def galwey_meff(marker_values: pd.DataFrame | np.ndarray, *, floor: bool = True) -> float:
    """Return the Galwey effective-test count used by the historical R workflow.

    The old manuscript scripts hard-coded the floored poolr Galwey values:
    28 for additive/uvLMM marker tests and 39 for epistasis tests. This helper
    is used by validation and sensitivity checks, while exact manuscript
    rebuilds default to the historical literal constants above.
    """
    values = np.asarray(marker_values, dtype=float)
    corr = np.corrcoef(values, rowvar=False)
    eigvals = np.linalg.eigvalsh(corr)
    eigvals = np.maximum(eigvals, 0)
    meff = float(np.square(np.sqrt(eigvals).sum()) / eigvals.sum())
    return float(np.floor(meff)) if floor else meff


def add_galwey_adjusted_pvalues(
    df: pd.DataFrame,
    *,
    raw_column: str = "pv20",
    adjusted_column: str = "pv20_adj_galwey",
    meff: float,
    clip: bool = False,
) -> pd.DataFrame:
    """Add an adjusted p-value column while preserving the raw p-value column."""
    if adjusted_column == raw_column:
        raise ValueError("raw and adjusted p-value columns must be distinct")
    if raw_column not in df.columns:
        raise KeyError(f"Missing raw p-value column: {raw_column}")
    out = df.copy()
    adjusted = pd.to_numeric(out[raw_column], errors="coerce") * float(meff)
    if clip:
        adjusted = adjusted.clip(upper=1.0)
    out[adjusted_column] = adjusted
    out["galwey_meff"] = float(meff)
    return out


def binary_marker_matrix(markers: pd.DataFrame) -> pd.DataFrame:
    """Validate and return a 0/1 marker matrix."""
    bad: list[str] = []
    for column in markers.columns:
        values = set(pd.Series(markers[column]).dropna().unique().tolist())
        if not values.issubset({0, 1, 0.0, 1.0, False, True}):
            bad.append(str(column))
    if bad:
        raise ValueError(f"Non-binary marker columns: {bad[:10]}")
    if markers.isna().any().any():
        missing = markers.columns[markers.isna().any()].astype(str).tolist()
        raise ValueError(f"Marker columns contain missing values: {missing[:10]}")
    return markers.astype(np.uint8)


def four_cell_interaction_metadata(
    markers: pd.DataFrame,
    *,
    min_fraction: float = 0.01,
) -> pd.DataFrame:
    """Generate pairwise marker products using the original four-cell rule.

    The historical R code used `model.matrix(Y ~ .^2)` and retained interaction
    products only when all four parent genotype cells were at least 1% of
    isolates. It did not hash-deduplicate interaction products for the exact
    replication output.
    """
    matrix = binary_marker_matrix(markers)
    names = list(matrix.columns.astype(str))
    min_count = float(len(matrix) * min_fraction)
    rows: list[dict[str, object]] = []

    for i, parent_a in enumerate(names):
        a = matrix[parent_a].to_numpy(dtype=np.uint8)
        for parent_b in names[i + 1 :]:
            b = matrix[parent_b].to_numpy(dtype=np.uint8)
            n_00 = int(((a == 0) & (b == 0)).sum())
            n_10 = int(((a == 1) & (b == 0)).sum())
            n_01 = int(((a == 0) & (b == 1)).sum())
            n_11 = int(((a == 1) & (b == 1)).sum())
            if min(n_00, n_10, n_01, n_11) < min_count:
                continue
            interaction = f"{parent_a}:{parent_b}"
            rows.append(
                {
                    "interaction": interaction,
                    "parent_a": parent_a,
                    "parent_b": parent_b,
                    "n_00": n_00,
                    "n_10": n_10,
                    "n_01": n_01,
                    "n_11": n_11,
                    "n_present": n_11,
                    "frequency": float(n_11 / len(matrix)),
                    "min_cell_count": min(n_00, n_10, n_01, n_11),
                    "one_percent_count": min_count,
                }
            )

    return pd.DataFrame(rows)


def normalized_relatedness_features(values: np.ndarray) -> np.ndarray:
    """Return the feature matrix whose cross-product equals limix.linear_kinship."""
    x = np.asarray(values, dtype=float)
    means = np.nanmean(x, axis=0)
    x = np.where(np.isnan(x), means, x)
    x = x - means
    std = np.std(x, axis=0)
    if np.any(std == 0):
        bad = int(np.sum(std == 0))
        raise ValueError(f"Cannot build kinship with {bad} zero-variance columns")
    x = x / std
    x = x / math.sqrt(x.shape[1])
    return x


def marker_order_from_effects(effect_file: Path) -> list[str]:
    effects = read_csv(effect_file)
    candidate = effects[(effects["effect_type"] == "candidate") & (effects["env"] == "env1_D1")]
    return candidate["effect_name"].astype(str).tolist()


def chunk_bounds(total: int, chunk_size: int, array_index: int) -> tuple[int, int]:
    if array_index < 1:
        raise ValueError("Chunk array indexes are 1-based; got array_index < 1")
    start = (array_index - 1) * chunk_size
    end = min(start + chunk_size, total)
    if start >= total:
        raise ValueError(f"Chunk start {start} exceeds total rows {total}")
    return start, end


def lowrank_multitrait_scan(
    *,
    y: np.ndarray,
    trait_design: np.ndarray,
    covariates: np.ndarray,
    relatedness_features: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, object]:
    """Run the historical epistatic mvLMM model using low-rank relatedness features.

    This is algebraically the same model form as Limix's multivariate scan with a
    linear kinship matrix, but avoids constructing/eigendecomposing a dense
    n-by-n kinship matrix for every interaction.
    """
    from glimix_core.lmm import Kron2Sum
    from numpy_sugar.linalg import ddot, economic_qs_linear

    qs = economic_qs_linear(normalized_relatedness_features(relatedness_features))
    kg = ddot(qs[0][0], np.sqrt(qs[1]))

    lmm = Kron2Sum(y, trait_design, covariates, kg, restricted=False)
    lmm.fit(verbose=False)
    scanner = lmm.get_fast_scanner()
    result = scanner.scan(np.eye(y.shape[1]), candidate)

    lml0 = float(lmm.lml())
    lml2 = float(result["lml"])
    lrt = max(0.0, 2.0 * (lml2 - lml0))
    return {
        "lml0": lml0,
        "lml2": lml2,
        "dof20": y.shape[1],
        "scale2": float(result["scale"]),
        "pv20": float(chi2.sf(lrt, y.shape[1])),
        "candidate_eff": np.asarray(result["effsizes1"], dtype=float).reshape(-1),
        "candidate_se": np.asarray(result["effsizes1_se"], dtype=float).reshape(-1),
        "covariate_eff": np.asarray(result["effsizes0"], dtype=float),
        "covariate_se": np.asarray(result["effsizes0_se"], dtype=float),
    }
