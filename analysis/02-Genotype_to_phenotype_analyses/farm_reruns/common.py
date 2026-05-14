#!/usr/bin/env python3
"""Shared utilities for corrected-panel mvLMM farm reruns."""

from __future__ import annotations

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

from glimix_core.lmm import Kron2Sum  # noqa: E402
from numpy import eye, sqrt  # noqa: E402
from numpy_sugar.linalg import ddot, economic_qs_linear  # noqa: E402


DEFAULT_NEW_THRESHOLD = 0.0009078488974311251
DEFAULT_UV_THRESHOLD = 0.001


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


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
        raise ValueError("LSF array indexes are 1-based; got array_index < 1")
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
    qs = economic_qs_linear(normalized_relatedness_features(relatedness_features))
    kg = ddot(qs[0][0], sqrt(qs[1]))

    lmm = Kron2Sum(y, trait_design, covariates, kg, restricted=False)
    lmm.fit(verbose=False)
    scanner = lmm.get_fast_scanner()
    result = scanner.scan(eye(y.shape[1]), candidate)

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
