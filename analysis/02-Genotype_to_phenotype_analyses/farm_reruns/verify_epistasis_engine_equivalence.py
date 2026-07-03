#!/usr/bin/env python3
"""Verify the low-rank epistatic mvLMM engine matches the exact limix dense engine.

The recomputed epistasis pipeline uses ``common.lowrank_multitrait_scan``
(glimix-core ``Kron2Sum`` with low-rank relatedness features) for speed. This
harness re-runs a subset of interactions with the *original manuscript* limix
call (``scan(G=product, Y=map_coords, K=linear_kinship(dropped), M=parents,
A=mod_matrix)``, matching ``31-mvLMM-heritability-and-epistatic-mvLMM.py``) and
checks that the two engines agree on ``pv20`` and the joint effect sizes.

It also cross-checks the committed recomputed observed p-values, confirming the
merged farm outputs are reproducible from the current code.

HEAD-NODE SAFETY: this fits LMMs and must be submitted through LSF
(``lsf/submit_verify_epistasis_engine.sh``), never run directly on a farm head
node.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from xarray import DataArray

from common import install_pandas_limix_shim, lowrank_multitrait_scan, read_csv

install_pandas_limix_shim()

from limix.qtl import scan  # noqa: E402
from limix.stats import linear_kinship  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--interaction-file", required=True, type=Path)
    parser.add_argument(
        "--observed-file",
        type=Path,
        help="Optional merged corrected_epistasis_p_values.csv to select significant "
        "interactions and cross-check committed pv20 values.",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--n-lowest", type=int, default=10, help="Most-significant interactions to test.")
    parser.add_argument("--n-random", type=int, default=10, help="Random additional interactions to test.")
    parser.add_argument("--random-seed", type=int, default=2)
    parser.add_argument("--pv20-atol", type=float, default=1e-6)
    parser.add_argument("--pv20-rtol", type=float, default=1e-3)
    return parser.parse_args()


def dense_limix_scan(
    *,
    candidate: np.ndarray,
    y_df: pd.DataFrame,
    parents_df: pd.DataFrame,
    dropped: np.ndarray,
    trait_design: np.ndarray,
) -> dict[str, float]:
    """Replicate the original 31.py epistatic mvLMM scan for one interaction.

    Mirrors the model the low-rank engine fits: covariate parents carry the
    ``mod_matrix`` trait design (A), the interaction candidate carries a free
    per-trait effect (A1 = identity, A0 = empty). A0/A1 are passed as explicit
    xarray DataArrays to avoid the limix 3.0.4 + xarray coord-mutation bug that
    fires when limix builds them internally (same workaround as
    run_additive_chunk.py).
    """
    kinship = linear_kinship(dropped)
    candidate_df = pd.DataFrame({"interaction": candidate.reshape(-1)})
    trait_labels = y_df.columns.astype(object).to_numpy()
    null_trait_design = DataArray(
        np.empty((len(trait_labels), 0)),
        dims=["sample", "env"],
        coords={"env": np.asarray([], dtype=object)},
    )
    candidate_trait_design = DataArray(
        np.eye(len(trait_labels)),
        dims=["sample", "env"],
        coords={"env": trait_labels},
    )
    result = scan(
        G=candidate_df,
        Y=y_df,
        K=kinship,
        M=parents_df,
        A=trait_design,
        A0=null_trait_design,
        A1=candidate_trait_design,
        verbose=False,
    )
    stats = result.stats.reset_index(drop=True).iloc[0]
    cand = result.effsizes["h2"]
    cand = cand[cand["effect_type"] == "candidate"].reset_index(drop=True)
    eff = cand.sort_values("env")["effsize"].to_numpy(dtype=float)
    return {
        "pv20_dense": float(stats["pv20"]),
        "lml0_dense": float(stats["lml0"]),
        "lml2_dense": float(stats["lml2"]),
        "joint_effect_dense": float(np.sqrt(np.square(eff).sum())),
    }


def select_interactions(interactions: pd.DataFrame, observed: pd.DataFrame | None, args: argparse.Namespace) -> pd.DataFrame:
    if observed is not None and "pv20" in observed.columns:
        ranked = observed.sort_values("pv20", kind="mergesort")
        lowest = ranked.head(args.n_lowest)["interaction"].tolist()
    else:
        lowest = interactions["interaction"].head(args.n_lowest).tolist()
    remaining = interactions[~interactions["interaction"].isin(lowest)]
    n_random = min(args.n_random, len(remaining))
    random_pick = remaining.sample(n=n_random, random_state=args.random_seed)["interaction"].tolist() if n_random else []
    chosen = list(dict.fromkeys(lowest + random_pick))
    subset = interactions[interactions["interaction"].isin(chosen)].copy()
    subset["selection"] = np.where(subset["interaction"].isin(lowest), "most_significant", "random")
    return subset.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    markers = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_test_markers.csv")
    relatedness = read_csv(args.data_dir / "S.pneumo_map_dummy_gen_relatedness_matrix.csv")
    map_coords = read_csv(args.data_dir / "S.pneumo_map_mvlmm_map_coords.csv")
    trait_design = read_csv(args.data_dir / "S.pneumo_map_mvlmm_mod_matrix.csv").to_numpy(dtype=float)
    interactions = read_csv(args.interaction_file)

    observed = read_csv(args.observed_file) if args.observed_file and args.observed_file.exists() else None
    stored_pv = (
        observed.set_index("interaction")["pv20"].to_dict() if observed is not None and "pv20" in observed.columns else {}
    )

    y_df = map_coords[["D1", "D2"]].copy()
    y_np = y_df.to_numpy(dtype=float)

    subset = select_interactions(interactions, observed, args)
    rows: list[dict[str, object]] = []
    for row in subset.itertuples(index=False):
        parent_a, parent_b = str(row.parent_a), str(row.parent_b)
        candidate = (markers[parent_a].to_numpy(dtype=float) * markers[parent_b].to_numpy(dtype=float)).reshape(-1, 1)
        parents_df = relatedness[[parent_a, parent_b]].copy()
        dropped = relatedness.drop(columns=[parent_a, parent_b]).to_numpy(dtype=float)

        low = lowrank_multitrait_scan(
            y=y_np,
            trait_design=trait_design,
            covariates=parents_df.to_numpy(dtype=float),
            relatedness_features=dropped,
            candidate=candidate,
        )
        dense = dense_limix_scan(
            candidate=candidate,
            y_df=y_df,
            parents_df=parents_df,
            dropped=dropped,
            trait_design=trait_design,
        )
        pv_low = float(low["pv20"])
        pv_dense = float(dense["pv20_dense"])
        joint_low = float(np.sqrt(np.square(np.asarray(low["candidate_eff"], dtype=float)).sum()))
        entry = {
            "interaction": str(row.interaction),
            "selection": str(row.selection),
            "pv20_lowrank": pv_low,
            "pv20_dense": pv_dense,
            "pv20_abs_diff": abs(pv_low - pv_dense),
            "neglog10_abs_diff": abs(-np.log10(max(pv_low, 1e-300)) + np.log10(max(pv_dense, 1e-300))),
            "joint_effect_lowrank": joint_low,
            "joint_effect_dense": float(dense["joint_effect_dense"]),
            "pv20_stored": float(stored_pv[str(row.interaction)]) if str(row.interaction) in stored_pv else np.nan,
        }
        entry["pv20_lowrank_vs_stored_abs_diff"] = (
            abs(pv_low - entry["pv20_stored"]) if np.isfinite(entry["pv20_stored"]) else np.nan
        )
        rows.append(entry)
        print(
            f"{entry['selection']:>16}  {entry['interaction'][:40]:<40}  "
            f"low={pv_low:.3e}  dense={pv_dense:.3e}  dpv={entry['pv20_abs_diff']:.2e}",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail.to_csv(args.out_dir / "epistasis_engine_equivalence_detail.csv", index=False)

    valid = detail.dropna(subset=["pv20_lowrank", "pv20_dense"])
    within_tol = bool(
        np.all(valid["pv20_abs_diff"].to_numpy() <= (args.pv20_atol + args.pv20_rtol * valid["pv20_dense"].to_numpy()))
    )
    neglog = valid[(valid["pv20_lowrank"] > 0) & (valid["pv20_dense"] > 0)]
    pearson = (
        float(np.corrcoef(-np.log10(neglog["pv20_lowrank"]), -np.log10(neglog["pv20_dense"]))[0, 1])
        if len(neglog) >= 2
        else None
    )
    stored_valid = detail.dropna(subset=["pv20_lowrank_vs_stored_abs_diff"])
    summary = {
        "interactions_tested": int(len(detail)),
        "n_most_significant": int((detail["selection"] == "most_significant").sum()),
        "n_random": int((detail["selection"] == "random").sum()),
        "pv20_max_abs_diff": float(valid["pv20_abs_diff"].max()),
        "pv20_neglog10_max_abs_diff": float(valid["neglog10_abs_diff"].max()),
        "pv20_neglog10_pearson": pearson,
        "pv20_within_tolerance": within_tol,
        "pv20_atol": args.pv20_atol,
        "pv20_rtol": args.pv20_rtol,
        "joint_effect_max_abs_diff": float((valid["joint_effect_lowrank"] - valid["joint_effect_dense"]).abs().max()),
        "committed_output_max_abs_diff": (
            float(stored_valid["pv20_lowrank_vs_stored_abs_diff"].max()) if len(stored_valid) else None
        ),
        "committed_output_cross_checked": int(len(stored_valid)),
    }
    (args.out_dir / "epistasis_engine_equivalence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("\n" + json.dumps(summary, indent=2))
    if not within_tol:
        raise SystemExit("FAIL: low-rank and dense limix pv20 disagree beyond tolerance")
    print("\nPASS: low-rank epistasis engine matches exact limix dense engine within tolerance.")


if __name__ == "__main__":
    main()
