"""Frozen finite-GUE audit runner."""
from __future__ import annotations

from dataclasses import asdict
from typing import Any

import numpy as np

from .arithmetic import build_cutoffs
from .channels import ChannelConfig, direct_binned_channels
from .controls import DEFAULT_CONTROL_MODES, make_control_channels
from .gue import local_unfolded_spacing_ks
from .operator import FrozenOperator, build_jacobi_spectrum, frozen_operators
from .utils import mean_std, write_csv


def run_frozen_audit(
    *,
    N_list: list[int],
    jacobi_dims: list[int],
    controls: int,
    seed: int,
    out_prefix: str,
    control_modes: list[str] | None = None,
    oversample: int = 4,
    detrend_degree: int = 3,
    operators: list[FrozenOperator] | None = None,
    verbose: bool = True,
) -> dict[str, list[dict[str, Any]]]:
    """Run the no-search frozen audit and write CSV outputs."""

    rng = np.random.default_rng(seed)
    modes = control_modes or DEFAULT_CONTROL_MODES
    ops = operators or frozen_operators()
    cfg = ChannelConfig(oversample=oversample, detrend_degree=detrend_degree)

    if verbose:
        print("=" * 132)
        print("ROT RH / GLOBAL GUE — FROZEN REFLECTED STIELTJES AUDIT")
        print("=" * 132)
        print(f"N_list       : {N_list}")
        print(f"jacobi_dims  : {jacobi_dims}")
        print(f"N denominator: {len(N_list)}")
        print(f"controls     : {controls}")
        print(f"control_modes: {modes}")
        print(f"operators    : {[op.key for op in ops]}")
        print(f"seed         : {seed}")
        print(f"out_prefix   : {out_prefix}")
        print("=" * 132)
        print("\nPrecomputing Mangoldt arrays...")

    cutoffs = build_cutoffs(N_list, verbose=verbose)

    real_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []

    if verbose:
        print("\nComputing frozen operators...")

    for op_id, op in enumerate(ops):
        if verbose:
            print(f"  operator {op_id}: {op.key}")
        for d in jacobi_dims:
            M = max(d * oversample, d + 8)
            for N in N_list:
                ch = direct_binned_channels(cutoffs[N], M, cfg)
                vals = build_jacobi_spectrum(ch, d, op)
                real_ks = local_unfolded_spacing_ks(vals, op.window_lo, op.window_hi, op.local_width)

                base = {
                    "operator_id": op_id,
                    "operator_key": op.key,
                    "source": "real",
                    "jacobi_dim": d,
                    "support_M": M,
                    "N": N,
                    "lo": op.window_lo,
                    "hi": op.window_hi,
                    "local_width": op.local_width,
                    **asdict(op),
                }
                base["ks"] = real_ks
                real_rows.append(base)

                for ctrl_id in range(controls):
                    for mode in modes:
                        cch = make_control_channels(ch, mode, rng)
                        cvals = build_jacobi_spectrum(cch, d, op)
                        cks = local_unfolded_spacing_ks(cvals, op.window_lo, op.window_hi, op.local_width)
                        control_rows.append({
                            "operator_id": op_id,
                            "operator_key": op.key,
                            "source": "control",
                            "control_id": ctrl_id,
                            "control_mode": mode,
                            "jacobi_dim": d,
                            "support_M": M,
                            "N": N,
                            "lo": op.window_lo,
                            "hi": op.window_hi,
                            "local_width": op.local_width,
                            "ks": cks,
                        })

    summary: list[dict[str, Any]] = []
    for rr in real_rows:
        op_id = int(rr["operator_id"])
        d = int(rr["jacobi_dim"])
        N = int(rr["N"])
        vals = [
            float(c["ks"])
            for c in control_rows
            if int(c["operator_id"]) == op_id and int(c["jacobi_dim"]) == d and int(c["N"]) == N
        ]
        mu, sd = mean_std(vals)
        ksz = (mu - float(rr["ks"])) / sd if sd > 1e-12 else float("nan")
        out = dict(rr)
        out["control_ks_mean"] = mu
        out["control_ks_std"] = sd
        out["KSz"] = ksz
        summary.append(out)

    by_dim: list[dict[str, Any]] = []
    aggregate: list[dict[str, Any]] = []
    denom = len(N_list)

    for op_id, op in enumerate(ops):
        dim_rows = []
        for d in jacobi_dims:
            rows = [r for r in summary if int(r["operator_id"]) == op_id and int(r["jacobi_dim"]) == d]
            vals = np.array([float(r["KSz"]) for r in rows], dtype=np.float64)
            byN = {int(r["N"]): float(r["KSz"]) for r in rows}
            br = {
                "operator_id": op_id,
                "operator_key": op.key,
                "jacobi_dim": d,
                "N_denominator": denom,
                "pos_count": int(np.sum(vals > 0)),
                "strong_count": int(np.sum(vals > 1)),
                "mean_KSz": float(np.mean(vals)),
                "min_KSz": float(np.min(vals)),
            }
            for N in N_list:
                br[f"KSz_{N}"] = byN.get(N, float("nan"))
            by_dim.append(br)
            dim_rows.append(br)

        all_vals = np.array([float(r["KSz"]) for r in summary if int(r["operator_id"]) == op_id], dtype=np.float64)
        ar = {
            "operator_id": op_id,
            "operator_key": op.key,
            "build_mode": op.build_mode,
            "lo": op.window_lo,
            "hi": op.window_hi,
            "local_width": op.local_width,
            "N_denominator": denom,
            "min_pos_count_across_dims": int(min(int(x["pos_count"]) for x in dim_rows)),
            "mean_pos_count_across_dims": float(np.mean([int(x["pos_count"]) for x in dim_rows])),
            "min_min_KSz_across_dims": float(np.min(all_vals)),
            "mean_KSz_across_dims": float(np.mean(all_vals)),
        }
        for N in N_list:
            ar[f"min_KSz{N}_across_dims"] = float(min(float(x.get(f"KSz_{N}", float("nan"))) for x in dim_rows))
        aggregate.append(ar)

    aggregate = sorted(
        aggregate,
        key=lambda r: (
            int(r["min_pos_count_across_dims"]),
            float(r["mean_pos_count_across_dims"]),
            float(r["min_min_KSz_across_dims"]),
            float(r["mean_KSz_across_dims"]),
        ),
        reverse=True,
    )

    if verbose:
        print("\nFROZEN AUDIT LEADERBOARD")
        print("-" * 132)
        for r in aggregate:
            print(
                f"{r['operator_key']:<42s} "
                f"mode={r['build_mode']:<10s} "
                f"win=[{float(r['lo']):.2f},{float(r['hi']):.2f}] "
                f"L={int(r['local_width']):<2d} "
                f"minPos={int(r['min_pos_count_across_dims'])}/{int(r['N_denominator'])} "
                f"meanPos={float(r['mean_pos_count_across_dims']):.2f}/{int(r['N_denominator'])} "
                f"minMinKSz={float(r['min_min_KSz_across_dims']):+.3f} "
                f"meanKSz={float(r['mean_KSz_across_dims']):+.3f}"
            )

    write_csv(f"{out_prefix}_real_rows.csv", real_rows)
    write_csv(f"{out_prefix}_control_rows.csv", control_rows)
    write_csv(f"{out_prefix}_summary.csv", summary)
    write_csv(f"{out_prefix}_leaderboard_by_dim.csv", by_dim)
    write_csv(f"{out_prefix}_aggregate.csv", aggregate)

    return {
        "real_rows": real_rows,
        "control_rows": control_rows,
        "summary": summary,
        "leaderboard_by_dim": by_dim,
        "aggregate": aggregate,
    }
