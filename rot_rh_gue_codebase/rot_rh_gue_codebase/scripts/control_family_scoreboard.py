#!/usr/bin/env python3
"""Control-family scoreboard for a frozen audit output prefix."""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from typing import Any

from rot_rh_gue.controls import ARITHMETIC_DISRUPTION_MODES
from rot_rh_gue.utils import read_csv, write_csv


def f(x: Any, default: float = float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default


def i(x: Any, default: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return default


def mean_sd(vals: list[float]) -> tuple[float, float]:
    if not vals:
        return float("nan"), float("nan")
    mu = sum(vals) / len(vals)
    var = sum((x - mu) ** 2 for x in vals) / len(vals)
    return mu, math.sqrt(var)


def fmt(x: float) -> str:
    return "nan" if math.isnan(x) else f"{x:+.3f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out-prefix", default="")
    args = ap.parse_args()

    prefix = args.prefix
    out_prefix = args.out_prefix or f"{prefix}_control_family"
    summary = read_csv(f"{prefix}_summary.csv")
    controls = read_csv(f"{prefix}_control_rows.csv")

    dims = sorted({i(r.get("jacobi_dim")) for r in summary})
    Ns = sorted({i(r.get("N")) for r in summary})
    denom = len(Ns)
    op_keys = sorted({r.get("operator_key", r.get("lock_key", "")) for r in summary})

    ctrl_lookup: dict[tuple[str, int, int, str], list[float]] = defaultdict(list)
    for c in controls:
        op = c.get("operator_key", c.get("lock_key", ""))
        d = i(c.get("jacobi_dim"))
        N = i(c.get("N"))
        mode = c.get("control_mode", "")
        val = f(c.get("ks"))
        ctrl_lookup[(op, d, N, mode)].append(val)
        ctrl_lookup[(op, d, N, "all_controls")].append(val)
        if mode in ARITHMETIC_DISRUPTION_MODES:
            ctrl_lookup[(op, d, N, "arithmetic_disrupted")].append(val)

    families = ["all_controls", "arithmetic_disrupted", "gaussian", "permuted", "phase_scramble", "signflip"]
    family_rows = []
    cell_rows = []

    print("=" * 132)
    print("ROT RH / CONTROL-FAMILY SCOREBOARD")
    print("=" * 132)
    print(f"prefix : {prefix}")
    print(f"dims   : {dims}")
    print(f"Ns     : {Ns}")
    print("=" * 132)

    for op in op_keys:
        print(f"\n{op}")
        print("-" * 132)
        op_summary = [r for r in summary if r.get("operator_key", r.get("lock_key", "")) == op]
        rows_for_op = []
        for fam in families:
            cells = []
            for r in op_summary:
                d = i(r.get("jacobi_dim"))
                N = i(r.get("N"))
                real_ks = f(r.get("ks"))
                vals = ctrl_lookup[(op, d, N, fam)]
                mu, sd = mean_sd(vals)
                ksz = (mu - real_ks) / sd if sd > 1e-12 else float("nan")
                cell = {"operator_key": op, "family": fam, "jacobi_dim": d, "N": N, "real_ks": real_ks, "control_mean": mu, "control_std": sd, "KSz": ksz}
                cells.append(cell)
                cell_rows.append(cell)
            dim_counts = []
            for d in dims:
                vals = [x["KSz"] for x in cells if x["jacobi_dim"] == d]
                dim_counts.append(sum(x > 0 for x in vals))
            ksz_vals = [x["KSz"] for x in cells]
            worst = min(cells, key=lambda x: x["KSz"])
            fr = {
                "operator_key": op,
                "family": fam,
                "N_denominator": denom,
                "min_pos_count_across_dims": min(dim_counts),
                "mean_pos_count_across_dims": sum(dim_counts) / len(dim_counts),
                "min_KSz": min(ksz_vals),
                "mean_KSz": sum(ksz_vals) / len(ksz_vals),
                "worst_dim": worst["jacobi_dim"],
                "worst_N": worst["N"],
                "worst_KSz": worst["KSz"],
            }
            rows_for_op.append(fr)
            family_rows.append(fr)
        rows_for_op.sort(key=lambda r: (r["min_pos_count_across_dims"], r["mean_pos_count_across_dims"], r["min_KSz"]), reverse=True)
        for r in rows_for_op:
            print(
                f"{r['family']:<22s} minPos={int(r['min_pos_count_across_dims'])}/{denom} "
                f"meanPos={float(r['mean_pos_count_across_dims']):.2f}/{denom} "
                f"minKSz={fmt(float(r['min_KSz']))} meanKSz={fmt(float(r['mean_KSz']))} "
                f"worst=d{int(r['worst_dim'])},N{int(r['worst_N'])},KSz={fmt(float(r['worst_KSz']))}"
            )

    write_csv(f"{out_prefix}_family_summary.csv", family_rows)
    write_csv(f"{out_prefix}_cell_rows.csv", cell_rows)
    print("\nWrote:")
    print(f"  {out_prefix}_family_summary.csv")
    print(f"  {out_prefix}_cell_rows.csv")


if __name__ == "__main__":
    main()
