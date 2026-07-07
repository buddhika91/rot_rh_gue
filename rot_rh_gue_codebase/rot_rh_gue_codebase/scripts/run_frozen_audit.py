#!/usr/bin/env python3
"""CLI runner for the frozen reflected Stieltjes finite-GUE audit."""
from __future__ import annotations

import argparse

from rot_rh_gue.audit import run_frozen_audit
from rot_rh_gue.utils import parse_int_list, parse_str_list


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-list", default="20000,24000,28000,32000,36000,40000,44000,48000,56000,64000,80000,96000")
    ap.add_argument("--jacobi-dims", default="64,96,128,160,192,224")
    ap.add_argument("--controls", type=int, default=32)
    ap.add_argument("--seed", type=int, default=2027)
    ap.add_argument("--control-modes", default="permuted,signflip,gaussian,phase_scramble")
    ap.add_argument("--oversample", type=int, default=4)
    ap.add_argument("--detrend-degree", type=int, default=3)
    ap.add_argument("--out-prefix", default="reflect_frozen_endpoint_96000_32_seed2027")
    args = ap.parse_args()

    run_frozen_audit(
        N_list=parse_int_list(args.N_list),
        jacobi_dims=parse_int_list(args.jacobi_dims),
        controls=args.controls,
        seed=args.seed,
        control_modes=parse_str_list(args.control_modes),
        oversample=args.oversample,
        detrend_degree=args.detrend_degree,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    main()
