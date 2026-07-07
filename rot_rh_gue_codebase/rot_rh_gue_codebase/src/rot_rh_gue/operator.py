"""Frozen reflected two-channel Mangoldt-Stieltjes operators."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np

from .stieltjes import (
    centered_diagonal_mix,
    cosine_support,
    jacobi_eigenvalues,
    normalize_offdiagonal,
    stieltjes_coefficients,
)
from .utils import zscore


@dataclass(frozen=True)
class FrozenOperator:
    """A fully specified finite-GUE operator and measurement choice."""

    key: str
    build_mode: str  # coeff_mix or bgeom_mix
    A: float
    beta: float
    dx_mix: float
    smooth_mix: float
    diag_scale: float
    offdiag_scale: float
    edge_taper: float
    gamma: float
    window_lo: float
    window_hi: float
    local_width: int


def frozen_operators() -> list[FrozenOperator]:
    """Return the two frozen finalist operators.

    The first is the strongest coefficient-mix local-width lock found in the
    64k audit.  The second is the cleaner 96k endpoint near-lock and the main
    report candidate because it separates arithmetic-disruption controls most
    cleanly.
    """

    return [
        FrozenOperator(
            key="LOCK_coeff_ba115_sm025_g045_W020_080_L7",
            build_mode="coeff_mix",
            A=2.20,
            beta=1.15,
            dx_mix=0.10,
            smooth_mix=0.25,
            diag_scale=0.04,
            offdiag_scale=0.85,
            edge_taper=0.25,
            gamma=0.45,
            window_lo=0.20,
            window_hi=0.80,
            local_width=7,
        ),
        FrozenOperator(
            key="LOCK_bgeom_best_W012_088_L13",
            build_mode="bgeom_mix",
            A=2.20,
            beta=1.35,
            dx_mix=0.10,
            smooth_mix=0.25,
            diag_scale=0.04,
            offdiag_scale=0.85,
            edge_taper=0.25,
            gamma=0.55,
            window_lo=0.12,
            window_hi=0.88,
            local_width=13,
        ),
    ]


def _measure_from_channel(ch: Dict[str, np.ndarray], op: FrozenOperator, reflect: bool) -> tuple[np.ndarray, np.ndarray]:
    M = len(ch["mass"])
    x = cosine_support(M)

    if reflect:
        y = ch["mass_reflect"]
        dx = ch["mass_reflect_dx"]
        sm = ch["mass_reflect_smooth"]
    else:
        y = ch["mass"]
        dx = ch["mass_dx"]
        sm = ch["mass_smooth"]

    driver = zscore(y + op.dx_mix * dx + op.smooth_mix * sm)

    # Positive Stieltjes measure.  The exponential puts arithmetic into the
    # weight; the edge factor regularizes the finite support.
    w = np.exp(np.clip(op.A * driver, -10.0, 10.0))
    w *= np.maximum(1e-12, 1.0 - x * x) ** op.beta
    w = np.maximum(w, 1e-14)
    w = w / np.sum(w)
    return x, w


def build_jacobi_spectrum(ch: Dict[str, np.ndarray], dim: int, op: FrozenOperator) -> np.ndarray:
    """Build eigenvalues for a frozen reflected two-channel operator."""

    xp, wp = _measure_from_channel(ch, op, reflect=False)
    xm, wm = _measure_from_channel(ch, op, reflect=True)

    ap, bp = stieltjes_coefficients(xp, wp, dim)
    am, bm = stieltjes_coefficients(xm, wm, dim)

    bp2 = normalize_offdiagonal(bp, dim, op.edge_taper, op.offdiag_scale)
    bm2 = normalize_offdiagonal(bm, dim, op.edge_taper, op.offdiag_scale)

    a = centered_diagonal_mix(ap, am, op.gamma, op.diag_scale)

    if op.build_mode == "coeff_mix":
        b = (1.0 - op.gamma) * bp2 + op.gamma * bm2
    elif op.build_mode == "bgeom_mix":
        b = np.exp((1.0 - op.gamma) * np.log(np.maximum(bp2, 1e-14)) + op.gamma * np.log(np.maximum(bm2, 1e-14)))
    else:
        raise ValueError(f"unknown build_mode={op.build_mode!r}")

    return jacobi_eigenvalues(a, b)
