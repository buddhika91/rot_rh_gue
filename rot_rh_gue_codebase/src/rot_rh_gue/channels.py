"""Arithmetic channels used to build reflected Stieltjes measures."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np

from .arithmetic import ArithmeticCutoff
from .utils import detrend, moving_average, zscore


@dataclass(frozen=True)
class ChannelConfig:
    oversample: int = 4
    detrend_degree: int = 3


def direct_binned_channels(cutoff: ArithmeticCutoff, support_size: int, cfg: ChannelConfig) -> Dict[str, np.ndarray]:
    """Build the direct-binned Mangoldt channel family.

    The main channel is the binned defect

        F_j = sum_{n in bin j} Λ(n) - bin_width_j.

    We also expose derivative/smoothed/reflected versions.  The reflected
    channel F(-x) was the crucial stabilizer in the successful audits.
    """

    N = cutoff.N
    lam = cutoff.mangoldt
    lam1 = lam[1:]
    M = int(support_size)

    edges = np.linspace(1, N + 1, M + 1).astype(int)
    mass = np.zeros(M, dtype=np.float64)
    expected = np.zeros(M, dtype=np.float64)
    support = np.zeros(M, dtype=np.float64)
    nz = lam > 0

    for j in range(M):
        a = max(1, min(int(edges[j]), N + 1))
        b = max(1, min(int(edges[j + 1]), N + 1))
        if b <= a:
            b = min(N + 1, a + 1)
        mass[j] = float(np.sum(lam[a:b]))
        expected[j] = float(b - a)
        support[j] = float(np.sum(nz[a:b]))

    mass_defect = zscore(detrend(mass - expected, cfg.detrend_degree))
    support_defect = zscore(detrend(support - np.mean(support), cfg.detrend_degree))

    psi = np.cumsum(lam1)
    n = np.arange(1, N + 1, dtype=np.float64)
    psi_defect = psi - n
    psi_grid = np.interp(
        np.linspace(0.0, 1.0, M),
        np.linspace(0.0, 1.0, len(psi_defect)),
        psi_defect,
    )
    psi_grid = zscore(detrend(psi_grid, cfg.detrend_degree))

    mass_reflect = zscore(mass_defect[::-1].copy())

    return {
        "mass": mass_defect,
        "mass_dx": zscore(np.gradient(mass_defect)),
        "mass_smooth": zscore(moving_average(mass_defect, 9)),
        "mass_reflect": mass_reflect,
        "mass_reflect_dx": zscore(np.gradient(mass_reflect)),
        "mass_reflect_smooth": zscore(moving_average(mass_reflect, 9)),
        "psi": psi_grid,
        "support": support_defect,
        "psi_minus_N": np.array([cutoff.psi_minus_N], dtype=np.float64),
    }
