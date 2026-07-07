"""Arithmetic primitives for the ROT-RH finite-GUE benchmark.

This module intentionally avoids dependencies beyond NumPy.  The core
arithmetic input is the von Mangoldt function Λ(n), which marks prime powers
with log(p) and all other integers with zero.  The benchmark then studies
finite cutoffs N through the arithmetic defect ψ(N)-N.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class ArithmeticCutoff:
    """Cached arithmetic data for a cutoff N."""

    N: int
    mangoldt: np.ndarray

    @property
    def psi_minus_N(self) -> float:
        return float(np.sum(self.mangoldt[1:]) - self.N)


def mangoldt_array(N: int) -> np.ndarray:
    """Return Λ(n) for 0 <= n <= N.

    Λ(p^k) = log(p) for prime p and k>=1; otherwise 0.
    """

    if N < 2:
        return np.zeros(N + 1, dtype=np.float64)

    lam = np.zeros(N + 1, dtype=np.float64)
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[:2] = False

    for p in range(2, int(math.isqrt(N)) + 1):
        if is_prime[p]:
            is_prime[p * p : N + 1 : p] = False

    for p in np.nonzero(is_prime)[0]:
        lp = math.log(float(p))
        q = int(p)
        while q <= N:
            lam[q] = lp
            if q > N // p:
                break
            q *= p

    return lam


def build_cutoffs(N_list: list[int], verbose: bool = True) -> Dict[int, ArithmeticCutoff]:
    """Build ArithmeticCutoff objects for each N."""

    out: Dict[int, ArithmeticCutoff] = {}
    for N in N_list:
        lam = mangoldt_array(N)
        item = ArithmeticCutoff(N=N, mangoldt=lam)
        out[N] = item
        if verbose:
            print(f"  N={N:<8d} psi_minus_N={item.psi_minus_N:+.6e}")
    return out
