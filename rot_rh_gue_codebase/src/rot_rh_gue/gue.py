"""GUE spacing measurement and local unfolding."""
from __future__ import annotations

import math

import numpy as np

from .utils import moving_average


def gue_wigner_cdf(s: np.ndarray) -> np.ndarray:
    """Approximate GUE nearest-neighbor spacing CDF using Wigner surmise.

    The GUE surmise density is p(s)=32/pi^2 s^2 exp(-4s^2/pi), normalized to
    mean spacing 1.  Its CDF is

        erf(2s/sqrt(pi)) - (4s/pi) exp(-4s^2/pi).
    """

    s = np.asarray(s, dtype=np.float64)
    return np.array([
        math.erf(2 * x / math.sqrt(math.pi)) - (4 * x / math.pi) * math.exp(-4 * x * x / math.pi)
        for x in s
    ])


def ks_gue(spacings: np.ndarray) -> float:
    """One-sample KS distance from GUE Wigner surmise."""

    s = np.asarray(spacings, dtype=np.float64)
    s = s[np.isfinite(s)]
    s = s[s >= 0]
    if len(s) < 8:
        return float("nan")

    s = np.sort(s)
    n = len(s)
    emp = np.arange(1, n + 1, dtype=np.float64) / n
    target = gue_wigner_cdf(s)
    return float(np.max(np.abs(emp - target)))


def local_unfolded_spacing_ks(vals: np.ndarray, lo: float, hi: float, local_width: int) -> float:
    """Bulk-window local-unfolded GUE KS distance.

    vals are Jacobi eigenvalues.  We restrict to a bulk window [lo,hi], take
    spacings, divide by a local moving mean, renormalize to mean one, and test
    against the GUE Wigner surmise.
    """

    vals = np.asarray(vals, dtype=np.float64)
    d = len(vals)
    a = max(0, min(int(lo * d), d - 3))
    b = max(a + 10, min(int(hi * d), d))

    sp = np.diff(vals[a:b])
    sp = sp[np.isfinite(sp)]
    sp = sp[sp > 0]
    if len(sp) < 8:
        return float("nan")

    loc = moving_average(sp, local_width)
    loc = np.where(loc <= 1e-12, np.mean(sp), loc)
    sp = sp / loc
    sp = sp / np.mean(sp)
    return ks_gue(sp)
