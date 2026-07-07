"""Stieltjes recurrence and Jacobi matrix construction."""
from __future__ import annotations

import math

import numpy as np

from .utils import zscore


def cosine_support(size: int) -> np.ndarray:
    """Cosine support x_j in [-1,1].

    The successful operators used a direct-binned measure on this compact
    support.  The edge factor (1-x^2)^β gives the measure a semicircle-like
    bias while leaving the arithmetic driver in the exponential.
    """

    u = np.linspace(0.0, 1.0, int(size))
    return -np.cos(math.pi * u)


def stieltjes_coefficients(x: np.ndarray, w: np.ndarray, dim: int) -> tuple[np.ndarray, np.ndarray]:
    """Return Jacobi recurrence coefficients (a,b) for a positive measure.

    The measure is discrete: sum_j w_j δ_{x_j}.  The recurrence builds
    orthonormal polynomials q_n satisfying

        x q_n = b_n q_{n+1} + a_n q_n + b_{n-1} q_{n-1}.

    This is the numerical S-fraction / Stieltjes-recursion step: the measure's
    Stieltjes transform has a continued fraction with coefficients a_n,b_n.
    """

    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    w = np.maximum(w, 1e-15)
    w = w / np.sum(w)

    q_prev = np.zeros_like(x)
    q = np.ones_like(x)
    q = q / math.sqrt(float(np.sum(w * q * q)))

    a = np.zeros(dim, dtype=np.float64)
    b = np.zeros(dim - 1, dtype=np.float64)
    beta_prev = 0.0

    for n in range(dim):
        r = x * q - beta_prev * q_prev
        alpha = float(np.sum(w * r * q))
        r = r - alpha * q

        # Light reorthogonalization keeps the small-dimensional finite audit
        # stable without invoking a heavy arbitrary-precision package.
        r = r - float(np.sum(w * r * q)) * q
        if n > 0:
            r = r - float(np.sum(w * r * q_prev)) * q_prev

        beta = math.sqrt(max(0.0, float(np.sum(w * r * r))))
        a[n] = alpha

        if n < dim - 1:
            if beta < 1e-13:
                beta = 1e-13
            b[n] = beta
            q_prev, q = q, r / beta
            beta_prev = beta

    return a, b


def normalize_offdiagonal(b: np.ndarray, dim: int, edge_taper: float, offdiag_scale: float) -> np.ndarray:
    """Normalize and taper off-diagonal Jacobi coefficients."""

    bb = np.maximum(np.asarray(b, dtype=np.float64), 1e-14)
    bb = bb / (np.mean(bb) + 1e-14)

    if edge_taper != 0.0:
        k = np.arange(1, dim, dtype=np.float64)
        taper = np.sqrt(k * (dim - k)) / (dim / 2)
        taper = np.maximum(taper, 1e-8)
        bb = bb * (taper ** edge_taper)

    return offdiag_scale * bb


def jacobi_eigenvalues(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Eigenvalues of a real symmetric tridiagonal Jacobi matrix."""

    J = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
    return np.linalg.eigvalsh(J)


def centered_diagonal_mix(ap: np.ndarray, am: np.ndarray, gamma: float, diag_scale: float) -> np.ndarray:
    """Mix plus/reflected diagonals using the frozen convention."""

    return diag_scale * zscore(ap - gamma * am)
