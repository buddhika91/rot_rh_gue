"""Adversarial control families for the ROT-RH finite-GUE benchmark."""
from __future__ import annotations

from typing import Dict

import numpy as np

from .utils import moving_average, zscore

ARITHMETIC_DISRUPTION_MODES = {"permuted", "phase_scramble", "signflip"}
DEFAULT_CONTROL_MODES = ["permuted", "signflip", "gaussian", "phase_scramble"]


def make_control_channels(ch: Dict[str, np.ndarray], mode: str, rng: np.random.Generator) -> Dict[str, np.ndarray]:
    """Return a control channel dictionary.

    Modes:
      * permuted: destroys position/order while preserving samples.
      * signflip: preserves magnitudes, randomizes signs.
      * phase_scramble: preserves approximate power spectrum, randomizes phase.
      * gaussian: generic smooth/universal null; hardest null because GUE spacing
        can arise generically through universality.
    """

    n = len(ch["mass"])
    keys = [
        "mass",
        "mass_dx",
        "mass_smooth",
        "mass_reflect",
        "mass_reflect_dx",
        "mass_reflect_smooth",
        "psi",
        "support",
    ]

    if mode == "permuted":
        out = {"psi_minus_N": ch["psi_minus_N"].copy()}
        for k in keys:
            vv = ch[k].copy()
            rng.shuffle(vv)
            out[k] = zscore(vv)
        return out

    if mode == "signflip":
        signs = rng.choice([-1.0, 1.0], size=n)
        out = {"psi_minus_N": ch["psi_minus_N"].copy()}
        for k in keys:
            out[k] = zscore(ch[k] * signs)
        return out

    if mode == "gaussian":
        m = zscore(rng.normal(size=n))
        mr = m[::-1].copy()
        return {
            "mass": m,
            "mass_dx": zscore(np.gradient(m)),
            "mass_smooth": zscore(moving_average(m, 9)),
            "mass_reflect": zscore(mr),
            "mass_reflect_dx": zscore(np.gradient(mr)),
            "mass_reflect_smooth": zscore(moving_average(mr, 9)),
            "psi": zscore(rng.normal(size=n)),
            "support": zscore(rng.normal(size=n)),
            "psi_minus_N": ch["psi_minus_N"].copy(),
        }

    if mode == "phase_scramble":
        def scramble(y: np.ndarray) -> np.ndarray:
            y = zscore(y)
            Y = np.fft.rfft(y)
            amp = np.abs(Y)
            phase = rng.uniform(0.0, 2.0 * np.pi, len(Y))
            phase[0] = 0.0
            if n % 2 == 0:
                phase[-1] = 0.0
            return zscore(np.fft.irfft(amp * np.exp(1j * phase), n=n))

        out = {"psi_minus_N": ch["psi_minus_N"].copy()}
        for k in keys:
            out[k] = scramble(ch[k])
        return out

    raise ValueError(f"unknown control mode: {mode}")
