"""Small utilities used across the finite-GUE codebase."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def parse_int_list(text: str) -> list[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_str_list(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def parse_windows(text: str) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"window must be lo:hi, got {part!r}")
        lo, hi = part.split(":", 1)
        out.append((float(lo), float(hi)))
    return out


def zscore(y: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)
    y = y - float(np.mean(y))
    sd = float(np.std(y))
    if sd < eps:
        return np.zeros_like(y)
    return y / sd


def moving_average(y: np.ndarray, width: int) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)
    width = int(width)
    if width <= 1:
        return y.copy()
    if width % 2 == 0:
        width += 1
    width = min(width, max(3, len(y) // 2 * 2 - 1))
    if width < 3:
        return y.copy()
    k = np.ones(width, dtype=np.float64) / float(width)
    return np.convolve(y, k, mode="same")


def detrend(y: np.ndarray, degree: int) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)
    if degree <= 0 or len(y) < degree + 3:
        return y - float(np.mean(y))
    x = np.linspace(-1.0, 1.0, len(y))
    coef = np.polyfit(x, y, degree)
    return y - np.polyval(coef, x)


def write_csv(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def mean_std(vals: Iterable[float]) -> tuple[float, float]:
    arr = np.asarray(list(vals), dtype=np.float64)
    if len(arr) == 0:
        return float("nan"), float("nan")
    return float(np.mean(arr)), float(np.std(arr))
