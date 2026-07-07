# ROT-RH Finite-GUE Reflected Stieltjes Operator

This repository contains a reproducible Python codebase for the current ROT-RH finite-GUE benchmark branch.

It implements the frozen reflected two-channel Mangoldt-Stieltjes Jacobi operators, the GUE spacing measurement, the adversarial controls, and the falsification/audit protocol that led to the current result.

## Current honest status

This is **not** a proof of the Riemann Hypothesis and it is **not** an asymptotic theorem.

The current claim is narrower and computational:

> A frozen reflected two-channel Mangoldt-Stieltjes Jacobi operator produces robust finite-GUE spacing on tested blind cutoffs up to `N = 96000`. Against arithmetic-disruption controls it achieves a full lock; against pooled controls including Gaussian universality it is a near-lock.

Best frozen operator:

```text
LOCK_bgeom_best_W012_088_L13
mode         = bgeom_mix
window       = [0.12, 0.88]
local_width  = 13
A            = 2.20
beta         = 1.35
gamma        = 0.55
```

Empirical scoreboard from the frozen endpoint audit:

```text
arithmetic_disrupted controls:
  minPos  = 12/12
  meanPos = 12.00/12
  minKSz  = +0.218
  meanKSz = +2.159

all pooled controls:
  minPos  = 11/12
  meanPos = 11.83/12
  minKSz  = -0.035
  meanKSz = +1.363

gaussian controls:
  minPos  = 6/12
  meanPos = 8.50/12
  minKSz  = -1.247
  meanKSz = +0.399
```

So the clean conclusion is:

```text
Finite-GUE arithmetic-disruption benchmark: passed.
Full Gaussian-universality benchmark: near-lock, not fully passed.
```

## Operator structure

The arithmetic field is built from the von Mangoldt function `Λ(n)`. For a cutoff `N`, direct bins define a defect field:

```text
F_j = sum_{n in bin j} Λ(n) - bin_width_j.
```

The reflected two-channel Stieltjes measures are:

```math
\mu_+(x)=\exp(A F(x))(1-x^2)^\beta dx
```

```math
\mu_-(x)=\exp(A F(-x))(1-x^2)^\beta dx
```

A Stieltjes recurrence converts each measure into Jacobi coefficients:

```text
mu_+ -> (a_+, b_+)
mu_- -> (a_-, b_-)
```

For the strongest current lock:

```math
a = a_+ - \gamma a_-
```

```math
b = \exp((1-\gamma)\log b_+ + \gamma \log b_-)
```

with:

```text
A     = 2.20
beta  = 1.35
gamma = 0.55
```

Then the finite Jacobi matrix is:

```math
J_d = \operatorname{tridiag}(b,a,b).
```

Its eigenvalue spacings are locally unfolded in a bulk window and compared to the GUE Wigner-surmise spacing CDF.

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e .
```

This project only requires NumPy.

## Run the frozen audit

Quick smoke run:

```powershell
python scripts/run_frozen_audit.py `
  --N-list 20000,24000 `
  --jacobi-dims 64,96 `
  --controls 1 `
  --out-prefix smoke_frozen_audit
```

Main endpoint run:

```powershell
python scripts/run_frozen_audit.py `
  --N-list 20000,24000,28000,32000,36000,40000,44000,48000,56000,64000,80000,96000 `
  --jacobi-dims 64,96,128,160,192,224 `
  --controls 32 `
  --seed 2027 `
  --out-prefix reflect_frozen_endpoint_96000_32_seed2027
```

Then run the control-family separation:

```powershell
python scripts/control_family_scoreboard.py `
  --prefix reflect_frozen_endpoint_96000_32_seed2027
```

## Outputs

The audit writes:

```text
<prefix>_real_rows.csv
<prefix>_control_rows.csv
<prefix>_summary.csv
<prefix>_leaderboard_by_dim.csv
<prefix>_aggregate.csv
```

The control-family analyzer writes:

```text
<prefix>_control_family_family_summary.csv
<prefix>_control_family_cell_rows.csv
```

## What counts as success?

For a finite benchmark with `K` target N-values:

```text
minPos = K/K
minMinKSz > 0
```

is a full finite lock against the selected control family.

A near-lock is:

```text
minPos >= (K-1)/K
minMinKSz > -0.25
meanKSz > +1.2
```

For the current endpoint-to-96000 run, the best operator is a full lock against arithmetic-disruption controls and a pooled-control near-lock.

## Repository layout

```text
src/rot_rh_gue/
  arithmetic.py   # von Mangoldt and ψ(N)-N
  channels.py     # direct-binned arithmetic channels F(x), F(-x)
  stieltjes.py    # Stieltjes recurrence / Jacobi coefficients
  operator.py     # frozen reflected two-channel operators
  gue.py          # local unfolding and GUE KS distance
  controls.py     # adversarial control families
  audit.py        # frozen audit runner

scripts/
  run_frozen_audit.py
  control_family_scoreboard.py

docs/
  OPERATOR_STRUCTURE.md
  PROOF_PROTOCOL.md
  FALSIFICATION_HISTORY.md
```

## Caveats

1. This is a finite numerical benchmark, not an RH proof.
2. Gaussian controls remain a competitive universality baseline.
3. The strongest defensible claim is about arithmetic-disruption controls and frozen endpoint near-locks.
4. Any new parameter search must be reported separately from frozen confirmation runs.
