# Falsification History

This file records the major attacks that shaped the current operator.

## 1. Global phase transport

Early phase-transport operators produced real correlations, but the signal was not stable globally. The dilation/log-periodic branch discovered a horizon-like behavior, but it was not itself the final finite-GUE operator. Boundary/clipping audits were necessary because dilation can create artifacts.

Outcome:

```text
useful diagnostic, not final operator
```

## 2. One-channel Stieltjes operators

One-channel Stieltjes/S-fraction Jacobi generators produced promising local rows but did not produce a robust global lock across dimensions and blind targets.

Outcome:

```text
failed as final global operator
```

## 3. Recursive-memory / gated-memory operators

Gated recursive-memory and two-channel recursive-memory scans produced candidates with decent raw KS but did not beat controls robustly enough across the full grid.

Outcome:

```text
useful search direction, not the final lock
```

## 4. Two-channel reflection

Adding the reflected channel `F(-x)` was the major structural improvement. It stabilized the operator and produced a strong finite-GUE candidate.

Outcome:

```text
critical discovery
```

## 5. Fixed mean unfolding artifact

A previous obstruction at `N=32000,d=128` disappeared after switching to local unfolding.

Outcome:

```text
not an operator failure; measurement artifact
```

## 6. d=96 endpoint resonance

Endpoint extension to `N=64000` exposed a `d=96` resonance. Scanning local width/window showed the resonance could be resolved without changing the operator family.

Outcome:

```text
measurement-scale mismatch; led to frozen local_width choices
```

## 7. Frozen endpoint confirmation

Frozen endpoint confirmation to `N=96000` showed:

```text
bgeom lock vs arithmetic-disrupted controls: 12/12 full lock
bgeom lock vs all pooled controls: 11/12 near-lock, minKSz=-0.035
```

The one remaining pooled weakness was Gaussian-driven, not permutation/phase-scramble driven.

Outcome:

```text
finite arithmetic-disruption benchmark passed;
Gaussian universality remains competitive
```
