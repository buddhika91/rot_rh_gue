# Operator Structure

This document describes the frozen operator implemented in `src/rot_rh_gue`.

## 1. Arithmetic input

The only arithmetic input is the von Mangoldt array:

```math
\Lambda(n)=\begin{cases}\log p,&n=p^k,\\0,&\text{otherwise.}\end{cases}
```

For each cutoff `N`, the code records:

```math
\psi(N)-N=\sum_{n\le N}\Lambda(n)-N.
```

## 2. Direct-binned arithmetic field

For each support grid size `M`, the interval `[1,N]` is divided into `M` bins. The direct arithmetic defect field is:

```math
F_j=\sum_{n\in B_j}\Lambda(n)-|B_j|.
```

The code detrends and z-scores this field. It also builds:

```text
F(x)
F'(x) proxy
smoothed F(x)
F(-x)
F'(-x) proxy
smoothed F(-x)
```

## 3. Two-channel reflection

The essential discovery is that the successful operator needs both channels:

```math
F(x),\quad F(-x).
```

The plus and reflected measures are:

```math
\mu_+(x)=\exp(A(F(x)+dF'(x)+s\,\mathrm{smooth}(F)(x)))(1-x^2)^\beta dx,
```

```math
\mu_-(x)=\exp(A(F(-x)+dF'(-x)+s\,\mathrm{smooth}(F)(-x)))(1-x^2)^\beta dx.
```

The current `bgeom` frozen lock uses:

```text
A          = 2.20
beta       = 1.35
dx_mix     = 0.10
smooth_mix = 0.25
```

## 4. Stieltjes recurrence / S-fraction

Each measure is converted into Jacobi recurrence coefficients by orthogonalizing polynomials:

```math
xq_n=b_nq_{n+1}+a_nq_n+b_{n-1}q_{n-1}.
```

This is the finite Stieltjes/S-fraction recursion stage:

```text
mu_+ -> (a_+, b_+)
mu_- -> (a_-, b_-)
```

## 5. Frozen mixing laws

Two frozen finalist operators are implemented.

### Coefficient-mix finalist

```text
LOCK_coeff_ba115_sm025_g045_W020_080_L7
```

```math
a=0.04\,z(a_+-0.45a_-),
```

```math
b=0.55b_+ + 0.45b_-.
```

Measurement:

```text
window      = [0.20, 0.80]
local_width = 7
```

### Geometric-b finalist

```text
LOCK_bgeom_best_W012_088_L13
```

```math
a=0.04\,z(a_+-0.55a_-),
```

```math
b=\exp(0.45\log b_+ + 0.55\log b_-).
```

Measurement:

```text
window      = [0.12, 0.88]
local_width = 13
```

This is the main endpoint-to-96000 candidate.

## 6. Spectral test

For each Jacobi dimension `d`, build:

```math
J_d=\operatorname{tridiag}(b,a,b).
```

The eigenvalue spacings are measured in the bulk window, locally unfolded, and compared to the GUE Wigner-surmise CDF using a KS distance.

The reported score is:

```math
KSz=\frac{\mathbb{E}[KS_\mathrm{control}]-KS_\mathrm{real}}{\mathrm{Std}(KS_\mathrm{control})}.
```

Positive `KSz` means the real arithmetic operator is closer to GUE than the selected control family.
