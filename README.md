# ROT RH / Reflected Mangoldt–Stieltjes Finite-GUE Operator

Project Status

This repository is still under active construction. The current ROT–RH–GUE shows promising structure and passes several adversarial/control tests, but it does not yet fully defeat Gaussian universality. That remains the main open benchmark I am working on.

Until that is resolved, researchers and reviewers can use this repository to test the rest of the framework, including the operator construction, prime/Mangoldt structure, S-fraction recursion, control-family comparisons, and finite GUE-spacing diagnostics.

---------------------------------------------------------------------------------------------------------------------------------------------------------

A reproducible Python codebase for constructing and falsifying a finite-dimensional arithmetic Jacobi operator whose local eigenvalue spacing matches the Gaussian Unitary Ensemble (GUE) benchmark across blind Riemann-prime scales.

> **Status:** finite-GUE computational benchmark candidate.  
> **Not claimed:** a proof of the Riemann Hypothesis, an asymptotic theorem, or a Hilbert–Pólya construction for the zeta zeros.

---

## 1. Executive summary

This repository implements a frozen arithmetic operator discovered through a sequence of adversarial tests. The construction begins with the Mangoldt function, builds a direct-binned arithmetic field, reflects that field, converts the two fields into positive Stieltjes measures, extracts Jacobi recurrence coefficients, and tests the resulting finite Jacobi matrix against GUE spacing.

The core discovery is that the one-channel Stieltjes operator is not stable enough, while the reflected two-channel operator is. The reflected pair

```math
F(x),\qquad F(-x)
```

acts as a self-dual arithmetic source. Each channel generates a positive measure. The two Stieltjes recurrences are then mixed into one finite Jacobi operator.

The strongest frozen endpoint audit so far gives the following summary for the geometric reflected lock:

```text
LOCK_bgeom_best_W012_088_L13
N = 20000,24000,28000,32000,36000,40000,44000,48000,56000,64000,80000,96000
d = 64,96,128,160,192,224
controls = 32 per control family seed panel

pooled controls:
  minPos    = 11/12
  meanPos   = 11.83/12
  minMinKSz = -0.035
  meanKSz   = +1.363

arithmetic-disrupted controls only:
  minPos    = 12/12
  meanPos   = 12.00/12
  minKSz    = +0.218
  meanKSz   = +2.159
```

Thus, the current honest claim is:

> We found a frozen reflected two-channel Mangoldt–Stieltjes Jacobi operator that produces robust finite-GUE spacing up to `N=96000`. It fully passes the arithmetic-disruption benchmark and remains a pooled-control near-lock. Gaussian controls remain a competitive universality baseline.

---

## 2. Mathematical object being tested

### 2.1 Mangoldt source

Let

```math
\Lambda(n)=
\begin{cases}
\log p, & n=p^k,\ k\ge 1,\\
0, & \text{otherwise}.
\end{cases}
```

For a cutoff `N`, the code forms a direct-binned arithmetic defect field from the cumulative prime-power mass. The field is not built from zeta zeros. It is built directly from the prime-side Mangoldt data.

A binning map sends integer mass on `1 <= n <= N` into a finite support grid

```math
x_j\in[-1,1].
```

The default support uses a cosine/Chebyshev-like grid

```math
x_j=-\cos\left(\pi\frac{j}{M-1}\right),\qquad j=0,\dots,M-1.
```

The raw mass in bin `j` is compared to the expected bin length. After detrending and normalization this gives the arithmetic field

```math
F_N(x_j).
```

In code this is implemented by direct binning, polynomial detrending, smoothing, and z-score normalization.

---

## 3. Reflected two-channel Stieltjes construction

The key empirical discovery is that a one-channel field

```math
F_N(x)
```

is not stable enough. The global GUE signal stabilizes when the arithmetic field is paired with its reflected dual

```math
F_N(-x).
```

This creates two positive measures:

```math
\mu_+(x)=\exp(A F_N(x))(1-x^2)^\beta dx,
```

```math
\mu_-(x)=\exp(A F_N(-x))(1-x^2)^\beta dx.
```

The best current frozen geometric lock uses

```math
A=2.20,\qquad \beta=1.35.
```

The secondary coefficient-mix lock uses

```math
A=2.20,\qquad \beta=1.15.
```

These are finite-dimensional numerical measures. They are used to generate orthogonal polynomial recurrence coefficients through a Stieltjes/Jacobi recurrence.

---

## 4. Stieltjes recurrence and Jacobi matrix

For each positive measure, compute the orthonormal polynomial recurrence

```math
xq_n(x)=b_{n+1}q_{n+1}(x)+a_n q_n(x)+b_n q_{n-1}(x).
```

This gives two coefficient systems:

```math
\mu_+\longrightarrow (a_+,b_+),
```

```math
\mu_-\longrightarrow (a_-,b_-).
```

The final operator is a finite Jacobi matrix

```math
J_d=
\begin{pmatrix}
a_0 & b_0 & 0 & \cdots & 0\\
b_0 & a_1 & b_1 & \cdots & 0\\
0 & b_1 & a_2 & \ddots & 0\\
\vdots & \vdots & \ddots & \ddots & b_{d-2}\\
0 & 0 & 0 & b_{d-2} & a_{d-1}
\end{pmatrix}.
```

The tested dimensions are typically

```math
d\in\{64,96,128,160,192,224\}.
```

---

## 5. Frozen finalist operators

This repository freezes two finalists. They are not searched inside the final audit.

### 5.1 Geometric reflected lock

This is the current best theory candidate.

```text
lock key      : LOCK_bgeom_best_W012_088_L13
operator mode : bgeom_mix
window        : [0.12, 0.88]
local width   : 13
A             : 2.20
beta          : 1.35
gamma         : 0.55
```

The diagonal mixing is

```math
a=a_+-0.55a_-.
```

The off-diagonal geometric mixing is

```math
b=\exp\left(0.45\log b_+ + 0.55\log b_-\right).
```

The geometric off-diagonal mix is important because Jacobi off-diagonals are positive scale variables. Geometric mixing treats them multiplicatively instead of additively.

### 5.2 Coefficient reflected lock

```text
lock key      : LOCK_coeff_ba115_sm025_g045_W020_080_L7
operator mode : coeff_mix
window        : [0.20, 0.80]
local width   : 7
A             : 2.20
beta          : 1.15
gamma         : 0.45
```

The diagonal mixing is

```math
a=a_+-0.45a_-.
```

The off-diagonal coefficient mixing is

```math
b=0.55b_+ + 0.45b_-.
```

This lock is useful because it is structurally different from the geometric lock, yet it survives many of the same tests.

---

## 6. GUE spacing benchmark

For each finite Jacobi matrix `J_d`, the eigenvalues are computed:

```math
\lambda_1\le\lambda_2\le\cdots\le\lambda_d.
```

A bulk window is selected, for example

```math
[0.12,0.88]
```

or

```math
[0.20,0.80].
```

Nearest-neighbor spacings are formed and locally unfolded:

```math
s_i = \frac{\lambda_{i+1}-\lambda_i}{\text{local mean spacing near }i}.
```

After normalization to mean one, the empirical spacing distribution is compared with the Wigner surmise for GUE:

```math
P_{\mathrm{GUE}}(s)=\frac{32}{\pi^2}s^2e^{-4s^2/\pi}.
```

The corresponding CDF used for the Kolmogorov-Smirnov comparison is

```math
F_{\mathrm{GUE}}(s)=\operatorname{erf}\left(\frac{2s}{\sqrt\pi}\right)-\frac{4s}{\pi}e^{-4s^2/\pi}.
```

For a given real operator and a family of null controls, the code computes

```math
KSz=\frac{\operatorname{mean}(KS_{control})-KS_{real}}{\operatorname{std}(KS_{control})}.
```

Thus:

```text
KSz > 0  means the real arithmetic operator is closer to GUE than controls.
KSz < 0  means controls are closer to GUE than the real operator.
```

---

## 7. Controls and falsification strategy

The benchmark is adversarial. The real arithmetic field is compared against several nulls.

### 7.1 Permuted control

Randomly permutes the arithmetic field values. This preserves the empirical amplitude distribution but destroys arithmetic order.

Purpose:

```text
Tests whether the signal is just the histogram of field values.
```

### 7.2 Phase-scramble control

Preserves approximate Fourier amplitudes while randomizing phases.

Purpose:

```text
Tests whether the signal is only the power spectrum, not the arithmetic phase/order.
```

### 7.3 Signflip control

Randomly flips signs across the field.

Purpose:

```text
Tests sensitivity to coherent sign structure.
```

### 7.4 Gaussian control

Replaces the arithmetic source with Gaussian fields.

Purpose:

```text
Tests against generic random-field universality.
```

This is the hardest control family. Gaussian random fields can naturally generate GUE-like spacing because GUE statistics are universal. Therefore, the Gaussian benchmark should be reported separately from the arithmetic-disruption benchmark.

---

## 8. Current result summary

For the endpoint audit up to `N=96000` with `32` controls, the geometric reflected lock gives:

```text
LOCK_bgeom_best_W012_088_L13

pooled controls:
  minPos    = 11/12
  meanPos   = 11.83/12
  minMinKSz = -0.035
  meanKSz   = +1.363

arithmetic-disrupted controls:
  minPos    = 12/12
  meanPos   = 12.00/12
  minKSz    = +0.218
  meanKSz   = +2.159

phase_scramble only:
  minPos    = 12/12
  meanPos   = 12.00/12
  minKSz    = +0.467
  meanKSz   = +2.477

permuted only:
  minPos    = 12/12
  meanPos   = 12.00/12
  minKSz    = +0.312
  meanKSz   = +2.273

gaussian only:
  minPos    = 6/12
  meanPos   = 8.50/12
  minKSz    = -1.247
  meanKSz   = +0.399
```

The single pooled negative cell is

```text
d = 96
N = 56000
KSz = -0.035
```

This is a near-tie, not a collapse. The same cell is positive against permutation and phase-scramble controls, but negative against Gaussian controls.

---

## 9. Interpretation of the result

The result should be interpreted in two layers.

### 9.1 Arithmetic-disruption benchmark

The operator fully passes the arithmetic-disruption benchmark.

```text
The reflected Mangoldt-Stieltjes operator is closer to GUE than controls
that destroy arithmetic order or arithmetic phase.
```

This supports the claim that the arithmetic ordering and reflected dual structure carry genuine finite-GUE information.

### 9.2 Gaussian-universality benchmark

The operator does not fully defeat Gaussian universality.

```text
Gaussian controls remain competitive and sometimes more GUE-like.
```

This does not falsify the arithmetic signal. It means the final claim must distinguish arithmetic-disruption controls from generic random-field universality controls.

---

## 10. Falsification history

This project did not begin with the current frozen lock. The current operator survived after several failures.

### 10.1 One-channel Stieltjes failure

A one-channel Stieltjes operator using only `F(x)` showed local promise but did not remain stable under grid changes and blind-scale extension.

Failure mode:

```text
one-channel architecture lacked global stability.
```

### 10.2 Direct Stieltjes/S-fraction search failure

Large one-channel candidate searches did not produce a robust global finite-GUE lock.

Failure mode:

```text
many candidates were locally GUE-like but failed blind N-values or dimensions.
```

### 10.3 Recursive-memory and gated-memory failure

Recursive memory and gated S-fraction variants generated interesting raw spectra but did not beat the reflected two-channel operator under adversarial controls.

Failure mode:

```text
extra recursion did not produce a stable global benchmark advantage.
```

### 10.4 Mean-unfolding artifact

An earlier reflected two-channel candidate appeared to fail at

```text
N = 32000, d = 128.
```

Bulk-window and local-unfolding audits showed that this was a measurement artifact. Local unfolding removed the apparent defect.

### 10.5 d=96 endpoint resonance

Endpoint extension exposed a weak line at

```text
d = 96.
```

A local-width audit showed that the issue was sensitive to the unfolding width. Correcting local width rescued the endpoint behavior and produced full positivity against arithmetic-disruption controls.

### 10.6 Gaussian-control near-failure

The final pooled benchmark retains one near-zero failure due to Gaussian controls. This is now treated as a separate universality null rather than hidden inside the arithmetic-disruption result.

---

## 11. How to run

Install:

```bash
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
pip install -e .
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e .
```

Run the frozen endpoint audit:

```powershell
python scripts/run_frozen_audit.py `
  --N-list 20000,24000,28000,32000,36000,40000,44000,48000,56000,64000,80000,96000 `
  --jacobi-dims 64,96,128,160,192,224 `
  --controls 32 `
  --seed 2027 `
  --out-prefix reflect_frozen_endpoint_96000_32_seed2027
```

Run the control-family scoreboard:

```powershell
python scripts/control_family_scoreboard.py `
  --prefix reflect_frozen_endpoint_96000_32_seed2027
```

Expected output pattern for the best lock:

```text
LOCK_bgeom_best_W012_088_L13
arithmetic_disrupted   minPos=12/12 meanPos=12.00/12 minKSz>0
aII_controls           near-lock, usually 11/12 with one near-zero Gaussian-driven cell
```

---

## 12. Repository layout

```text
src/rot_rh_gue/
  arithmetic.py      # Mangoldt function and prime-power source
  channels.py        # direct-binned arithmetic fields and reflected fields
  controls.py        # permuted, signflip, gaussian, phase-scramble nulls
  stieltjes.py       # Stieltjes/Jacobi recurrence
  operator.py        # frozen reflected Jacobi operators
  gue.py             # GUE spacing and KS metric
  audit.py           # frozen audit runner
  utils.py           # parsing and CSV helpers

scripts/
  run_frozen_audit.py
  control_family_scoreboard.py

docs/
  OPERATOR_STRUCTURE.md
  PROOF_PROTOCOL.md
  FALSIFICATION_HISTORY.md
```

---

## 13. What would falsify the current claim?

The current claim would be weakened or falsified if any of the following happen reproducibly:

1. The arithmetic-disruption benchmark falls below full positivity under larger seeds and controls.
2. Endpoint extension to larger `N` produces systematic negative cells across many dimensions.
3. The reflected operator loses its advantage under independent implementations.
4. The same performance appears for non-arithmetic fake fields under the arithmetic-disruption controls.
5. The local-unfolding rule is shown to be overfitted or dependent on post-hoc tuning.

The repository therefore freezes the locks and separates discovery scripts from confirmation scripts.

---

## 14. What is not proven

This project does not prove:

```text
RH
Hilbert-Polya
asymptotic GUE for zeta zeros
uniqueness of the operator
```

The result is finite and computational:

```text
A frozen reflected Mangoldt-Stieltjes Jacobi operator produces robust
finite-GUE spacing across tested blind scales and survives arithmetic-disruption
controls.
```

---

## 15. Research meaning

The discovery suggests that the prime-side Mangoldt field contains a reflected Stieltjes structure that can generate finite-GUE-like spectra without directly using zeta zeros. The most important architectural features are:

```text
1. direct arithmetic source Λ(n)
2. reflected duality F(x), F(-x)
3. positive Stieltjes measures
4. Jacobi recurrence coefficients
5. geometric or coefficient mixing
6. local unfolding
7. adversarial control separation
```

The strongest conceptual point is the reflection principle:

```math
F(x)\longleftrightarrow F(-x).
```

In this benchmark, reflection is the difference between unstable local signals and robust finite-GUE behavior.

---

## 16. Suggested citation language

If referencing this repository, use cautious language:

```text
This repository reports a finite-dimensional computational benchmark in which a
frozen reflected Mangoldt-Stieltjes Jacobi operator reproduces GUE-like local
spacing across tested prime cutoffs. The result is robust against arithmetic-
disruption controls, while Gaussian universality remains a competitive null.
No proof of RH or asymptotic theorem is claimed.
```

---

## 17. License

MIT License, unless otherwise specified.
