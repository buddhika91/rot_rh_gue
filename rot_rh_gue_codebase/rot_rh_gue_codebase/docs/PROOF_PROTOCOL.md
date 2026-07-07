# Finite-GUE Proof Protocol

This is a finite numerical benchmark protocol, not a mathematical proof of RH.

## Objects fixed before confirmation

A confirmation run is valid only when these are frozen before the run:

```text
operator formula
operator parameters
N-list
Jacobi dimensions
bulk window
local-unfolding width
control modes
control count
random seed
```

Parameter searches and frozen confirmations must be reported separately.

## The finite claim

For a fixed set of target cutoffs `N_1,...,N_K` and dimensions `d_1,...,d_m`, a control family is passed when:

```text
minPos = K/K
minMinKSz > 0
```

where `minPos` is the minimum, over dimensions, of the number of positive `KSz` target cutoffs.

A near-lock is:

```text
minPos >= (K-1)/K
minMinKSz > -0.25
meanKSz > +1.2
```

## Control families

The repository separates control families because they test different hypotheses.

### Permuted

Destroys spatial arithmetic order while preserving the sample distribution.

### Phase scramble

Preserves approximate power spectrum but destroys phase coherence.

### Sign flip

Preserves magnitudes but randomizes signs.

### Gaussian

Generic random universality baseline. This is the hardest null because local GUE-like spacing can arise generically in random matrix-like constructions.

### Arithmetic-disrupted pooled control

Defined as:

```text
permuted + phase_scramble + signflip
```

This is the best null for asking whether the arithmetic ordering/reflection structure matters.

## Current result summary

For the endpoint-to-96000 frozen run, the geometric-b lock reported:

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
```

Interpretation:

```text
Arithmetic-disruption benchmark: passed.
Pooled Gaussian-including benchmark: near-lock.
```

## What would falsify the claim?

Any of the following would weaken or falsify the current finite benchmark claim:

1. A frozen rerun with the same settings and larger controls produces strongly negative `minMinKSz`.
2. Endpoint extension causes many dimensions to fail at new N values.
3. Arithmetic-disruption controls become competitive across many cells.
4. The result vanishes when the reflected channel `F(-x)` is removed.
5. The result depends on changing parameters after seeing the target output.

## What would strengthen it?

1. A 64-control or 128-control frozen endpoint confirmation remains near-lock or full lock.
2. Extension to larger endpoints, e.g. `N=128000,160000`, remains stable.
3. The same operator works across more Jacobi dimensions.
4. A theoretical explanation is found for why reflection plus Stieltjes recursion carries GUE structure.
