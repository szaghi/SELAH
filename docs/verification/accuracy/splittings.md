---
title: A4 Splittings
---

# A4 Splittings

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A4 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P) | yes | `selah_test_splittings`, `selah_test_stepper` |

## Case and scope

A splitting composes flows that each advance one partition of the components. Its order depends on the composition's
coefficients and on the sub-flows being exact. The case checks:
- the compositions themselves: each partition sums to one step, the sequence is a palindrome, the partitions
  alternate;
- the observed order, with either partition first, on two separable Hamiltonian problems whose sub-flows are exact;
- the same orders through the stepper (`selah_test_stepper`, the oscillator).

## Setup

| Item | Value |
|---|---|
| Problems | the harmonic oscillator $q' = p$, $p' = -q$ (exact); the pendulum $q' = p$, $p' = -\sin q$ (reference: classic RK4 with 40 000 steps in quadruple precision) |
| Sub-flows | exact: $q \mathrel{+}= c\,\Delta t\, p$ with $p$ frozen, $p \mathrel{-}= c\,\Delta t\, V'(q)$ with $q$ frozen |
| Steps | 4 halvings from 8 steps per unit time (2 for the order-6 Blanes–Moan S10, which reaches the R8P floor sooner) |
| Bound | observed order of the last halving at least $p - 0.2$ |

## Program

```bash
./exe/selah_test_splittings ; ./exe/selah_test_stepper
SELAH_VV_DIR=exe/vv-data ./exe/selah_test_splittings   # splitting-orders.dat
```

## Results

![Splitting orders](/verification/splitting-orders.png)

*The pendulum, partition 1 first.*

<!--@include: ../_tables/splitting-orders.md-->

## Findings

- Strang 2, Yoshida 4, Blanes–Moan S6 4 and S10 6, with either partition first, on both problems.
- The stepper's sub-flow is a forward-Euler step of one partition, $y_p \mathrel{+}= c_k\Delta t\,R_p$. It is the exact
  flow of that partition only when $R_p$ does not depend on partition $p$. That holds for the separable Hamiltonians
  here, not for a general splitting, where the composition falls to first order
  ([Tutorial 6](/tutorial/06-splittings)).
