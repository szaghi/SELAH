---
title: A6 Additive splitting
icon: order-conditions
---

# A6 Additive splitting

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A6 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P, the queries in R4P); GPU | no | `selah_test_additive_splitting`, `selah_test_xfail_additive_backward` |

## Case and scope

`additive_splitting` composes a step of the stepper (the field step, $F$) with a flow the caller integrates ($G$)
([guide](/guide/composition)). The case checks the composition's orders with exact sub-flows, with the stepper as the
field step, with both terms time-dependent (the substep times), and with a stiff flow of the radiative-cooling form;
and that a composition with backward substeps is refused unless declared reversible.

## Setup

| Item | Value |
|---|---|
| Problem | $u' = A u + B u$, $A$ the rotation $x' = y$, $y' = -x$, $B$ the decay $x' = -x$, $y' = -2y$ (non-commuting); exact solution $e^{(A+B)t} u_0$ in closed form, $u_0 = (1, 1/2)$, $t = 1$ |
| Exact sub-flows | the rotation and the decay in closed form |
| Field step | the stepper, rk4 and ssp-54 in the Butcher form |
| Time-dependent | $F$ + (sin 2t, 0), $G$ with the rate $(1 + t)$: the stepper's stages at $t_m + c_k h$, the decay integrated exactly from $t_m$ = `substep_time`; reference rk4 on the whole system, 20 000 steps |
| Stiff flow | $u' = a u - k u^2$, $a = 1$: the field step $a u$ (rk4), the exact flow $u/(1 + k h u)$ of $-k u^2$; exact solution the logistic curve; $k u_0$ = 0.01 to 1000 |
| Steps | 5 levels from 4 steps (8 for the stiff case) to $t = 1$ |
| Bound | the last two observed orders whose finer error is above $10^{-12}$ within 0.15 of the order |

## Program

```bash
./exe/selah_test_additive_splitting ; ./exe/selah_test_xfail_additive_backward   # the second must stop
```

## Results

| Composition | Errors at 4 … 64 steps | Last orders |
|---|---|---|
| lie, exact flows (either first) | 2.45e-2 … 1.53e-3 | 1.00, 1.00 |
| strang, flow first / field first | 2.99e-3 … 1.17e-5 / 1.58e-3 … 6.11e-6 | 2.00, 2.00 |
| yoshida-4 (reversible) | 1.03e-4 … 1.63e-9 | 4.00, 4.00 |
| blanes-moan-s6 (reversible) | 1.59e-7 … 2.41e-12 | 4.00, 4.00 |
| blanes-moan-s10 (reversible) | 3.66e-10 … round-off | 6.01, 5.97 (above round-off) |
| lie / strang + the stepper (rk4, ssp-54) | as the exact field flow to 3 digits | 1.00 / 2.00 |
| lie / strang, both terms time-dependent | 6.93e-2 … 4.53e-3 / 9.31e-3 … 3.63e-5 | 0.99, 1.00 / 2.00, 2.00 |
| lie, stiff flow, $k u_0 = 10^{-2} \dots 10^3$ | relative 1.08e-3 … 3.89e-3 at the finest step | 1.00 |
| strang, stiff flow, flow first / field first, $k u_0 = 10^3$ | 1.30e-3 … 5.08e-6 / 6.51e-4 … 2.54e-6 | 2.00, 2.00 |

## Findings

- **No order reduction on the cooling model.** With the stiff flow integrated exactly, Lie and Strang keep orders 1
  and 2 up to $k u_0 \Delta t = 125$, and the error constant saturates as $k$ grows. The literature's order reduction
  in the stiff limit (Verwer & Sportisse 1998) concerns other models; SELAH documents it as possible, not excluded.
- **Each operator carries its own clock.** With both terms time-dependent, substep times from the same operator's
  fractions keep both orders; the stepper's splitting rule (the other partition's fractions) is for sub-flows that are
  Euler steps of one partition.
- The rk4 and ssp-54 field steps give the exact field flow's errors to three digits: the field step's own error
  ($O(h^4)$) is invisible next to the splitting's.
