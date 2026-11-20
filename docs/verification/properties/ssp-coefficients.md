---
title: P2 SSP coefficients
icon: ssp
---

# P2 SSP coefficients

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P2 | Properties | P | CPU (gfortran, nvfortran, ifx; R8P and R4P); the limiter also on the device | no | `selah_test_schemes`, `selah_test_tableau_forms`, `selah_vv_stability`, `selah_test_stage_euler` |

## Case and scope

SELAH never types an SSP coefficient: `butcher(scheme)%ssp_coefficient` is computed from the scheme's Shu–Osher form,
$C = \min_{\beta_{ik} > 0} \alpha_{ik}/\beta_{ik}$, or 0 if any coefficient is negative. The case checks the computed
values against the published optimal ones. It also checks the function on hand-made forms, and the disc property
[A5](../accuracy/stability) derives from $C$.

## Setup

| Scheme | Published $C$ | Source | Effective $C/s$ |
|---|---|---|---|
| `euler`, `ssp-22`, `ssp-33` | 1 | Shu & Osher (1988); Gottlieb & Shu (1998) | 1, 0.5, 0.333 |
| `ssp-54` | 1.508 | Spiteri & Ruuth (2002) | 0.302 |
| `ssp-104` | 6 | Ketcheson (2008) | 0.6 |
| `rk4`, `lsrk-33`, `lsrk-54` | not SSP (0) | | |

Bounds: $|C - 1| \le 64\,\varepsilon$ for euler, ssp-22, ssp-33; $|C - 6| \le 64\,\varepsilon \cdot 6$ for ssp-104;
$|C - 1.508| < 5\cdot10^{-4}$ for ssp-54 (Spiteri & Ruuth publish $C$ to four digits; SELAH computes 1.5082 from their
15-digit coefficients); exactly 0 for rk4 and the 2N schemes.

## Program

```bash
./exe/selah_test_schemes ; ./exe/selah_test_tableau_forms
```

## Results

The values are those of the [stability table](../accuracy/stability) (column $C$), and every one matches.

## Per-stage forward-Euler steps

`stage_euler_fraction(k)` $= f_k = \max_{i>k} \beta_{ik}/\alpha_{ik}$, the forward-Euler step (a fraction of
$\Delta t$) that the residual of stage $k$ enters in the Shu–Osher form: the step a positivity limiter must make
admissible at stage $k$ ([guide](/guide/stepper#limiters-and-forward-euler-steps)). `selah_test_stage_euler` checks:

| Check | Bound |
|---|---|
| $f_k$ = 1 (euler, ssp-22, ssp-33), 1/6 (ssp-104), the row ratios of Spiteri & Ruuth's coefficients (ssp-54: 0.3918, 0.6631, 0.6631, 0.6631, 0.5844) | 100 ε |
| the same values in every storage form of the scheme (two-register, Ketcheson) | 100 ε |
| $\max_k f_k = 1/C$ | 100 ε |
| 0 for rk4, lsrk-33, lsrk-54; a planted form with a negative coefficient, or $\beta_{ik} > 0$ with $\alpha_{ik} = 0$ | exact |
| a stage entering two steps (planted: ratios 1 and 1/2) takes the larger | exact |
| the decay $u_i' = -\lambda_i u_i$, $\lambda_i \Delta t$ = 0.5, 2, 5, 20, the residual of stage $k$ scaled by $\min(1, 1/(\lambda \tau_k))$: every stage state $\ge 0$ with $\tau_k = f_k \Delta t$; some $< -0.1$ with $f_k \Delta t/2$; $\ge 0$ with $\Delta t$ | 100 ε |

Five mutants of the fractions (α/β for β/α; $f_k = 1$, ADAM's choice; the last row dropped; the last ratio instead of
the largest; $C = 1$) each fail the program.

## Findings

- **$C$ is a property of a representation.** The minimum is taken over the Shu–Osher coefficients SELAH stores, and
  another Shu–Osher form of the same Butcher tableau can give a larger value. For the catalogue the stored forms are the
  optimal ones, so $C$ is the scheme's SSP coefficient. For a tableau you build yourself, 0 means "not computed" as well
  as "not SSP".
- SSP(10,4) has the best effective coefficient in the catalogue: per residual evaluation it advances twice as far as
  SSP(5,4) at its strong-stability limit ([Cookbook](/cookbook/schemes#the-largest-strong-stability-preserving-time-step)).
