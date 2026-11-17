---
title: C2 Storage forms
icon: selah
---

# C2 Storage forms

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C2 | Contracts | C | CPU (gfortran, nvfortran, ifx; R8P and R4P) | yes | `selah_vv_forms` |

## Case and scope

A scheme with a low-storage form computes the same step in two registers. The two-register, Williamson 2N and
Ketcheson forms store other coefficients and add the stage terms in another order. Over a run they must stay on the
Butcher form's solution to round-off: otherwise a coefficient of the low-storage form, or the form's recipe in the
stepper, is wrong.

## Setup

| Item | Value |
|---|---|
| Pairs | euler, ssp-22, ssp-33: two-register; lsrk-33, lsrk-54: Williamson 2N; ssp-104: Ketcheson |
| Problem | periodic advection, 64 cells, central residual, CFL 0.5, a sine plus a jump |
| Steps | 40, through the stepper |
| Bound | $\max |q_{form} - q_{Butcher}| \le 1000\,\varepsilon \max|q|$ at every step |

## Program

```bash
./exe/selah_vv_forms
SELAH_VV_DIR=exe/vv-data ./exe/selah_vv_forms      # forms.dat
```

## Results

![Forms](/verification/forms.png)

*The relative difference between each low-storage form and the Butcher form of the same scheme, step by step.*

<!--@include: ../_tables/forms.md-->

## Findings

- Every low-storage form stays within $10^{-14}$ of its Butcher form. The Ketcheson form, with 10 stages, is the
  furthest at $8\cdot10^{-15}$. euler's two forms are bitwise identical.
- The golden file pins every form's values separately ([C4](./golden)), so a change that kept the forms equal but
  moved both is caught there.
