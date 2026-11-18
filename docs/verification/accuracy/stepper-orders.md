---
title: A3 Orders through the stepper
icon: order-conditions
---

# A3 Orders through the stepper

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A3 | Accuracy | P | CPU (gfortran, nvfortran, ifx); GPU (A30, `tests-nvf-oac`) | yes | `selah_test_stepper`, `selah_test_stepper_ranks` |

## Case and scope

The same order, now through everything a solver uses: the stepper, its storage forms, the box kernels, a multi-block
field and a caller that fills ghost cells and evaluates the residual. Every scheme runs in every form it has: Butcher,
Butcher with stages on q, two-register, Williamson 2N and Ketcheson.

## Setup

| Item | Value |
|---|---|
| Problem | $u_t + u_x = 0$, periodic on $[0, 1]$, $u_0 = \sin 2\pi x$ |
| Grid | 4 blocks of 16 cells, 1 ghost cell, layout V_FIRST; FUNDAL arrays (device memory in a device build) |
| Residual | central flux form, $R_i = -(u_{i+1} - u_{i-1})/(2\Delta x)$, written by the caller on the device |
| Reference | the exact semi-discrete solution $\sin 2\pi(x - ct)$, $c = \sin(2\pi\Delta x)/(2\pi\Delta x)$: the error is the time integrator's alone |
| Steps | to $t = 0.25$; 16, 32, 64, 128 steps (8, …, 64 for order 4) |
| Bound | observed order of the last halving at least $p - 0.2$ |

### Fields of rank 1 to 4

`selah_test_stepper_ranks` takes the same orders through the stepper's rank-k fronts
([Stepper](/guide/stepper#fields-of-rank-1-to-4)).

| Item | Value |
|---|---|
| Fields | seven shapes: scalar `(i)`, `(i, j)`, `(i, j, k)`, `(i, j, k, b)` (two blocks, a `box_range` box), and two components `(v, i)`, `(v, i, j)`, `(v, i, j, k)`; 3 interior cells per dimension, 1 ghost cell; a pointer of the field's rank over contiguous FUNDAL storage |
| Problem | $y' = -y^2$ cell by cell (no space operator), $y_0 = 1 + 0.01\,n$; exact $1/(1/y_0 + t)$ |
| Steps | to $t = 0.5$; 16, 32, 64, 128 steps (8, …, 64 for order 4) |
| Splittings | the oscillator $q' = p$, $p' = -q$ cell by cell on the shapes with components, to $t = 1$ |
| Bound | observed order of the last halving at least $p - 0.2$ (R8P builds) |

## Program

```bash
./exe/selah_test_stepper
SELAH_VV_DIR=exe/vv-data ./exe/selah_test_stepper      # stepper-orders.dat, telescoping.dat
./exe/selah_test_stepper_ranks                         # the rank-k fronts
```

## Results

![Orders through the stepper](/verification/stepper-orders.png)

*Error at $t = 0.25$ against the number of steps, Butcher form. Every form of every scheme is in the table.*

<!--@include: ../_tables/stepper-orders.md-->

## Findings

- Every scheme reaches its order in every form, and the forms of a scheme agree in their errors to round-off
  ([C2](../contracts/forms)).
- `stages_on_q` changes no digit ([C3](../contracts/stage-contract)).
- At every rank from 1 to 4, every scheme reaches its order in every form, and every splitting its order on the
  shapes with components. The rank-k run is bitwise the rank-5 run of the same data.
- The same program checks the stage contract and the conservation identity (C3, [P4](../properties/conservation)), on
  the host and on the A30.
