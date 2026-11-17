---
title: P4 Conservation
icon: ssp
---

# P4 Conservation

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P4 | Properties | P | CPU (gfortran, nvfortran, ifx); GPU (A30, `tests-nvf-oac`) | yes | `selah_test_stepper` |

## Case and scope

A conservative finite-volume solver changes the mass only through its boundary fluxes. With a Runge–Kutta step, the
step's flux is the weighted sum of the stage fluxes:

$$
m^{n+1} - m^n = -\Delta t \sum_k b_k \left(F_{\text{right}} - F_{\text{left}}\right)_k .
$$

An AMR flux register accumulates exactly this sum at coarse–fine faces, using `stage_weight(k)`. The case checks the
identity in every Runge–Kutta form. That includes the low-storage ones, whose coefficients are not $b_k$, so their
`stage_weight` is derived from the Butcher form. As a control it repeats the sum with uniform weights $1/s$.

## Setup

| Item | Value |
|---|---|
| Problem | $u_t + u_x = 0$, inflow $u = 1$ on the left, outflow on the right, $u_0 = 1 + \frac12\sin 2\pi x$ |
| Grid | 4 blocks of 16 cells, central flux form; FUNDAL arrays |
| Steps | 4, every scheme in every form (Butcher, two-register, 2N, Ketcheson) |
| Bound | $|\Delta m + \Delta t \sum_k w_k F_k| \le 256\,\varepsilon \max(|m|, 1)$ at every step |
| Control | the same with $w_k = 1/s$ (reported, not asserted) |

## Program

```bash
./exe/selah_test_stepper
SELAH_VV_DIR=exe/vv-data ./exe/selah_test_stepper      # telescoping.dat
```

## Results

![Telescoping](/verification/telescoping.png)

*The conservation defect with `stage_weight` (coloured) and with uniform weights (grey), largest over 4 steps.*

<!--@include: ../_tables/telescoping.md-->

## Findings

- The identity holds to round-off in every form of every scheme.
- Uniform weights miss by up to $3.5\cdot10^{-2}$ (ssp-33): a refluxed AMR solver with the wrong weights silently
  stops conserving. For euler, ssp-22 and ssp-104 the weights *are* uniform, and both columns agree.
- The [tutorial](/tutorial/08-flux-register) and the [AMR cookbook](/cookbook/amr) show the register in use, and a
  correction at a face with mismatched fluxes restoring conservation to $2\cdot10^{-16}$.
