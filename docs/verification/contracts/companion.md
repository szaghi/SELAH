---
title: C8 Companion states
icon: selah
---

# C8 Companion states

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C8 | Contracts | C | CPU (gfortran, nvfortran, ifx; R8P and R4P, the orders in R8P); GPU | no | `selah_test_companion`, `selah_test_rk_particles` |

## Case and scope

A companion (`initialize_companion`) is a stepper with the field stepper's scheme, form and stages on q, for a state
that advances with the field in the same step, its residual coupled both ways inside each stage: PRISM's particles and
PML auxiliary fields ([guide](/guide/stepper#companion-states)). The case checks that the pair is the coupled system's
Runge–Kutta step, bitwise; that the companion's active extent and component range are the box's, whatever its
capacity and layout; and the properties of particles pushed this way.

## Setup

| Item | Value |
|---|---|
| Coupled system | a field $y$ and 4 particles $(z_1, z_2, c)$: $y' = -y + \mathrm{mean}_i(c_i z_{1,i}^2)$, $z_1' = z_2$, $z_2' = -c(1 + y)z_1$, $c$ constant |
| Forms | every Runge–Kutta scheme in the Butcher form (and with `stages_on_q`), ssp-33 two-register, lsrk-54 2N, ssp-104 Ketcheson |
| One state | the same steps with $y$, $z_1$, $z_2$ as one state of one stepper |
| Extent | capacity 8 particles, 8, 5, 6 active in turn; the companion rebuilt between steps 1 and 2 |
| Layouts | particles as `(3, np)` (components first) and as `(np, 3)` (a scalar rank-2 field) |
| Particles | $x' = v$, $v' = E + v\times B$, $B = \omega\hat z$: gyration ($E = 0$), $E\times B$ drift ($E \perp B$), on a companion of a field stepper |

## Program

```bash
./exe/selah_test_companion ; ./exe/selah_test_rk_particles
```

## Results

| Check | Result |
|---|---|
| field and companion against one state, 8 steps | bitwise, every scheme and form |
| `(np, 3)` against `(3, np)` | bitwise |
| $n$ of 8 particles active against a companion of exactly $n$; the rebuild | bitwise; the particles past $n$ and the constant component untouched |
| the coupled system's order | 1.01, 2.00, 2.99, 3.99, 4.00, 3.99, 3.00, 4.00 (euler … lsrk-54), the same in the low-storage forms |
| gyration: $w_n = R(-i\omega\Delta t)^n w_0$ ($R$ by a complex triangular solve of the tableau) | to $10^{-14}$ |
| stability: $\omega\Delta t = \beta_I(1 - 10^{-6})$ / $\beta_I(1 + 10^{-3})$, $10^4$ steps | no growth / growth (R8P); euler and ssp-22 grow at 0.1 |
| $E\times B$: started at $v_E$ | $v_E$ kept bit for bit (4.7e-15 in the Ketcheson form, whose $(u^n + 9q_1)/25$ rounds); the position advances by $\Delta t\,v_E$ per step |

## Findings

- **A companion needs no code of its own in SELAH**: the stepper's box already holds an active extent, a component
  range and either layout. What is new is the guarantee that the pair is the coupled step (bitwise) and the call that
  binds a companion to its field's scheme.
- **Gyration energy decays in every stable scheme** ($|R(i\omega\Delta t)| < 1$: 0.887 over 50 steps of
  $\omega\Delta t = 0.5$ for ssp-33, 0.998 for ssp-54); euler and ssp-22 gain energy at any step. A two-level pusher
  (Boris) keeps it exactly.
