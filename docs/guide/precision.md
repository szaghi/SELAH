---
title: Precision
icon: order-conditions
---

# Precision

SELAH computes in one real kind per build: `R_P`, PENF's working kind, re-exported by `use selah`.

| Build | `R_P` | How |
|---|---|---|
| default | R8P (double) | every mode |
| single | R4P | the `*-r4` modes: `-D_R_P_IS_R4P` for SELAH and PENF alike, into their own directories (`libselah_r4.a`, `exe-r4/`) |
| quadruple | refused | `-D_R_P_IS_R16P` stops the build with an `#error` |

A build in which SELAH and PENF disagree on the kind does not compile: `selah_kind_parameters` holds a constant whose
initialisation divides by zero unless SELAH's `R_P` is PENF's. There are no per-kind generics: a program uses one
precision, the one the library was built with.

## What single precision changes

- **The coefficients are rounded to R4P**, the literals and the converted ones alike. The rational coefficients of
  lsrk-54 (numerators near $10^{12}$) are divided in `R_P`.
- **Orders cannot be observed**: the rounded coefficients satisfy the order conditions only to $10^{-7}$, which caps
  the error before the asymptotic order shows. The R4P builds check the coefficients and the contract, not the orders.
- **Some properties hold only to the rounding**: rk4's imaginary-axis stability, for instance, which is a tangency at
  the origin ([A5](/verification/accuracy/stability)).
- **Long-time results below round-off are out of reach**: the energy of the high-order splittings over 62 832 steps
  is not checked in R4P ([P3](/verification/properties/splitting-energy)).

What runs in single precision, and with which tolerances: [C7](/verification/contracts/precision).

## The tests' own precision

The order conditions and the scalar observed orders are evaluated in quadruple precision on the stored coefficients
(`selected_real_kind(33, 4931)`), so the residual measured is the tableau's, not the arithmetic's. nvfortran has no
quad kind, and there they fall back to R8P. The tolerances of the kernels and the stepper scale with $\varepsilon(R_P)$.
