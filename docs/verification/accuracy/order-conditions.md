---
title: A1 Order conditions
icon: order-conditions
---

# A1 Order conditions

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A1 | Accuracy | P | CPU (gfortran, nvfortran, ifx; debug and optimised; R8P and R4P) | yes | `selah_test_schemes` |

## Case and scope

A Runge–Kutta scheme has order $p$ when its Butcher coefficients satisfy one algebraic condition per rooted tree of up
to $p$ nodes. SELAH's coefficients come from several sources:
- typed literals (euler, ssp-22, ssp-33, rk4);
- Shu–Osher coefficients converted at run time (ssp-54, ssp-104);
- 2N coefficients converted at run time (lsrk-33, lsrk-54).

A digit lost in any of them shows up here first.

What the case catches:
- **a wrong coefficient**: a residual far above round-off, up to the declared order;
- **a scheme better than declared**: when $p < 4$, some condition of order $p+1$ must fail. Otherwise the declared
  order is not the scheme's order;
- **an inconsistent tableau**: $c = A\mathbf{1}$, and A strictly lower triangular (explicit).

## Setup

| Item | Value |
|---|---|
| Conditions | the 8 trees of 1–4 nodes: $\sum b = 1$, $\sum bc = 1/2$, $\sum bc^2 = 1/3$, $\sum bAc = 1/6$, $\sum bc^3 = 1/4$, $\sum bcAc = 1/8$, $\sum bAc^2 = 1/12$, $\sum bAAc = 1/24$ |
| Arithmetic | the stored coefficients converted to quadruple precision (`selected_real_kind(33)`; R8P with nvfortran, which has none) |
| Bound | 64 ε(R_P) on every condition up to $p$; at least one condition of order $p+1$ above $10^{-6}$ |
| Code | `order_residuals` in [`selah_test_tools.F90`](https://github.com/szaghi/SELAH/blob/main/src/tests/selah_test_tools.F90) |

## Program

```bash
fobis build --mode tests-gnu && ./exe/selah_test_schemes
SELAH_VV_DIR=exe/vv-data ./exe/selah_test_schemes      # report mode: order-conditions.dat
```

## Results

![Order conditions](/verification/order-conditions.png)

*The largest residual of the conditions up to the declared order (coloured), and of the first conditions of the next
order (grey).*

<!--@include: ../_tables/order-conditions.md-->

## Findings

- Every scheme satisfies its conditions to round-off of the stored coefficients. The coefficients converted at run
  time (ssp-54, ssp-104, the 2N schemes) are within a few ε, the literal ones exact.
- No scheme of order below 4 is accidentally of higher order.
- The same program also checks the SSP coefficients ([P2](../properties/ssp-coefficients)), the equivalence of the
  coefficient forms, the consistency of the Shu–Osher forms (each row of α sums to 1) and the published values.
