---
title: A2 Orders on scalar ODEs
icon: order-conditions
---

# A2 Orders on scalar ODEs

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A2 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P) | yes | `selah_test_schemes` |

## Case and scope

The order conditions hold for the coefficients; this case checks the order of the *integrator*. Every scheme
integrates two problems to $t = 1$ with 8, 16, 32, 64 and 128 steps, and the observed order between the last two must
reach the declared one. The nonlinear problem catches a scheme that is only of order $p$ on linear problems: its
linear order can exceed its order.

## Setup

| Item | Value |
|---|---|
| Problems | $y' = -y$, $y(0) = 1$ (exact $e^{-t}$); $y' = -y^2$, $y(0) = 1$ (exact $1/(1+t)$) |
| Integration | the Butcher tableau, in quadruple precision (`integrate_scalar`) |
| Bound | observed order of the last halving at least $p - 0.2$ |
| Precision | R8P builds only: R4P coefficients cannot show the order of the exact ones |

## Program

```bash
./exe/selah_test_schemes
SELAH_VV_DIR=exe/vv-data ./exe/selah_test_schemes      # scalar-orders.dat
```

## Results

![Scalar orders, nonlinear problem](/verification/scalar-orders-2.png)

*The error at $t = 1$ of every scheme on $y' = -y^2$. The linear problem: [figure](/verification/scalar-orders-1.png).*

<!--@include: ../_tables/scalar-orders.md-->

## Findings

- Every scheme reaches its declared order on both problems.
