---
title: A5 Linear stability
---

# A5 Linear stability

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A5 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P and R4P); GPU (A30, `tests-nvf-oac`) | yes | `selah_vv_stability` |

## Case and scope

For $y' = \lambda y$ a step multiplies $y$ by $R(z)$, $z = \lambda\Delta t$. For an explicit tableau $R$ is the
polynomial $R(z) = 1 + \sum_{j=1}^{s} b^T A^{j-1}\mathbf{1}\,z^j$. The case computes $R$ from `butcher(scheme)` itself,
so it verifies the library's coefficients, not a textbook's.

What the case catches:
- **a coefficient that breaks the order**: the first $p+1$ coefficients of $R$ must be $1/j!$;
- **a wrong scheme**: the limits on the axes must match the closed forms where they are known;
- **a wrong SSP coefficient**: an SSP scheme with coefficient $C$ must be stable on the disc $|z + C| \le C$.

## Setup

| Item | Value |
|---|---|
| Arithmetic | quadruple precision (R8P with nvfortran); the polynomial by Horner |
| Limits | bisection to $|R| = 1$ along the negative real and the imaginary axis, from 0 |
| Closed forms | euler, ssp-22: $-2$ and 0; ssp-33: $-2.512745$ and $\sqrt 3$; rk4: $-2.785294$ and $2\sqrt 2$ |
| SSP disc | 720 points of the circle $|z + C| = C$ (the maximum principle covers its inside); real limit $\le -2C$ |
| Tolerance | coefficients and limits: $\max(10^{-9}, 1000\,\varepsilon(R_P))$; on the imaginary axis $|R| \le 1 + 64\,\varepsilon(R_P)$, the rounding of the coefficients; "no interval" (euler, ssp-22) allows the one this rounding opens, up to $(8 \cdot 64\,\varepsilon)^{1/4}$ |

## Program

```bash
./exe/selah_vv_stability
SELAH_VV_DIR=exe/vv-data ./exe/selah_vv_stability     # stability-<scheme>.dat (|R| on a grid), stability-summary.dat
```

## Results

![Stability regions](/verification/stability.png)

*The boundary $|R(z)| = 1$ of every scheme. Left: the SSP schemes; ssp-104's region is the largest by far. Right: rk4
and the 2N schemes.*

<!--@include: ../_tables/stability.md-->

## Findings

- Every polynomial agrees with $e^z$ up to $z^p$, every closed-form limit matches, and every SSP scheme is stable on
  its disc.
- lsrk-33 has the stability polynomial of ssp-33: three stages of third order leave no freedom.
- **In single precision the imaginary axis is a matter of rounding.** rk4 is tangent to $|R| = 1$ at the origin, and its
  R4P coefficients cross that line by about $10^{-8}$: strictly, no interval is left. Conversely, euler and ssp-22, which
  have no interval, gain one at the rounding scale (0.004 and 0.088 in R4P). Measured with the tolerance of the
  coefficients' rounding, the R4P limits agree with the R8P ones to $3\cdot10^{-6}$ ([C7](../contracts/precision)).
