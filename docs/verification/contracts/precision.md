---
title: C7 Precision
icon: selah
---

# C7 Precision

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C7 | Contracts | C | CPU (gfortran, nvfortran; `*-r4` modes); GPU (`tests-nvf-oac-r4`) | no | every program of the R4P builds |

## Case and scope

SELAH computes in one real kind per build, `R_P`: double precision by default, single precision with `-D_R_P_IS_R4P`
(the `*-r4` modes, for SELAH and PENF alike). Quadruple precision is refused at compile time. This case says what the
single-precision builds run and assert, and what they cannot.

## Setup

| Program | In R4P |
|---|---|
| kernels, device kernels, storage, tableau forms, stepper options, golden | the same checks; tolerances scale with $\varepsilon(R_P)$; the golden file `selah_golden_r4_v1.txt` |
| schemes, splittings | coefficient checks within 64 ε(R4P); no observed orders (R4P coefficients cannot show the order of the exact ones) |
| stepper | contract and telescoping; no observed orders |
| stability | the polynomial, the closed-form limits and the SSP discs, to the rounding of R4P coefficients |
| forms, TVD | the same bounds (1000 ε) |
| splitting energy | excluded: the high-order errors (about $10^{-8}$) are below R4P round-off over 62 832 steps |
| ADAM parity, ADAM kernels | excluded: ADAM's coefficients and kernels are R8P |

## Program

```bash
fobis build --mode tests-gnu-r4 && for t in exe-r4/selah_*; do $t > /dev/null || echo "FAIL $t"; done
```

## Results

17 of 17 programs pass in R4P with gfortran, the three expected failures included. The golden values are bitwise on
gfortran, and within $4.2\cdot10^{-7}$ with nvfortran ([C4](./golden)).

## Findings

- The stability limits in R4P agree with the R8P ones to $3\cdot10^{-6}$. On the imaginary axis the rounding decides,
  though: rk4, which is tangent to $|R| = 1$ at the origin, crosses it by $10^{-8}$, and euler and ssp-22 gain tiny
  intervals ([A5](../accuracy/stability)).
- The order conditions are evaluated in quadruple precision on the stored coefficients. With nvfortran, which has no
  quad kind, they fall back to R8P, which still resolves R8P coefficients to their own rounding.
