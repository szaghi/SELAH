---
title: C4 Golden regression
---

# C4 Golden regression

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C4 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised; R8P and R4P) | no | `selah_test_golden` |

## Case and scope

The verification cases assert orders, properties and bounds. A change of a few ulp, or an error that grows inside its
bound, passes them. The golden test pins the values themselves: what every scheme, form and kernel computes on fixed
inputs.

## Setup

| Item | Value |
|---|---|
| Cases | 41: every Runge–Kutta scheme in every form it has (3 steps of a periodic linear system, 32 values); every splitting (3 steps of a partitioned linear system); every kernel, box V_FIRST unmasked and V_LAST masked, and flat |
| Inputs | exact in binary: dyadic values $((37i + 11s) \bmod 61 - 30)/32$ and dyadic coefficients and steps (no libm, no random numbers), so the inputs are the same bits on every compiler and option |
| Files | `src/tests/golden/selah_golden_v1.txt` (R8P) and `selah_golden_r4_v1.txt` (R4P), 1312 values each, written by gfortran 16.2 `-O3` (compiler and options in the header) |
| Tolerance | $1000\,\varepsilon(R_P)$ relative to $\max(|v|, 1)$; every computed case must be in the file and every line of the file computed |
| Rule | a golden file is never edited. `selah_test_golden write FILE` refuses to overwrite one. A deliberate change of results gets a new file (v2) and a new default path in the test |

## Program

```bash
./exe/selah_test_golden                                    # from the source root, or with SELAH_SRC set
./exe/selah_test_golden write src/tests/golden/selah_golden_v2.txt   # only for a deliberate change, new file
```

## Results

Measured on 2026-11-17:

| Build | Bitwise equal | Worst relative difference |
|---|---|---|
| gfortran 16.2, optimised (the writer) and debug | 1312 / 1312 | 0 |
| ifx 2026.2, optimised and debug | 1312 / 1312 | 0 |
| nvfortran 26.9, optimised and debug | 1246 / 1312 | $6.7\cdot10^{-16}$ |
| gfortran 16.2, R4P (the writer of the R4P file) | 1312 / 1312 | 0 |
| nvfortran 26.9, R4P | 1232 / 1312 | $4.2\cdot10^{-7}$ |

## Findings

- gfortran and ifx agree bit for bit at every optimisation level: the kernels' private-sum-then-store order leaves the
  compilers nothing to reorder. nvfortran contracts some products into fused multiply–adds, a rounding per affected
  value.
- In the mutation run ([C6](./mutation)) the golden test alone kills 21 of the 28 mutants. The seven it misses change
  what its cases do not record:
  - the SSP coefficient, the stage weight, the register count and the stages-on-q flags;
  - a mask branch its masked case does not take;
  - a 13th-digit change of a 2N coefficient;
  - the unread Shu–Osher coefficient.

  It is the cheapest net, and the least explanatory one: the case pages say *what* broke.
