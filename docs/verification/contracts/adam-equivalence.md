---
title: C1 Equivalence with ADAM
icon: selah
---

# C1 Equivalence with ADAM

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C1 | Contracts | C | CPU (gfortran, nvfortran, ifx); GPU (A30) | no | `selah_test_adam_parity`, `selah_test_kernels_adam`, `selah_test_device_adam` |

## Case and scope

SELAH was carved from ADAM's Runge–Kutta layer, and it will be wired back into ADAM. Wiring is safe only if SELAH
computes what ADAM computes, so that ADAM's regression goldens still hold. The oracle is ADAM's own code:
- `src/tests/adam_fnl_rk_kernels.F90`, a byte-exact copy of ADAM's FNL kernels (sha256 `a693ca35…`);
- `src/tests/selah_adam_reference.F90`, ADAM's coefficients;
- `src/tests/selah_adam_cpu_reference.F90`, ADAM's CPU loops;
- `src/bench/adam_flume_rk_residual.F90`, FLUME's step-residual kernel, for the benchmarks.

## Setup

| Program | Compares | Bound |
|---|---|---|
| `selah_test_adam_parity` | the Butcher tableaux of euler, ssp-22, ssp-33, ssp-54 against ADAM's `ssp_tableau`; the two-register coefficients of euler, ssp-22, ssp-33 against ADAM's `ark`, `brk`, `crk` | bitwise |
| `selah_test_kernels_adam` | whole SSP(5,4) and two-register SSP(3,3) steps built from SELAH's box kernels (stage state, step update with the step residual, two-register stage), masked and unmasked: on the FNL layout against ADAM's FNL kernels, on the CPU layout against ADAM's CPU loops | FNL: bitwise; CPU: 4 ε (the summation order differs, ADAM#128 D6) |
| `selah_test_device_adam` | the same on device memory, FNL kernels on the device | bitwise |

## Program

```bash
./exe/selah_test_adam_parity ; ./exe/selah_test_kernels_adam
fobis build --mode tests-nvf-oac --varset iac && ./exe/selah_test_device_adam   # on a GPU node
```

## Results

- Coefficients: bitwise, 11 of 11 checks.
- FNL kernels: bitwise, on gfortran, nvfortran and ifx.
- CPU loops: within 1.00 ε measured, against the 4 ε asserted.
- Device: bitwise on the A30, 18 of 18 checks (medooza jobs 1368099 to 1368177).

## Findings

- SELAH adopts ADAM's FNL summation order everywhere: a private sum, then one store. ADAM's CPU backend sums in another
  order, so ported CPU results move by round-off, and the FNL ones do not move at all
  ([Porting](/cookbook/porting)).
- Equivalence is per kernel. Whole steps differ by construction (SELAH builds stage states from q, ADAM in place), with
  the same mathematical result: see [Performance](/guide/performance) for the traffic and [C2](./forms) for the forms.
