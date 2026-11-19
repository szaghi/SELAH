---
title: Test programs
---

# Test programs

`bash scripts/run_tests.sh` runs every program in `exe/`, built by a `tests-*` mode of FoBiS (`--exe-dir exe-r4` for
the single-precision `tests-*-r4` modes). Each program exits
non-zero on any failed check, and prints one line per check. The `*_xfail_*` programs must fail: the runner passes them
when they exit with an error. The suite runs:
- with gfortran 16.2 (optimised, debug with every runtime check, OpenMP with 4 threads), nvfortran 26.9 (optimised,
  debug) and ifx 2026.2;
- under fpm;
- in single precision ([C7](./contracts/precision));
- on an NVIDIA A30 ([C5](./contracts/device)).

| Program | What it asserts | Checks (R8P) | Page |
|---|---|---|---|
| `selah_test_kinds` | the kinds reach a client through `use selah` | 5 | — |
| `selah_test_fundal` | FUNDAL allocation, bounds, initial values; a device is present in a device build | 3 | — |
| `selah_test_tableau_forms` | the conversions between coefficient forms; the order checker; `is_explicit`; the SSP coefficient on hand-made forms | 12 | [P2](./properties/ssp-coefficients) |
| `selah_test_schemes` | order conditions, exact order, $c = A\mathbf 1$, explicitness, SSP coefficients, Shu–Osher consistency, form equivalence, published values, observed orders | 78 | [A1](./accuracy/order-conditions), [A2](./accuracy/scalar-orders), [P2](./properties/ssp-coefficients) |
| `selah_test_adam_parity` | ADAM's tableaux and two-register coefficients, bitwise | 11 | [C1](./contracts/adam-equivalence) |
| `selah_test_splittings` | partition sums, palindromes, alternation, observed orders on two problems | 37 | [A4](./accuracy/splittings) |
| `selah_test_storage` | `registers_needed` for every scheme and form | 23 | — |
| `selah_test_kernels` | every kernel, both layouts, masked and unmasked, component ranges, active prefixes: values and untouched cells | 316 | [C3](./contracts/stage-contract) |
| `selah_test_kernels_adam` | whole steps from SELAH's kernels against ADAM's FNL kernels (bitwise) and CPU loops (4 ε) | 27 | [C1](./contracts/adam-equivalence) |
| `selah_test_device_kernels` | every kernel on FUNDAL device memory | 57 | [C5](./contracts/device) |
| `selah_test_device_adam` | SELAH's against ADAM's kernels on the device, bitwise | 18 | [C1](./contracts/adam-equivalence), [C5](./contracts/device) |
| `selah_test_stepper` | orders in every form, the stage contract, stages on q, stage coefficients, telescoping, splittings through the stepper | 314 | [A3](./accuracy/stepper-orders), [C3](./contracts/stage-contract), [P4](./properties/conservation) |
| `selah_test_stepper_options` | `stage_time` (the abscissae; a splitting's from the other partition's flows), caller registers, `dq`, masks through the stepper | 24 | [C3](./contracts/stage-contract) |
| `selah_test_forest` | the stepper driven as ADAM's forest drives its integrator: two realms of two blocks, seams filled per stage or lagged (#113), every Runge–Kutta scheme and form: per-stage seams give the single-block solution bitwise; lagged seams the lagged scheme's definition; the reflux register telescopes; end-of-step and per-stage reflux conserve; caller registers are scratch between steps (the D_flux diagnostic) | 220 | [C3](./contracts/stage-contract) |
| `selah_test_nonautonomous` | time-dependent residuals through the stepper, every scheme and form: $u' = q\,t^{q-1}$ exact for $q \le p$ (splittings with either partition time-dependent); observed orders on the Riccati equation and on a forced oscillator | 44 (22 in R4P, without the orders) | [C3](./contracts/stage-contract) |
| `selah_test_stepper_ranks` | fields of rank 1 to 4, seven shapes: rank k bitwise rank 5, orders, the stage contract, stage coefficients, masks, splittings | 848 (694 in R4P, without the RK orders) | [A3](./accuracy/stepper-orders), [C3](./contracts/stage-contract), [C5](./contracts/device) |
| `selah_test_golden` | 1312 pinned values | 42 | [C4](./contracts/golden) |
| `selah_vv_stability` | stability polynomials, axis limits, SSP discs | 22 | [A5](./accuracy/stability) |
| `selah_vv_ssp_tvd` | TVD up to $\nu = C$ | 5 | [P1](./properties/ssp-tvd) |
| `selah_vv_splitting_energy` | bounded energy error for splittings, drift for Runge–Kutta | 5 | [P3](./properties/splitting-energy) |
| `selah_vv_forms` | every low-storage form on its Butcher form's solution | 6 | [C2](./contracts/forms) |
| `selah_test_xfail_dq_form` | `dq` for a non-Butcher form stops the program | — | [C3](./contracts/stage-contract) |
| `selah_test_xfail_registers` | registers numbered from 0 stop the program | — | [C3](./contracts/stage-contract) |
| `selah_test_xfail_form` | a scheme in a form it does not have stops the program | — | [C3](./contracts/stage-contract) |

The benchmarks (`src/bench/`, modes `bench-*`) are not tests: they time, they do not assert.
