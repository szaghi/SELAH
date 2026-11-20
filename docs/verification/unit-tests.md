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
| `selah_test_kernels` | every kernel, both layouts, masked and unmasked (and the stage state with `masked_q`), component ranges, active prefixes: values and untouched cells | 348 | [C3](./contracts/stage-contract) |
| `selah_test_kernels_adam` | whole steps from SELAH's kernels against ADAM's FNL kernels (bitwise) and CPU loops (4 ε) | 27 | [C1](./contracts/adam-equivalence) |
| `selah_test_device_kernels` | every kernel on FUNDAL device memory | 57 | [C5](./contracts/device) |
| `selah_test_device_adam` | SELAH's against ADAM's kernels on the device, bitwise | 18 | [C1](./contracts/adam-equivalence), [C5](./contracts/device) |
| `selah_test_stepper` | orders in every form, the stage contract, stages on q, stage coefficients, telescoping, splittings through the stepper | 314 | [A3](./accuracy/stepper-orders), [C3](./contracts/stage-contract), [P4](./properties/conservation) |
| `selah_test_stepper_options` | `stage_time` (the abscissae; a splitting's from the other partition's flows), caller registers, `dq`, masks through the stepper, the mask policy `masked_q` | 49 | [C3](./contracts/stage-contract) |
| `selah_test_forest` | the stepper driven as ADAM's forest drives its integrator: two realms of two blocks, seams filled per stage or lagged (#113), every Runge–Kutta scheme and form: per-stage seams give the single-block solution bitwise; lagged seams the lagged scheme's definition; the reflux register telescopes; end-of-step and per-stage reflux conserve; caller registers are scratch between steps (the D_flux diagnostic) | 220 | [C3](./contracts/stage-contract) |
| `selah_test_nonautonomous` | time-dependent residuals through the stepper, every scheme and form: $u' = q\,t^{q-1}$ exact for $q \le p$ (splittings with either partition time-dependent); observed orders on the Riccati equation and on a forced oscillator | 44 (22 in R4P, without the orders) | [C3](./contracts/stage-contract) |
| `selah_test_stage_euler` | the per-stage forward-Euler fractions and the SSP coefficient against the literature, in every form; planted Shu–Osher forms; a positivity limiter through the stepper: positive with $f_k\,\Delta t$, negative with half of it | 76 | [P2](./properties/ssp-coefficients) |
| `selah_test_additive_splitting` | the field step composed with the caller's flow: queries, orders of lie, strang and the palindromes with exact sub-flows (either operator first), with the stepper as the field step, with both terms time-dependent, with a stiff flow (the cooling model) | 48 (14 in R4P, the queries) | [A6](./accuracy/additive-splitting) |
| `selah_test_particle_flows` | kick and drift through the splitting form: KDK's discrete solution and phase error on the oscillator, Kepler orders and bounded energy for every splitting, gas and sinks coupled through the additive splitting | 16 (6 in R4P, the oscillator) | [A7](./accuracy/particle-flows) |
| `selah_test_companion` | a field and a companion (`initialize_companion`) coupled both ways: bitwise one state in every form, both particle layouts, active extent within a capacity, rebuild, the coupled orders | 53 (42 in R4P, without the orders) | [C8](./contracts/companion) |
| `selah_test_rk_particles` | `imag_axis_bound`, `real_axis_bound` against closed forms and computed values; gyration $w_n = R^n w_0$ on a companion; stability at the bound; the $E\times B$ drift kept | 71 (51 in R4P) | [C8](./contracts/companion), [A5](./accuracy/stability) |
| `selah_test_pushers` | `selah_push`: Boris's angle, speed, drift and layouts; every pusher's defining equation; force-free; relativistic drift; volume; the non-relativistic limit | 23 (19 in R4P, without the volume values) | [A8](./accuracy/pushers) |
| `selah_test_stage_current` | `lookahead` and `stage_current`: Gauss's law ($D + \sin x$) at every stage and step, every scheme; the look-ahead against the next stage; admissibility | 27 | [C9](./contracts/stage-currents) |
| `selah_test_stepper_ranks` | fields of rank 1 to 4, seven shapes: rank k bitwise rank 5, orders, the stage contract, stage coefficients, masks, splittings | 848 (694 in R4P, without the RK orders) | [A3](./accuracy/stepper-orders), [C3](./contracts/stage-contract), [C5](./contracts/device) |
| `selah_test_golden` | 1312 pinned values | 42 | [C4](./contracts/golden) |
| `selah_vv_stability` | stability polynomials, axis limits (and the library's `imag_axis_bound`, `real_axis_bound` against them), SSP discs | 30 | [A5](./accuracy/stability) |
| `selah_vv_ssp_tvd` | TVD up to $\nu = C$ | 5 | [P1](./properties/ssp-tvd) |
| `selah_vv_splitting_energy` | bounded energy error for splittings, drift for Runge–Kutta | 5 | [P3](./properties/splitting-energy) |
| `selah_vv_forms` | every low-storage form on its Butcher form's solution | 6 | [C2](./contracts/forms) |
| `selah_test_xfail_dq_form` | `dq` for a non-Butcher form stops the program | — | [C3](./contracts/stage-contract) |
| `selah_test_xfail_registers` | registers numbered from 0 stop the program | — | [C3](./contracts/stage-contract) |
| `selah_test_xfail_additive_backward` | an additive splitting with backward substeps, not declared reversible, stops the program | — | [A6](./accuracy/additive-splitting) |
| `selah_test_xfail_form` | a scheme in a form it does not have stops the program | — | [C3](./contracts/stage-contract) |

The benchmarks (`src/bench/`, modes `bench-*`) are not tests: they time, they do not assert.
