# FOODIE audit: the integrators and tests SELAH does not have

*2026-11-18, night foodie-audit. FOODIE at `6f4140b7` (`~/fortran/FOODIE`), SELAH at `c4e5c58`. Every coefficient
claim below was checked by `notes/foodie-audit-check.py`, which parses the coefficients from FOODIE's Fortran (nothing
transcribed by hand) and checks them in 60-digit arithmetic. Its output is `notes/foodie-audit-check.txt`. Statements
about code cite `file:line` of `src/lib/` (integrators) or `src/tests/` (tests).*

## Summary

- FOODIE implements **13 families, 86 named schemes** (`supported_schemes_` of each `foodie_integrator_*.F90`). SELAH
  covers 7 of those names, which are 5 distinct schemes: euler (three FOODIE names), SSP(2,2), SSP(3,3), SSP(5,4) and
  LSRK(5,4). SELAH also has four kinds of scheme FOODIE
  lacks: rk4, LSRK(3,3), SSP(10,4) with its Ketcheson form, and the four splittings.
- **Every coefficient table is right to the precision it is typed in.** All explicit RK tables meet their declared
  order, and all AB/AM/BDF/LMM-SSP tables have their exact order. There is no default-real literal and no integer
  division in any coefficient block.
- **Two FOODIE defects are confirmed; port neither scheme as it is.**
  - The variable-step LMM-SSP schemes do not take the step they are given.
  - The linear-SSP RK `order_s` variant overflows a default-integer factorial at 13 stages.
- **One SELAH test gap is exposed: no non-autonomous problem.** FOODIE's convergence problem is the Riccati equation
  U_t = a t U². No SELAH test evaluates a time-dependent residual, so `stage_time(k)` has never been checked through an
  observed order.
- **Recommended first port:** the five 2N low-storage schemes (6, 7, 12, 13 and 14 stages, order 4). They need only
  coefficients, since `FORM_WILLIAMSON_2N` exists, and they run in 2 registers with optimised stability regions.

## 1. Integrators

`Fits` describes how a port fits SELAH's design (P1–P7 of `notes/PLAN.md`):
- **coefficients**: an existing form and kernel suffice;
- **new form**: needs an L2 form, and possibly a kernel;
- **out of scope**: excluded by PLAN §0.

| Family (file) | FOODIE schemes | Order | Type | Source (as cited in the code) | In SELAH | Fits |
|---|---|---|---|---|---|---|
| explicit Euler (`euler_explicit.F90:29`) | `euler_explicit` | 1 | explicit 1-step | — | yes (`euler`) | — |
| RK-SSP (`runge_kutta_ssp.F90:115`) | `_stages_1_order_1`, `_2_2`, `_3_3`, `_5_4` | 1, 2, 3, 4 | explicit RK | Butcher 1963; Gottlieb, Ketcheson & Shu. The (5,4) coefficients come from NodePy, refined to 80 digits (`:317`) | yes (`euler`, `ssp-22`, `ssp-33`, `ssp-54`) | — |
| RK low storage 2N (`runge_kutta_low_storage.F90:166`) | `_stages_1_order_1`, `_5_order_4` | 1, 4 | explicit RK, 2N | Williamson 1980; Carpenter & Kennedy 1994 (solution 3) | yes (`euler`, `lsrk-54`) | — |
| | `_stages_6_order_4`, `_7_order_4` | 4 | explicit RK, 2N | Allampalli et al., HALE-RK (aeroacoustics); the published 12 digits refined to satisfy the conditions exactly (`:57`, `:71`, ref `:149`) | **no** | coefficients (`FORM_WILLIAMSON_2N`) |
| | `_stages_12_order_4`, `_13_order_4`, `_14_order_4` | 4 | explicit RK, 2N | Niegemann, Diehl & Busch 2012, optimised stability regions (`:86`–`:123`, ref `:152`) | **no** | coefficients (`FORM_WILLIAMSON_2N`) |
| RK embedded (`runge_kutta_embedded.F90:280`) | `_stages_2_order_2` | 2(1) | explicit RK + estimator | Heun–Euler (`:62`) | **no** | coefficients + an error-estimate kernel |
| | `_stages_6_order_5` | 5(4) | | Cash & Karp 1990 (`:72`) | **no** | same |
| | `_stages_7_order_4` | **5**(4), see D3 | | Dormand & Prince 1980 (`:86`) | **no** | same |
| | `_stages_9_order_6` | 6(5) | | Calvo, Montijano & Rández 1990 (`:102`) | **no** | same |
| | `_stages_17_order_10` | 10(8) | | Feagin (`:132`) | **no** | same (17 stages: low value for CFD) |
| RK linear SSP (`runge_kutta_lssp.F90:83`) | `_stages_s_order_s_1`, `_stages_s_order_s`, any s | s−1, s (linear autonomous problems only) | explicit RK, s stages | Gottlieb–Shu linear SSP | **no** | new form (an accumulating stage chain) |
| Leapfrog (`leapfrog.F90:51`) | `leapfrog`, `leapfrog_raw` | 2 | explicit 2-step, Robert–Asselin–Williams filter | Robert 1966; Asselin 1972; Williams (RAW) (`:32`–`:37`) | **no** (deferred at M4) | new form (two time levels) |
| Adams–Bashforth (`adams_bashforth.F90:41`) | `_1` … `_16` | = steps | explicit LMM | — | **no** | new form (residual history) |
| Adams–Moulton (`adams_moulton.F90:43`) | `_0` … `_15` | steps + 1 | implicit LMM, fixed-point iterations | — | **no** | out of scope (implicit) |
| Adams–Bashforth–Moulton (`adams_bashforth_moulton.F90:93`) | `_1` … `_16` | = steps | predictor–corrector (AB_k, AM_{k−1}, fixed-point iterations) | — | **no** | new form; effectively explicit PE(CE)^m |
| BDF (`backward_differentiation_formula.F90:48`) | `back_df_1` … `_6` | = steps | implicit LMM | — | **no** | out of scope (implicit) |
| LMM-SSP (`lmm_ssp.F90:40`) | `_steps_3_order_2`, `_4_3`, `_5_3` | 2, 3, 3 | explicit SSP LMM, 1 residual per step | Gottlieb, Ketcheson & Shu (`:26`) | **no** | new form (solution + residual history) |
| LMM-SSP variable step (`lmm_ssp_vss.F90:52`) | `_steps_2_order_2`, `_3_2`, `_3_3`, `_4_3`, `_5_3` | 2, 2, 3, 3, 3 | explicit SSP LMM, variable dt | Hadjimichael, Ketcheson, Lóczi & Németh, SINUM 54(5) 2016 (`:38`) | **no** | new form; **port the published formula, not FOODIE's (D1)** |
| Multistep–multistage RK-SSP (`ms_runge_kutta_ssp.F90:52`) | `_steps_2_stages_2_order_3`, `_3_2_3`, `_4_5_8` | 3, 3, 8 | explicit MSRK-SSP | Bresten, Gottlieb, Grant, Higgs, Ketcheson & Németh, Math. Comp. (cited as 2016, `:39`) | **no** | new form (history of u and F); cache the residuals (D4) |

Not in FOODIE, in SELAH: `rk4`, `lsrk-33` (Williamson), `ssp-104` (Ketcheson SSPRK(10,4) and its 2-register form), the
splittings `strang`, `yoshida-4`, `blanes-moan-s6` and `blanes-moan-s10`. FOODIE has no splitting family.

### Coefficient verification (`notes/foodie-audit-check.txt`)

| Family | Result |
|---|---|
| RK-SSP 1–4, LS 1/5/6/7 | the declared order exactly. Residuals are at the typed precision: 0 to 1e-26, the tables being typed as rationals or refined to 25+ digits. They fail at declared + 1, so the order is not underestimated. c = A·1 holds |
| LS 12/13/14 (NDB) | order 4, residuals 1e-16 to 1e-15: typed as doubles, so not recomputable in quad precision from FOODIE. Re-derive them from the paper before porting (SELAH's gate is ≤ 1e-14; these pass it, narrowly) |
| embedded, advancing row `beta(:,1)` (`runge_kutta_embedded.F90:387`, `U = self%U1`) / estimator row `beta(:,2)` | Heun–Euler 2/1, Cash–Karp 5/4, Dormand–Prince **5**/4, Calvo 6/5 (typed as doubles, 1e-15), Feagin 10/8 (21 digits, 1e-21) |
| AB 1–16, AM 0–15, BDF 1–6 | the exact order on polynomials: k, k + 1, k |
| LMM-SSP (3,2), (4,3), (5,3) | exact order; a ≥ 0, b ≥ 0, Σa = 1 (SSP form) |
| LMM-SSP-VSS | at constant dt: the declared order. With varying dt: see D1 |
| MSRK-SSP (2,2,3), (3,2,3) | order 3 (2.98, 2.99 at 256 steps) |
| MSRK-SSP (4,5,8) | 6.43, 7.10, 7.69 on y' = −4y² (16 to 128 steps), then a floor at 7e-15 set by its 16–19 typed digits: consistent with order 8, not proven at this precision |
| linear SSP | `order_s` reaches order s up to 12 stages; `order_s_1` reaches s − 1 up to 14 stages. Beyond: D2 |

## 2. Defects in FOODIE

**D1. The variable-step LMM-SSP schemes do not use the step they are given.** Confirmed.
- **The code:** `integrate_order_2`/`_3` (`lmm_ssp_vss.F90:248-267`) take the increment from `self%Dt(self%steps)`. That
  is the *previous* step, stored by `update_previous` (`foodie_integrator_multistep_object.F90:189`). They also build
  ω from `dt_ratio` = `Dt(s)/Dt(last)` over the previous steps only (`lmm_ssp_vss.F90:317-337`). The `Dt` argument reaches only
  `update_previous`, which stamps the new time t + Dt.
- **The effect:** at constant dt this is the right formula. When dt changes, the step advances by dt_{n−1} with ratios
  off by one index, while recording dt_n.
- **The test:** y' = −y², with dt cycling through ×1, ×1.5, ×2/3.
  - FOODIE's formula: observed orders 0.3 to 1.6, irregular.
  - The published formula, with ω_j = dt_{n−j}/dt_n for j = 1..k−1 and increment dt_n: the declared 2 and 3.
  - `steps_3_order_3` diverges in both, because this pattern violates its ω ≥ 2 requirement. That is a limit of the
    test, not of FOODIE.
- **Why FOODIE's suite misses it:** no FOODIE test drives a variable-step scheme with a varying dt. The suite's
  `integrate` takes one `Dt` per run (`tests/suite/foodie_test_convergence_order.f90:60`).

**D2. Linear SSP RK, `order_s` variant, is not consistent beyond 12 stages.** Confirmed.
- **The code:** `initialize_order_s` (`runge_kutta_lssp.F90:300`) computes `alpha(i) = 1._R_P / product([(j, j=1,i)])`
  as a default-integer factorial, which overflows int32 at 13! = 6 227 020 800.
- **The effect:** gfortran 16.2 wraps it to 1 932 053 504, and `-ftrapv` aborts. With the wrapped value R(z) − e^z =
  O(z), so the scheme is no longer consistent.
- **Why FOODIE's suite misses it:** it runs this family with `STAGES = 3` only (`tests/suite/foodie_test_suite_utils.f90:33`).
- **The fix:** compute the factorial in `R_P`.

**D3. Dormand–Prince is named `order_4` and advances at order 5.** It advances with the 5th-order row, which is local
extrapolation, as its comment says (`runge_kutta_embedded.F90:86`). FOODIE's suite knows about it and expects 5
(`foodie_test_suite_utils.f90:117`). A naming defect, not a numerical one.

**D4. The multistep–multistage RK evaluates residuals again at every use.** The stage loop calls `%t` (the residual)
inside the loops over earlier stages and steps (`ms_runge_kutta_ssp.F90:136-154`). For (4,5,8) that is several times
the s + (k − 1) residuals the scheme needs per step. It is correct, but it is the cost model SELAH exists to avoid.
A port caches each residual once.

**D5. Linear SSP stage times are `t + (s−2) dt`** (`runge_kutta_lssp.F90:323`, `:345`). These are not stage abscissae.
That is harmless only because the family is restricted to linear autonomous problems, which FOODIE's suite notes
(`foodie_test_suite_utils.f90:11`).

No defect was found in the leapfrog RAW filter (`leapfrog.F90:128-133`, Williams' filter: d = ν/2 (x_{n−1} − 2x_n +
x_{n+1}); x_n += α d; x_{n+1} −= (1 − α) d). The 2N recurrence (`runge_kutta_low_storage.F90:244-249`) is SELAH's.

## 3. Tests

FOODIE's active tests are the `tests-*` modes of its fobos. `src/tests/quarantena/` is not built.

| FOODIE test | What it verifies | SELAH equivalent | For SELAH |
|---|---|---|---|
| `suite/foodie_test_convergence_order.f90` | every scheme, the Riccati equation U_t = a t U² (**nonlinear, non-autonomous**). dt halved 10 times; error L∞ over all steps (the final-time error can vanish); the order from the finest pair with errors in [1e-12, 1e-3], ±0.3, verified to order 6 | `selah_test_schemes` (y' = −y, y' = −y², both autonomous), `selah_test_stepper` (advection, autonomous) | **gap**: no SELAH test evaluates a time-dependent residual, so `stage_time(k)` has never been checked through an order. Port the Riccati case through the stepper (T1) |
| `suite/foodie_test_polynomial_exactness.f90` | U_t = q t^{q−1}: exact for q ≤ p (quadrature conditions and stage times), not exact for q = p + 1 in an LMM (order not underestimated) | `selah_test_schemes` (order conditions in quad, stronger for RK); the stage-coefficient check of `selah_test_stepper` (q = 1 only) | **gap**: polynomial exactness through the stepper with t-dependent residuals, q ≤ p (T2) |
| `tester/` (`foodie_tester.f90`, integrands `lcce`, `oscillation`, `polynomial`, `riccati`, `ladvection`) | a CLI driver: any scheme on any problem, solution export, order estimate | `scripts/run_vv.sh` with the `selah_vv_*` programs | no: SELAH's V&V programs cover the same ground non-interactively |
| `lcce/foodie_test_lcce.f90` | linear constant-coefficient equation, CLI | `selah_vv_stability` (stability polynomials), `selah_test_schemes` | covered |
| `oscillation/foodie_test_oscillation.f90` | inertial oscillation (nonlinear ODE system), CLI | `selah_vv_splitting_energy`, the harmonic oscillator of `selah_test_stepper` | covered |
| `ladvection` integrand | 1-D linear advection PDE | `selah_test_stepper` (multi-block advection, central flux) | covered |
| `scripts/audit_rk_tableaus.py` | RK order conditions | `selah_test_schemes` (quad precision, in the test suite) | covered (SELAH's runs in CI) |
| `quarantena/regression/{lorenz,burgers,euler-1D}`, `quarantena/parallel/{openmp,caf}` | regressions and parallel drivers, not built | `selah_test_golden`, the device tests | no |

What the FOODIE suite lacks, and SELAH has: stability regions, SSP coefficients and TVD, splitting energy,
conservation and telescoping, bitwise ADAM parity, golden values, device tests, R4P builds, mutation testing.

What neither has, and D1 and D2 show is needed: **variable-dt tests** for any scheme that accepts a variable dt, and
**tests at the ends of a parametrised family** (the largest stage count, not only the default).

## 4. Recommended work, in priority order

For a block-structured CFD and plasma solver on GPUs (ADAM), one line of reasoning each:

1. **T1 + T2, the non-autonomous tests through the stepper** (Riccati, polynomial quadrature). These come first: they
   are cheap and they check `stage_time` on the path every consumer with time-dependent boundary conditions or
   sources uses. Before T1 is written, decide what a splitting's `stage_time` means: it currently returns 0. For a
   non-autonomous splitting, time is a third partition, and that is undefined in SELAH today.
2. **The 2N low-storage schemes 6/7 (HALE-RK) and 12/13/14 (Niegemann–Diehl–Busch).** They need coefficients only, run
   in 2 registers, and have stability regions optimised for advection-dominated problems: the cheapest way to raise
   the stable dt in an explicit CFD code. Re-derive 12/13/14 from the paper (FOODIE's are typed as doubles).
3. **Embedded weights `b_hat` for Heun–Euler, Cash–Karp, Dormand–Prince and Calvo.** PLAN §0 already allows embedded
   pairs as coefficients. The only new code is a kernel variant writing dt Σ (b − b̂) F; step control stays the
   caller's. Feagin 10(8) is not worth 17 stages here.
4. **Leapfrog + RAW**, already deferred from M4: the time integrator of Yee-type plasma solvers (PRISM).
5. **LMM-SSP and LMM-SSP-VSS, with the published variable-step formula:** one residual per step. AMR time refinement
   changes dt, so the variable-step form is the useful one. This needs the history form shared with 4 and 6: a new
   L2 form and a start-up (the caller's, or an RK starter).
6. **MSRK-SSP**, on the same history machinery, with every residual computed once (D4).
7. **Not recommended now:**
   - Adams–Bashforth: no SSP property and no consumer;
   - the linear-SSP RK: only for a linear autonomous consumer (a constant-coefficient Maxwell solve); compute the
     factorial in `R_P` (D2);
   - AM, ABM and BDF: implicit, out of scope for v0.x, behind the `stage_gamma` design of PLAN "Beyond v0.x".

## Reproduce

```bash
python3 notes/foodie-audit-check.py --foodie ~/fortran/FOODIE > notes/foodie-audit-check.txt   # mpmath, ~1 min
```
