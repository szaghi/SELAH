# Plan: one source of truth for time operators (SELAH), one for space operators (GIDEON, SETH)

*Written 2026-11-18. **Status: done (2026-11-19), tracked in issue #2.** The plan was written against SELAH `79c57ba`,
GIDEON `e84ac7b` (v0.3.0) and SETH `0bccbb2`; §9 records what each repository holds now.*

## 0. Goal and the owner's decisions

**Goal.** Time operators live in SELAH only, space operators in GIDEON (FD/FV stencils) and SETH (WENO/TENO) only.
No library's test suite re-implements the other's operators. Each library is verified in its own repository, so an
operator is verified once, not once per library plus once per test suite that copies it.

**Decided by the owner (2026-11-18):**
1. A test failing because the *other* library regressed is acceptable: each library is verified independently.
2. No test-only dependency mechanism in FoBiS. `fobis fetch` clones the direct dependencies only (no recursion), so
   listing a dependency used only by the tests costs a clone and nothing else.
3. The goldens change. New files get new names (`*_v2`); the v1 files are deleted only when nothing reads them.
4. A stepper for 1-D, 2-D, 3-D and 4-D fields is mandatory, besides the 5-D one ADAM uses.

**Decided by the owner (2026-11-18), §6:** D1, the docs examples use GIDEON's operators; D2, `gideon_vv_advection`
is excluded from GIDEON's `tests-gnu-r16`.

## 1. What changes, by repository (the 2026-11-18 audit)

| Repository | Programs | Today | After |
|---|---|---|---|
| **SELAH** | `selah_test_stepper` | 1-D FV advection: central face (u_L+u_R)/2 (`:311-321`), ghost fill (`:287`), boundary flux (`:326`) | `gideon_fv_face_c2`; the ghost fill stays (a boundary condition, the test's own) |
| | `selah_test_stepper_options` | the same residual, copied (`:252-262`) | `gideon_fv_face_c2` |
| | `selah_vv_forms` | central FD −(u_{i+1}−u_{i−1})/2h (`:68`) | `gideon_fd_d1_c2` |
| | `selah_vv_ssp_tvd` | first-order upwind (`:96`) | `gideon_fv_face_ul1` (the face value of the upwind cell) |
| | `selah_test_golden` | 3-point linear stencil 0.25(u_{i−1}+u_{i+1}) − 0.375 u_i (`:105`) | GIDEON's c2 second derivative + a reaction term; golden **v2** |
| | docs examples | central FD in `stepper`, `low_storage`, `production`, `gpu`; central + upwind faces and reflux in `cookbook_amr`; boundary flux in `flux_register`; `masks` | GIDEON, if D1 is decided yes |
| **SETH** | `seth_test_advection:138`, `seth_test_convergence:135`, `seth_vv_scalar:55` (`scalar_step`), `seth_vv_euler:425` (`euler_step`), `seth_vv_euler2d:177`, `seth_vv_adr:182`, `seth_vv_nonpoly:359,496,642` | 10 hand-written SSP-RK3 steps. `seth_vv_cases`, `_euler_order`, `_advection`, `_sod` and `_burgers` reach them through `euler_step` and `scalar_step` | SELAH `SCHEME_SSP_33`, `FORM_SHU_OSHER_2R`, through the 1-D/2-D/3-D stepper |
| **GIDEON** | `gideon_vv_advection:57-61` | a hand-written classic RK4 | SELAH `SCHEME_RK4`, `FORM_BUTCHER`, 1-D stepper |

What stays in each library because it is not the other library's operator:
- **In SELAH:** the scalar ODEs (`selah_test_schemes`), the oscillators (`osc_residual`, splittings, splitting
  energy), the stability polynomials, every kernel and device test (pointwise; their ghost cells are sentinels of the
  box bounds), and the benchmarks (the residual is left in the register).
- **In every test:** boundary conditions and ghost fills. A ghost fill is the test problem's boundary condition, not
  a discretisation operator. GIDEON and SETH do not provide them, and nobody else should.

The field shapes the stepper must take:
- SETH: `u(0:n-1)`, rank 1, periodic by index arithmetic, no ghosts; `q(3, 1-NG:n+NG)`, rank 2, components first;
  `q(4, 1-NG:n+NG, 1-NG:n+NG)`, rank 3 (`seth_vv_euler2d:104`).
- GIDEON: `u(N)`, rank 1.

## 2. Phase 1: SELAH, the stepper on fields of rank 1 to 4 (prerequisite of phases 2 and 3)

**Design (recommended): rank-k fronts over the rank-5 core, by pointer bounds remapping.**
- **One kernel source, as today.** The rank-5 box kernels and the stepper stay the only implementation.
- **The front.** A rank-k field (k = 1..4) is remapped to rank 5 with singleton dimensions, by F2008 bounds
  remapping of a simply contiguous target. That means a `contiguous` dummy or a whole array: no copy.
- **The mapping of dimensions to (v, i, j, k, b):**

  | Field | As rank 5 |
  |---|---|
  | rank 1, scalar `(i)` | `(1, i, 1, 1, 1)` |
  | rank 2 `(v, i)` | `(v, i, 1, 1, 1)` |
  | rank 3 `(v, i, j)` | `(v, i, j, 1, 1)` |
  | rank 4 `(v, i, j, k)` | `(v, i, j, k, 1)` |
  | scalar 2-D/3-D `(i, j[, k])` | `(1, i, j[, k], 1)` |

  The rank-k box constructors say which case applies (an explicit `components` argument), never a guess.
  `LAYOUT_V_LAST` stays rank 5 only (ADAM's FNL pool), until a caller needs it.
- **Interfaces:**
  - `begin_step`, `begin_stage`, `end_stage` and `finish_step` become generic over the rank of `q`, `dq` and `mask`;
  - `initialize` takes `lb`/`ub` of size k;
  - the box gains rank-k constructors (`box_interior_1d` … `_4d`, or one generic `box_interior` resolved by its
    arguments).
- **Pointer queries.** A function cannot be generic on the rank of its result, so `stage_state()` and
  `rhs_register()` keep their rank-5 meaning. Rank-k subroutine forms are added, generic on the rank of their pointer
  argument: `call st%stage_state_view(u)` and `call st%rhs_register_view(r)`. Both are null outside a stage window,
  as now.
- **The registers stay rank 6 internally,** with singleton dimensions.
- **Why not assumed rank (`dimension(..)`, `select rank`):** nvfortran's support inside OpenACC regions is the weak
  point, and it would put a `select rank` in front of every kernel launch.
- **Why not per-rank kernels:** 7 operations × masked/unmasked × 4 ranks of loop nests. That is a second source of
  truth inside SELAH itself.
- **The risk:** a remapped pointer to device data must resolve in OpenACC's present table, which looks data up by
  address, so it should. This is to be proved on the device before anything else (first gate below). If it fails for
  some kernel, that kernel alone gets a rank-specific loop nest from the same `.INC` body.

**Steps and gates:**
1. **A spike: the remap on the device.** One rank-1 field through `begin_stage`/`end_stage` with `tests-nvf-oac` on
   the RTX 4070 and the A30.
   - *Gate:* bitwise equal to the same field declared rank 5.
2. **The fronts, box constructors and view subroutines;** the stage contract restated for rank k.
   - *Gate:* the whole suite on gnu, gnu-r4, nvf, intel, nvf-oac.
3. **Tests per rank.**
   - Observed orders through the stepper on rank 1, 2, 3 and 4 fields, every scheme and form, with a problem whose
     semi-discrete solution is closed-form. Per-rank `selah_test_stepper` cases first use the current central
     residual; phase 4 replaces it with GIDEON's.
   - The contract checks (ghost sentinels, stage window, stage coefficients) per rank.
   - The result is bitwise equal to the rank-5 run of the same data.
   - *Gate:* as step 2, plus the mutation run.
4. **Docs.** Guide (stepper, box), the API reference, a tutorial section on a 1-D field.
   - *Gate:* docs build, `check_docs_links.py`, `check_api_doc.py`.

## 3. Phase 2: SETH's tests take their time integrator from SELAH

1. **fobos.**
   - `[dependencies]` gains `FUNDAL :: branch=main` and `SELAH :: branch=main`.
   - `$EXDIRS` gains `FUNDAL/src/tests/ FUNDAL/docs/ SELAH/src/tests/ SELAH/src/bench/ SELAH/docs/ SELAH/notes/`.
   - The test templates get `include = src/third_party/FUNDAL/src/lib/` and `-DFUNDAL_NO_MPI`.
   - The library modes are unchanged: SETH's library does not use SELAH.
2. **fpm.** `[dev-dependencies] SELAH = { git = "https://github.com/szaghi/SELAH", branch = "main" }`. In fpm this is
   the native form; FUNDAL comes through SELAH's own manifest.
3. **Each program in §1:** the hand-written SSP-RK3 is replaced by the stepper (`SCHEME_SSP_33`,
   `FORM_SHU_OSHER_2R`, rank 1, 2 or 3). Where a stage checks positivity (`seth_vv_euler`, `_euler2d`: `ok = .false.`
   when ρ ≤ 0 or p ≤ 0), the check moves into the stage window, on `stage_state` or q, unchanged in meaning.
   `euler_step` and `scalar_step` keep their signatures, so their callers do not change.
4. **The effect on results.** The arithmetic order changes:
   - SETH computes 0.75 u + 0.25 (u1 + dt r);
   - SELAH's 2-register stage computes a u⁰ + b q + dt c r, with a = 3/4, b = c = 1/4.

   Every SETH result moves by round-off, larger in the long runs. The SETH goldens that contain a time integration
   get v2 files (the others do not change: they are spatial kernels alone), and the V&V report data and figures are
   regenerated.
5. **Gate.**
   - SETH's whole suite in every mode (gnu, gnu-debug, gnu-r4, nvf, intel, device where SETH has one).
   - Every V&V verdict unchanged: observed orders within their tolerances, positivity, TV bounds.
   - The report-data diff reviewed: round-off, except where a published figure changes visibly.
   - SETH's Matrix CI green on main.

## 4. Phase 3: GIDEON's tests take their time integrator from SELAH

1. **fobos and fpm** as for SETH (phase 2, steps 1–2).
2. **`gideon_vv_advection`:** the RK4 loop becomes the 1-D stepper (`SCHEME_RK4`, `FORM_BUTCHER`). The energy check
   (Σu² never grows at CFL 0.25) moves into the step loop unchanged.
3. **Quad precision (D2).** `tests-gnu-r16` builds every GIDEON test in R16P, and SELAH refuses R16P
   (`selah_kind_parameters`: `#error`).
4. **Gate:** GIDEON's suite in every mode; the advection verdicts unchanged; GIDEON's CI green.

## 5. Phase 4: SELAH's tests take their space operators from GIDEON

1. **fobos.** `GIDEON :: branch=main` in `[dependencies]`; `$EXDIRS` gains `GIDEON/src/tests/ GIDEON/docs/
   GIDEON/scripts/`. This matters: GIDEON's tests will `use selah`, so they must never be compiled inside SELAH's
   build. fpm: `[dev-dependencies] GIDEON`.
2. **The five programs of §1**, each residual written as a device loop over the box calling GIDEON's pure kernels.
   GIDEON's kernels are `!$acc routine seq`, so the device builds keep working.
3. **The time error stays isolated.** Each replacement operator is the same discretisation the test uses today, so
   the closed-form semi-discrete solutions (central: wave speed sin(2πh)/(2πh)) and the TVD argument (first-order
   upwind) still hold. The test documentation cites the GIDEON kernel by name.
4. **Golden v2:** `selah_golden_v2.txt` and `selah_golden_r4_v2.txt`, written by the golden program's write mode,
   with `selah_test_golden` reading v2. v1 is deleted in the same commit, since nothing reads it any more.
5. **The verification pages** (A3 stepper orders, P1 SSP-TVD, P4 conservation and telescoping, C3, C4 golden): the
   setup tables name the GIDEON operators; the report data and figures are regenerated.
6. **Gate.**
   - SELAH's whole suite on gnu, gnu-r4, gnu-omp, nvf, intel and nvf-oac (RTX 4070), and an A30 job.
   - The mutation run: 28 or more killed, the stepper tests' sensitivity unchanged.
   - Docs build and checks; CI green.

## 6. Phase 5: SELAH's docs examples (D1)

**D1, decided 2026-11-18: yes.** Should the tutorial and cookbook examples (`docs/examples/src/`) use GIDEON's operators?
- **Yes (recommended).** The tutorial then shows the ecosystem as it is meant to be used: the caller's residual
  (principle P2) written with GIDEON, the time step with SELAH. A reader who copies an example copies the real
  pattern. The cost: the docs-examples workflow fetches GIDEON (already fetched for the tests after phase 4), and
  chapter 1 ("a step by hand") gains a GIDEON call next to the hand-written step.
- **No.** The examples keep minimal inline residuals, labelled as the caller's code. That is two sources for a
  central difference in SELAH's repository, the situation this plan removes.

**D2, decided 2026-11-18: (a).** GIDEON's R16P test mode and SELAH:
- **(a), recommended.** Exclude `gideon_vv_advection` from `tests-gnu-r16`: the mode verifies the orders of the
  operators, which needs no time integration.
- **(b).** Allow R16P in SELAH's host builds. The kernels would compile at R16P (no device), at the cost of a mode
  nobody else needs.

## 7. Phase 6: keeping it separated

A check in each repository's test suite, run in CI, written once (`scripts/check_separation.py`, the same file in
the three repositories):
- **SELAH:** no stencil in `src/tests/` (an index of the form `i±1`, `j±1` or `k±1` inside an assignment to a
  residual), except in ghost fills, which are named `fill_ghosts` and checked by name.
- **SETH and GIDEON:** no Runge–Kutta combination in `src/tests/` (`0.75…*u + 0.25…`, `k1 … k4`), and no stage
  loop not driven by `stepper`.

The regular expressions are allowed false positives, with a named waiver list in the script. The goal is to stop a
copy coming back unnoticed, not to parse Fortran.

## 8. Order, dependencies between phases, risks

| Phase | Needs | Repositories |
|---|---|---|
| 1, N-D stepper | — | SELAH |
| 2, SETH ← SELAH | 1 | SETH |
| 3, GIDEON ← SELAH | 1 (and D2) | GIDEON |
| 4, SELAH ← GIDEON | — (can run in parallel with 2–3; sequenced after 1 to change one thing at a time) | SELAH |
| 5, docs examples | 4, D1 | SELAH |
| 6, guard | 2–5 | all three |

**Risks:**
- **Cross-repository breakage.** With `branch=main` everywhere, an API change in SELAH breaks SETH's and GIDEON's CI
  at once. Rule: an API change to SELAH lands together with the commits adapting SETH and GIDEON, checked with
  `fobis ecosystem` (a cross-repository build and test) before any of them is pushed.
- **The dependency cycle at test level.** SELAH's tests use GIDEON, and GIDEON's tests use SELAH. It is safe because
  `fobis fetch` is not recursive and each repository's `$EXDIRS` excludes the other's tests. One check, in phase 4's
  gate: a clean clone of each repository builds its test modes from scratch.
- **The device remap** (phase 1, step 1). If it fails on nvfortran, the fallback is per-kernel rank-specific loop
  nests from the same `.INC` bodies, so the cost is bounded.
- **The SETH V&V figures** regenerate with round-off changes. A visibly changed published figure gets its page
  updated in the same commit.

**After this plan:** M7, the ADAM wiring, on a SELAH whose own tests use GIDEON's operators: the same pairing ADAM
will run in production.

## 9. Done: the final state of each repository (2026-11-19)

| Phase | Repository | Commits | What holds |
|---|---|---|---|
| 1 | SELAH | `a0dc4a0` `fa6dd51` `1ea5403` | the stepper on fields of rank 1 to 4 (one rank-5 implementation, explicit-shape bridges, `stage_state_view`/`rhs_register_view`); 848 checks, A30 jobs 1368451, 1368457 |
| 2 | SETH | `11ebd46` `dbaae69` `2328cb0` `703e66d` `a76bd80` | the 10 hand-written SSP-RK3 steps are SELAH's `SCHEME_SSP_33` (rank 1, 2, 3); the report utilities split out so that only the 16 time-stepping programs depend on SELAH; those programs leave the OpenACC modes |
| 3 | GIDEON | `01232a0` `204e075` | `gideon_vv_advection` advances with SELAH's RK4 (rank 1); out of `tests-gnu-r16` (D2) and of the OpenACC modes |
| 4 | SELAH | `16d006b` `1e85348` `78b40f9` `bdd84c6` | the tests' space operators are GIDEON's (`fv_face_c2`, `fd_d1_c2`, `fv_face_ul1`, `fd_d2_c2`); golden v2 (12 R8P values 1 ulp, 21 R4P) |
| 5 | SELAH | `5e75c51` | the docs examples' residuals are GIDEON's (D1); `static-gideon-gnu`, `docs_examples.sh` links `libgideon.a`; no example output changed |
| 6 | SELAH, SETH, GIDEON | (this commit's night) | `scripts/check_separation.py`, the same file in the three repositories, run by every Matrix job: 0 findings in each, planted copies caught |

What did not go as planned, and why:
- **The OpenACC test modes** of GIDEON and SETH exclude the programs that use SELAH: SELAH and FUNDAL do not compile
  under `-acc` without SELAH's device build (`fundal.H`'s host clause), and those programs hold host fields. Running
  them on the device would need them rewritten for FUNDAL memory: not planned, not needed for the separation.
- **Golden v2** came from 1-ulp changes, inside the tolerance: the C4 rule asks for a new file for any deliberate
  change, so v2 it is.
- **SETH's report data** (tables, figures) of the 16 programs were not regenerated: round-off changes only; the 2-D
  cases' report mode takes hours.
- **SETH's docs examples** (2026-11-19, after the plan): t06–t10 and `cookbook_euler` take their step from SELAH
  (`static-selah-gnu`, linked by SETH's `docs_examples.sh`); t03 keeps its hand-written step (chapter 3 teaches what a
  step is) and t11_gpu its own until it is wired to SELAH's device build: both are named waivers of
  `check_separation.py`, which now scans GIDEON's and SETH's `docs/examples/src/` too.
- **GIDEON's `cookbook_dissipation`** example writes its SSP-RK3 by hand, a copy the plan's audit had not listed (found
  when the check's scan reached the docs examples): a named waiver, pending the same migration (the owner's decision).

**Changelog.** No hand-written entries: each repository's `CHANGELOG.md` is generated from the commit messages by
git-cliff at the next release (`scripts/release.sh`), and the commits above are written for it.

**Next: M7, the ADAM wiring** (`notes/PLAN.md` §4), is unblocked: SELAH's own tests now pair SELAH's time operators
with GIDEON's space operators, as ADAM will in production.
