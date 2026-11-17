# SELAH: implementation plan

*Name:* **SELAH**, *Stage-wise Evolution Library for Accelerators and HPC* (decided 2026-11-17). Repository `szaghi/SELAH`, module prefix `selah_`, library files `libselah.a` / `libselah_r4.a`; **SELAH** in prose, *Selah* (italic) for the biblical word. *Selah* is the mark that appears 71 times in the Psalms and 3 in Habakkuk, usually read as a pause where the music stops and the singer breathes (the exact meaning is debated): the stage boundaries, where the orchestrator acts before the integrator continues. This file is the design record of the project (moved into the repository at M0).
*Origin:* the ADAM time-integrator audit of 2026-11-17 and ADAM issue [#128](https://github.com/szaghi/adam/issues/128) (its defects are fixed in ADAM independently of this plan).
*Status:* draft, nothing implemented.

---

## 0. Goal, non-goals, consumers

**Goal.** Time-evolution operators for structured, block-based field solvers on CPU and GPU. v0.x is explicit only; the design must not exclude implicit stages (see "Beyond v0.x" below). The operators must:
1. are at least as fast as ADAM's current RK layer, kernel for kernel, and reduce the integrator's memory traffic where the stage storage allows it;
2. expose a **stage contract** that an AMR / multi-realm orchestrator (ADAM's forest) can drive: per-stage begin/end, a documented stage-state buffer, effective stage weights for reflux;
3. carry a verified scheme catalogue: order conditions in quad precision, observed orders, SSP coefficients.

**Non-goals for v0.x:**
- implicit or IMEX stages, dual time stepping, Newton/Krylov;
- adaptive time-step *control*: embedded pairs ship as coefficients only;
- a polymorphic ODE abstraction (FOODIE's `integrand` with overloaded operators): it is the performance anti-pattern this library exists to avoid;
- particle pushers (Boris and similar). The particle *stage arithmetic* is in scope, through the flat kernels; the push physics is not.

**Beyond v0.x: implicit stages fit the contract, they do not change it.** A diagonally implicit stage (DIRK/SDIRK) is `U_s = q + dt Σ_{j<s} a_sj F_j + γ F(U_s)`, with `γ = dt a_ss`:
- the library builds the explicit part with the same `stage_state` kernel and exposes `stage_gamma(k)`;
- the caller solves the stage equation (Newton or JFNK) in its own code, as it already owns the residual (P2);
- IMEX (ARK) schemes are two tableaux over the same stage loop; multistep methods (BDF, Adams) are register rotations plus the same combinations.

Fully implicit RK, with stage solves coupled across stages, is the one real extension. To keep these doors open from M1:
- the L0 types do not assume `A` strictly lower triangular, except where the scheme declares it;
- the stage API returns, besides the state, the coefficient of the implicit term (0 for explicit stages).

**Consumers.** The first consumer is ADAM (FLUME, PRISM, RIPPLE). The second consumer is **undefined**, and that is the main strategic risk (§7). Xall_GPU would need dual time stepping, which is out of scope for v0.x.

---

## 1. Design principles (each one checkable)

| # | Principle | How it is checked |
|---|-----------|-------------------|
| P1 | No polymorphism, procedure pointers or callbacks inside or around a kernel. Scheme dispatch is a host `select case` once per kernel launch. | Code review; `grep` for `class(` and `procedure(` in the kernel modules returns nothing. |
| P2 | The caller owns the time loop and the residual evaluation; the library never calls back. | API (§3, L2). |
| P3 | Kernels act on plain contiguous arrays passed to dummies with explicit lower bounds, not assumed-shape rebased to 1. | Ghost-sentinel tests (§4, M2). |
| P4 | Box kernels update the **interior only**; the integrator never writes a ghost cell. | Ghost-sentinel tests: NaN sentinels in the ghosts survive bitwise. |
| P5 | One working precision per build (PENF `R_P`, `-D_R_P_IS_R4P` / R8P), as in SETH. | `r4` modes compile and pass. |
| P6 | One kernel source for host OpenMP, OpenACC and OpenMP offload, through FUNDAL's `fundal.H` macros (`DEVICEVAR`, `OMPLOOP`, `DEVICEPTR`), as ADAM's `adam_fnl_rk_kernels` already does. | Each backend mode builds from the same `.INC` bodies. |
| P7 | The floating-point summation order of every kernel is specified and documented, so golden tests can be bitwise. | Bitwise parity tests against ADAM's kernels (§4, M2/M3). |

---

## 2. Repository and infrastructure

- **Repository:** `~/fortran/SELAH`, GitHub `szaghi/SELAH`, private until v0.1.0 (SETH precedent).
- **Build:** FoBiS, `fobos` modelled on SETH's (`[project]`, `[scaffold]` matrix linux-gnu / linux-nvhpc / linux-oneapi, `[ecosystem] check`, `[dependencies]` through `fobis fetch`). Modes:
  - library: `static-{gnu,nvf,intel}`, `shared-{gnu,nvf,intel}`, the same with `-r4`;
  - host backend: `tests-gnu`, `tests-gnu-debug`, `tests-gnu-coverage`, `tests-gnu-omp` (host OpenMP: `OMPLOOP = parallel do`);
  - device backends: `tests-nvf-oac` (OpenACC, `-gpu=$GPU_CC`), `tests-gnu-omp-offload` (compile only until an AMD/NVIDIA offload runner exists);
  - benchmarks: `bench-gnu`, `bench-nvf-oac`;
  - varsets: `cc80` (A30, medooza), `cc89` (WSL RTX), `cc90`.
- **Dependencies:** PENF (kinds) from M0; FUNDAL (the `fundal.H` macros; `dev_alloc`/`dev_free` only in L2, when the library owns the registers) joins at M2, with the first kernels. Nothing else: no FiNeR, no MPI. The library is per rank, and MPI belongs to the caller.
- **Toolchains:** gfortran 16.2.0 (+ Open MPI 5.0.11 for the ADAM-wiring tests only) and nvfortran 26.9.
  - **CI** (GitHub Actions, fortran-toolchains images): CPU modes only.
  - **Device tests and benches:** run by hand on the WSL box (cc89) and on medooza (A30, cc80, `dev-gn`), with `nvhpc-job.slurm`-style templates in `scripts/`. Every result recorded in the docs states the job, the node and the commit.
- **Docs:** VitePress, mirroring StringiFor/PENF and SETH's `guide/` + `verification/` split. ⚠ VitePress fleet reminder: pin the vitepress/mermaid stack in a committed `package-lock.json` (mermaid ≥ 11.16 needs the es2022 target fix already in harp); do not start on vitepress 2 alpha.
- **Licence and version numbers:** your decision; not proposed here.

---

## 3. Architecture: L0–L3

### L0: coefficients (`selah_tableaux`, pure, no dependencies beyond PENF)

Types (plain data, no type-bound behaviour beyond pure queries):

```fortran
type :: butcher_tableau            ! RK: A (strictly lower triangular for explicit schemes, a flag says so), b, c = A·1
   integer(I4P)           :: stages, order
   real(R_P), allocatable :: a(:,:), b(:), c(:), b_hat(:)   ! b_hat: embedded weights, optional
   real(R_P)              :: ssp_coefficient = 0           ! 0 when the scheme is not SSP
endtype
type :: shu_osher_form             ! u(i) = sum_k alpha(i,k) u(k) + dt beta(i,k) F(u(k))
type :: low_storage_2n             ! Williamson: dQ = A_s dQ + dt F ; q = q + B_s dQ
type :: low_storage_2s             ! Ketcheson 2S / 2S* / 3S* (gamma, beta, delta)
type :: splitting_composition      ! sequence of (partition, coefficient) flows, palindrome expanded at init
```

Every scheme knows its admissible storage forms and its **stored-field count** per form. ADAM's #128 D1 bug came from confusing that count with the number of stages; here it is a first-class query, `registers_needed(scheme, form)`.

Catalogue for v0.1, each with its reference and its digits' provenance:

| Scheme | Forms | Fields | Source |
|---|---|---|---|
| forward Euler | Butcher | 1 | — |
| SSP(2,2), SSP(3,3) | Butcher, Shu–Osher (2 registers) | 2/3 or 2 | Shu & Osher 1988 |
| SSP(5,4) | Butcher (5 + U) | 6 | Spiteri & Ruuth 2002; 15-digit Shu–Osher form of Isherwood, Grant & Gottlieb 2018, Butcher derived (as ADAM #104) |
| SSP(10,4) | low-storage, 2 registers | 2 | Ketcheson 2008, SISC 30: C = 6, C_eff = 0.6 (vs 0.30 for SSP(5,4)) |
| classic RK4 | Butcher | 4 + U | — |
| LSRK(5,4) | 2N | 2 | Carpenter & Kennedy 1994 |
| LSRK3 | 2N | 2 | Williamson 1980 |
| Yoshida 4 | splitting | 0 | Yoshida 1990 |
| Blanes–Moan S6 (order 4) | splitting, palindrome | 0 | Blanes & Moan 2002 (the correct composition; see #128 D3) |
| leapfrog + RAW filter | multistep | 2 | Williams 2011 |

CFM is excluded until a matrix-exponential action exists (#128 D2).

**L0 tests** (pure, no arrays; `test_rk_ssp_tableau` is the template):
- order conditions up to the declared order, all rooted trees, in quad precision on the stored doubles (≤ 1e-14);
- `c = A·1`, strict lower triangularity, `Σ b = 1`;
- SSP coefficient recomputed from the Shu–Osher form (non-negative coefficients, the min ratio) against the declared value;
- equivalence of forms: the Butcher tableau derived from each Shu–Osher / 2N / 2S form equals the stored one (≤ 1e-15);
- observed order on `y' = -y` and on `y' = -y²` (±0.1); splittings on the harmonic oscillator. This is the replica test that exposed #128 D2/D3.

### L1: kernels over arrays the caller owns (`selah_kernels_*`)

The data model is the **box**, a small integer descriptor passed **per call** (never stored at init, because ADAM's active-block count changes at each regrid):

```fortran
type :: box
   integer(I4P) :: lb(5), ub(5)   ! allocation bounds of the rank-5 field (ghosts included)
   integer(I4P) :: lo(5), hi(5)   ! loop range: interior cells, components, active blocks 1..blocks_number
endtype
```

The arithmetic is pointwise, so it is layout-agnostic as long as every operand shares the layout. The **mask** is the only operation that must know which dimension holds the components. Hence two layout variants:
- `V_FIRST`: `(v,i,j,k,b)`, ADAM CPU;
- `V_LAST`: `(b,i,j,k,v)`, ADAM FNL.

The mask is a rank-5 array whose component dimension is selected by `mask_comp`; a cell is active where `mask(mask_comp,…) < 0`. That is exactly ADAM's `phi(all_solids,…)`, so ADAM needs no copy and no strided slice that would force copy-in temporaries.

Kernel families. Each is generated from one `.INC` body and instantiated for {`V_FIRST`, `V_LAST`} × {unmasked, masked}, plus a **flat** rank-1 instance for contiguous whole arrays (particles, PML slabs: PRISM's `rk_pic`/`rk_pml`/`rk_bc` clones):

| Kernel | Operation (per point) | Used by |
|---|---|---|
| `stage_state` | `U = q + dt Σ_{j<s} a_sj F_j` (private sum, then one store) | Butcher forms |
| `step_update` | `d = Σ_s b_s F_s ; q = q + dt d ; [dq_out = d]` (fused) | Butcher forms (fixes #128 D4 by construction) |
| `shu_osher_stage` | `q = α q0 + β q + γ dt R` | 2-register SSP |
| `lsrk_2n_stage` | `dQ = A_s dQ + dt R ; q = q + B_s dQ` (fused) | Williamson / Carpenter–Kennedy |
| `lsrk_2s_stage` | Ketcheson 2S/2S*/3S* update | SSP(10,4) etc. |
| `partition_update` | `q(v∈P) = q(v∈P) + c dt R(v∈P)` | splittings; `P` is a component range |
| `leapfrog_raw` | leapfrog step + RAW filter | leapfrog |
| `copy_interior` | `U = q` over the box | stage initialisation when needed |

Each kernel's summation order is written in its header comment (P7). The ADAM-compatible order is FNL's: private sum, then one update.

### L2: the stepper the caller drives (`selah_stepper`)

The integrator never owns the loop. It answers "what do I compute now, where do I write it":

```fortran
call stepper%initialize(scheme='ssp-54', form=FORM_BUTCHER, bx=bx, storage=STORAGE_OWNED)  ! or STORAGE_EXTERNAL
call stepper%begin_step(q)
do k = 1, stepper%stages()
   call stepper%begin_stage(k, dt, q)                 ! builds the stage state (L1 stage_state)
   call app_residual(stepper%stage_state(), stepper%rhs_register(k))   ! caller writes F_k in place: no assign copy
   call stepper%end_stage(k)
enddo
call stepper%finish_step(q, dt, dq_out=dq)            ! fused update + step residual
```

Queries:
- `stages()`, `stage_time(k)` (c_k), `stage_weight(k)`: the effective Butcher b_k for any form, which is what the reflux register needs;
- `is_ssp()`, `ssp_coefficient()`;
- `registers_needed(scheme, form)` (pure, for memory budgets before allocation);
- `state_aliases_q()`: true for 2N/2S forms, where the stage state *is* `q`, which the orchestrator must know when it fills seam ghosts.

**The stage contract** (the document ADAM's forest implements against):
1. between `begin_stage(k)` and `end_stage(k)` the buffer returned by `stage_state()` holds U⁽ᵏ⁾ in the interior;
2. the library never writes ghost cells of any buffer;
3. the caller may write ghost cells of `stage_state()` (seam fill, physical BCs) during that window;
4. `rhs_register(k)` is written only by the caller, and only in the interior;
5. after `finish_step`, all registers are scratch, which makes reuse explicit (FLUME's D_flux uses idle registers today).

For splitting schemes `begin_stage` also returns the **partition** to evaluate, so the caller can compute half residuals. That removes PRISM Yoshida's discarded halves (#128 D7.3).

Storage: `STORAGE_OWNED` allocates registers with FUNDAL `dev_alloc` (device builds) or `allocate` (host builds); `STORAGE_EXTERNAL` takes caller arrays (ADAM's fill-the-node pool). The backend is one per build, as in ADAM and SETH: no runtime host/device switch.

### L3: adapters (in the consumer, not in the library)

ADAM's `rk_object` / `rk_fnl_object` become thin adapters: INI parsing (`[runge_kutta]`), `mpih` messages, `grid`/`field` → `box`, `phi` → mask. They live in ADAM. The library ships one **example adapter** (a minimal multi-block advection driver, §4 M4) as the documented integration pattern.

---

## 4. Milestones and gates

| M | Content | Exit gate |
|---|---|---|
| **M0** | Repo scaffold (`fobis scaffold`), fobos, CI (CPU), docs skeleton, licence | `fobis build` default mode green in CI. *Locally green 2026-11-17:* gfortran 16.2 (static, shared, tests, debug, R4P), nvfortran 26.9 (static, tests, debug, R4P; `tests-nvf-oac` built, not run), ifx 2026.2 (static, tests); `selah_test_kinds` passes on each; the docs build (mermaid 11.17.2). **Closed 2026-11-17** at `aba09bf`: CI, fpm, Matrix (linux-gnu, linux-nvhpc, linux-oneapi) and Docs (build + Pages deploy) green. The Matrix needed `fobis fetch` in `scripts/ci/matrix.sh` (FoBiS#219). |
| **M1** | L0 catalogue + tests | all order conditions ≤ 1e-14; observed orders ±0.1; Blanes–Moan S6 at order 4; SSP(5,4) tableau bitwise equal to ADAM's `ssp_tableau` *Done in night 2026-11-17-m1-tableaux (not yet committed):* the catalogue as planned plus Strang; residuals ≤ 5.1e-16 (ssp-54), Blanes–Moan S6/S10 at 4.000/6.00, SSP(5,4) bitwise ADAM; leapfrog + RAW moved to M4, Ketcheson 2S types not needed (`.night/` log). |
| **M2** | L1 host kernels (box V_FIRST/V_LAST, flat, masked/unmasked) | (a) **bitwise parity** with ADAM's kernels, copied verbatim into `src/bench/` as SETH did with `adam_fnl_weno_kernels.F90`: vs the FNL kernels bitwise, vs the CPU `rk_object` loops within the documented round-off envelope (#128 D6); (b) ghost sentinels intact; (c) active-prefix growth: `blocks_number` grows between calls and nothing past it is touched |
| **M3** | L1 device: OpenACC (nvfortran 26.9); OpenMP offload compile-only | (a) bitwise parity with ADAM FNL kernels on device; (b) bench on A30 and cc89 with `!$acc wait` + `system_clock` (never `cpu_time`): each kernel within 2% of the ADAM original for the same operation, and effective bandwidth ≥ 80% of a device STREAM triad |
| **M4** | L2 stepper (Butcher, Shu–Osher, 2N, 2S, splitting, leapfrog), owned and external storage; the multi-block example adapter | (a) observed orders through the stepper on a periodic multi-block advection (WENO-free central scheme: the time error dominates); (b) contract tests: ghost writes by the caller between begin/end survive; `rhs_register` aliasing; `state_aliases_q`; (c) **reflux telescoping test**: two blocks, a fine/coarse face, fluxes accumulated with `stage_weight(k)`; the conservation defect equals the flux-register correction to round-off |
| **M5** | Fused paths measured | integrator memory traffic per step for SSP(5,4): predicted **49 → 28** field words (initialisation copy and assign copy removed, update + step residual fused). The gate is a measured time ratio consistent with the traffic ratio (bench, A30); SSP(10,4) 2-register vs SSP(5,4) at equal simulated time |
| **M6** | Docs (guide, API, verification pages with every measurement and its job/commit); `v0.1.0` (tag yours) | docs build; every catalogue entry has its verification page |
| **M7** | ADAM wiring, in ADAM, own issue: adapter; FLUME first (CPU, then FNL), then PRISM, RIPPLE; delete `adam_fnl_rk_kernels` and the PRISM clones; `seam_ghosts_to_stages` reduced to the single stage-state buffer | FLUME/PRISM regressions bitwise on FNL; CPU within the D6 envelope; FLUME `profile-couette3d` loop time not worse |

M1–M2 are independent of every open decision except the name. M3 needs GPU time: ADAM GPU tests have priority, so schedule the benches around them.

---

## 5. Testing strategy (summary)

- **Unit, L0:** pure, quad-precision order conditions, form equivalence, scalar observed orders.
- **Unit, L1:** each kernel on small hand-checkable boxes; sentinel ghosts; masked/unmasked; both layouts; R4P and R8P.
- **Parity:** verbatim copies of ADAM's kernels as the oracle; bitwise where the summation order matches, an envelope where it does not, never "close enough" by eye.
- **Contract, L2:** the five clauses of the stage contract, each an executable test.
- **Verification:** observed order through the full stepper on a multi-block problem; reflux telescoping.
- **Performance:** bench harness with an explicit device sync, the node and driver recorded, three repeats, the spread reported.

The acceptance tests are written before the kernels (TDD); an agent's claim of "verified" counts only after a rerun.

---

## 6. Decisions you need to take (in this order)

1. ~~Name~~: **SELAH**, *Stage-wise Evolution Library for Accelerators and HPC* (2026-11-17).
2. **Relation to FOODIE:** a new repository, or a FOODIE 2.0 major release that replaces the `integrand` API with L0–L2. Argument for new: FOODIE's public API and its users stay untouched. Argument for FOODIE 2.0: one integrator library in your ecosystem instead of two with overlapping catalogues.
3. **Second consumer:** name one, or accept that v0.x is shaped by ADAM alone and keep the API small until a second consumer exists.
4. **Licence.**
5. Whether OpenMP offload is a tested backend in v0.x, or compile-only until the AMD path of ADAM (#117) is settled.

---

## 7. Risks

| Risk | Effect | Mitigation |
|---|---|---|
| No second consumer | the API overfits ADAM, or generalises speculatively | keep L2 minimal; make each generalisation wait for a concrete user |
| Bitwise parity across summation orders | ADAM CPU goldens shift at round-off on wiring | adopt FNL order (#128 D6) before M7; envelope tests documented |
| Device LTO | a prebuilt `.a` loses ADAM's `-gpu=lto` cross-module inlining | consumers build the library from source through `fobis fetch` (as ADAM does with its dependencies); document `-gpu=lto` for prebuilt builds |
| OpenACC + assumed-shape with explicit lower bounds | descriptor or copy surprises in nvfortran | keep the dummy style ADAM already runs in production; parity tests on device catch regressions |
| Stage contract incomplete for the forest | M7 needs API changes | write the contract tests against ADAM's actual forest usage (#113 seam fill, reflux `b_s`, D_flux scratch) during M4, before wiring |
| Scope creep (implicit, IMEX, adaptivity) | v0.1 never ships | non-goals of §0; separate issues |
