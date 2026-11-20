---
title: Stepper
icon: stepper
---

# Stepper

SELAH's L2 layer is the stepper (`type(stepper)`): one time step that the caller drives one stage at a time. The stepper
owns the registers and the coefficients. The caller owns the solution, the time loop, the ghost cells and the residual.
The stepper never calls back.

## A step

```fortran
use selah
type(stepper) :: st

call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
! ... every step:
call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   ! the stage state: q itself when st%stage_state_is_q(k), else st%stage_state()
   ! fill its ghost cells (exchange, seams, boundaries), then write its residual into st%rhs_register()
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)        ! the Butcher form can also return the step residual: dq=...
! ... at the end:
call st%destroy
```

Every call takes the box (see [Concepts](./concepts#the-box)), and all four step calls take an optional `mask` and
`mask_comp`. Pass the same mask to all four. A splitting also needs its two partitions (`partitions=`), and for each
flow `stage_partition(k)` says which partition to compute the residual of. The rules the two sides keep are the
[stage contract](./concepts#the-stage-contract).

## Procedures

| Procedure | What it does | Preconditions and notes |
|---|---|---|
| `initialize(scheme, form, lb, ub [, partitions] [, registers] [, stages_on_q] [, components] [, masked_q])` | sets the scheme and form for a field of rank `size(lb)` (1 to 5); allocates `registers_needed(scheme, form)` registers with the bounds `lb:ub` through FUNDAL (device memory in a device build), or takes the caller's | stops if the scheme has no such form, or a splitting comes without `partitions`. `partitions(:,p)` = `(v1, v2)`, the components of partition `p`. `components` is required below rank 5 ([fields of rank 1 to 4](#fields-of-rank-1-to-4)). `registers` (see below), for rank 5 only. `masked_q`: the [mask policy](#masks). Calls `destroy` first, so a stepper can be re-initialised |
| `begin_step(bx, q [, mask, mask_comp])` | opens a step: the two-register and Ketcheson forms keep a copy of $q^n$ | |
| `begin_stage(k, bx, dt, q [, mask, mask_comp])` | opens stage `k`: the Butcher form builds its state into its register (unless `stage_state_is_q(k)`) | `k` from 1 to `stages()`, in order |
| `stage_state()` | pointer to the stage state, with the field's bounds | Butcher form only; **null** on a stage evaluated on `q`, and always null in the other forms (their state is `q`). **Null** outside a stage window (before `begin_stage`, after `end_stage`) |
| `rhs_register()` | pointer to the register the caller writes the residual of the active stage into, with the field's bounds | write the cells of the box; for a splitting, the components of `stage_partition(k)`. **Null** outside a stage window (before `begin_stage`, after `end_stage`) |
| `stage_state_view(u)`, `rhs_register_view(r)` | the same two pointers at the field's rank: `u` and `r` are pointers of rank 1 to 5 | null exactly when `stage_state()` and `rhs_register()` are; the rank must be the field's |
| `end_stage(k, bx, dt, q [, mask, mask_comp])` | closes stage `k`: the two-register, 2N, Ketcheson and splitting forms advance `q` here | |
| `finish_step(bx, dt, q [, dq] [, mask, mask_comp])` | closes the step: the Butcher form computes $q^{n+1}$, and `dq` $= \sum_k b_k F_k$ in the same pass if asked | `dq` with any form other than Butcher stops the program; `dq` is 0 where the mask is inactive |
| `stage_state_is_q(k)` | the state of stage `k` is `q` itself: evaluate its residual on `q` | true for every stage of the forms other than Butcher; for the zero rows of A with `stages_on_q` |
| `stages()` | the stages of a step (the flows of a splitting) | |
| `stage_time(k)` | the fraction of the step at which the residual of stage `k` is evaluated: the abscissa $c_k$ | for a splitting, the fractions of the other partition's flows before flow `k`: see [Time-dependent residuals](#time-dependent-residuals) |
| `stage_weight(k)` | the effective Butcher weight $b_k$; for a splitting, the fraction of the step of flow `k` | in every Runge–Kutta form, the low-storage ones included |
| `stage_partition(k)` | the partition flow `k` advances (1 or 2) | 0 for a Runge–Kutta form |
| `stage_euler_fraction(k)` | the forward-Euler step stage `k`'s residual takes, as a fraction of `dt`: $f_k = \max_i \beta_{ik}/\alpha_{ik}$ of the scheme's Shu–Osher form | the step a limiter must make admissible: see [Limiters](#limiters-and-forward-euler-steps). 0 when the scheme is not SSP; a splitting's flow: $c_k$ (0 if negative) |
| `imag_axis_bound()`, `real_axis_bound()` | the largest $\beta$ with $\lvert R(iy)\rvert \le 1$ on $[0, \beta]$ and with $\lvert R(-x)\rvert \le 1$ on $[0, \beta]$: a mode of frequency $\omega$ (a gyration, a wave) is stable for $\omega\,\Delta t \le$ `imag_axis_bound()`, a decay of rate $\lambda$ for $\lambda\,\Delta t \le$ `real_axis_bound()` | from the stability polynomial of the tableau; 0 for a splitting. euler and ssp-22 have no imaginary-axis interval (0) |
| `ssp_coefficient()` | the SSP coefficient $C$ of the scheme, from its Shu–Osher form | 0 when the scheme is not SSP, and for a splitting; $\max_k f_k = 1/C$ |
| `state_aliases_q()` | the form keeps the stage state in `q` (every form but Butcher) | the per-stage `stage_state_is_q(k)` generalises it |
| `lookahead(k, bx, dt, q, u)` | in the window of stage `k`, after its residual: $u = q + \Delta t \sum_{r\le k} a_{k+1,r} R_r$ ($b_r$ at the last stage), the state of stage `k+1` (the step's result at the last stage) | Butcher form, any rank; see [Charge-conserving stage currents](/guide/particles#charge-conserving-stage-currents) |
| `stage_current(k, bx, h, jr, jk)`, `stage_current_admissible()` | $J_k = (H - \sum_{r<k} a_{k+1,r} J_r)/a_{k+1,k}$ ($b$ at the last stage): the stage current whose weighted sum with the earlier ones (`jr`, the caller's registers) is `h`; whether every $a_{k+1,k}$ and $b_s$ is non-zero | Butcher form, rank 5; stops on an inadmissible scheme |
| `initialize_companion(field, lb, ub [, partitions] [, registers] [, components])` | a stepper of the scheme, form and stages on q of the stepper `field`, for another state that advances with the field | see [Companion states](#companion-states) |
| `destroy()` | frees what the stepper owns | no finalizer: call it |

## Fields of rank 1 to 4

`initialize` takes the bounds of a field of any rank from 1 to 5. Below rank 5 it also takes `components`: whether
the field's first dimension holds its components. The field is then viewed as rank 5 with singleton dimensions
([the box](./kernels#fields-of-rank-1-to-4)). The four step procedures are generic over the rank of `q`, and `dq` and
`mask` have `q`'s rank. The box must be built at the field's rank:

```fortran
real(R_P), pointer :: q(:), u(:), r(:)        ! a scalar 1-D field q(0:n+1), one ghost cell each side
type(box)          :: bx

bx = box_interior(lb=[0], ub=[n+1], ngc=1, components=.false.)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[0], ub=[n+1], components=.false.)
call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   if (st%stage_state_is_q(k)) then ; u => q ; else ; call st%stage_state_view(u) ; endif
   call st%rhs_register_view(r)
   call residual(u, r)                       ! yours
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)
```

- **No copy.** A lower-rank `q` reaches the rank-5 procedures as the rank-5 array it is in memory, through an
  explicit-shape dummy, with no copy when it is contiguous. In a device build it must be contiguous.
- **Checked at every call:** the box and the stepper have `q`'s rank, and `q` has the bounds given to `initialize`.
- **The views.** `stage_state()` and `rhs_register()` keep their rank-5 result. A function cannot be generic on the
  rank of its result, so the rank-k pointers are the subroutines `stage_state_view(u)` and `rhs_register_view(r)`,
  with the field's bounds.
- **The stage contract is unchanged at every rank.** Its clauses apply to the cells of the box, whatever the rank.
- **Caller registers** (`registers=`) are for rank-5 fields only.

## The caller's registers

`initialize(..., registers=pool)` takes the registers from the caller instead of allocating them. The stepper checks
two things and stops on a mismatch:
- the array is rank 6 and contiguous, with the field's bounds `lb:ub` on dimensions 1–5;
- its slots are numbered from 1 on dimension 6, and there are at least `registers_needed(scheme, form)` of them.

It cannot check the rest:
- the actual must be a pointer or a target;
- in a device build it must be device memory;
- two steppers in use at the same time need two pools.

The stepper never frees a caller's pool, and does not need it zeroed: no form reads a register before writing it
([C3](/verification/contracts/stage-contract)).

## Stages on q

The state of a stage whose row of A is zero is `q` itself: the first stage of every explicit scheme (in SELAH's
catalogue, the only one). By default the Butcher form still copies it into its stage-state register. With
`stages_on_q=.true.` it does not:

```fortran
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q), stages_on_q=.true.)
! ...
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   if (st%stage_state_is_q(k)) then
      u => q                ! stage 1: fill the ghost cells of q, evaluate the residual on q
   else
      u => st%stage_state()
   endif
   ! ... ghost fill and residual of u into st%rhs_register()
   call st%end_stage(k, bx, dt, q)
enddo
```

During such a stage `stage_state()` is null, so a caller that does not ask `stage_state_is_q(k)` fails at once
instead of reading a stale register. The caller then writes the ghost cells of `q` at stage 1 (contract clause 3),
which the plain Butcher form leaves alone. The solution is bitwise that of the plain Butcher form. On the A30 the option
takes 6% off the integrator's SSP(5,4) step ([Performance](./performance)).

## Companion states

A state that advances with the field in the same step, its residual built inside the field's at each stage (particles
pushed by the fields of each stage, the auxiliary fields of a perfectly matched layer), is the same scheme applied to
another state. A **companion** is a second stepper with the field stepper's scheme, form and stages on q, driven by
the same calls in the same stage window:

```fortran
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb, ub)                      ! the field
call sp%initialize_companion(st, lb=[1, 1], ub=[8, np_max], components=.true.) ! the particles, (8, np_max)
! each step, n particles active, components 1:6 evolving
bp = box_range(lb=[1, 1], lo=[1, 1], hi=[6, n], components=.true.)
call st%begin_step(bx, q) ; call sp%begin_step(bp, p)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q) ; call sp%begin_stage(k, bp, dt, p)
   ! the residuals of both: the particles' from the stage fields, the field's from the stage particles
   call st%end_stage(k, bx, dt, q) ; call sp%end_stage(k, bp, dt, p)
enddo
call st%finish_step(bx, dt, q) ; call sp%finish_step(bp, dt, p)
! then what the caller does on the final state: particle wrap and migration
```

The companion's bounds are a **capacity**; the box of each call is the active part: the particles `1..n` of a pool
whose `n` changes between steps, the components that evolve (charge and mass, components 7:8, stay as they are). The
kernels never touch the rest. A layout with the particles first, `(np, 8)` (the device's), is a rank-2 field with
`components=.false.` whose box ranges over `1:n` and `1:6`. Re-initialise the companion when its shape changes (a
regrid): between steps its registers are scratch. Stepped together, field and companion are bitwise the two stepped as
one state (`selah_test_companion`, every form). The stage window of the two must coincide: a companion's stage state
and residual register are, like the field's, null outside it.

## Masks

With `mask` present the kernels update the active cells only: a cell is active where `mask(..., mask_comp) < 0`
(ADAM's immersed-boundary distance, negative in the fluid; `mask_comp` defaults to the last component). The masked
cells of `q` keep their values through the step, in every form.

The stage state is another matter. The forms other than Butcher keep it in `q`, so its masked cells hold $q^n$. The
Butcher form builds it in a register, and a masked cell of that register is left as it was: the previous stage's,
the previous step's, or the register's initial value. A caller whose residual reads masked cells (an immersed
boundary extrapolating from the solid side, a stencil that reaches into a solid block) asks for the **mask policy**
at initialisation:

```fortran
call st%initialize(scheme, FORM_BUTCHER, lb, ub, masked_q=.true.)
```

and the masked cells of every stage state then hold $q^n$, written by the stage-state kernel in the same pass that
builds the active cells (no extra pass over the field). This is the Runge–Kutta solution of $u' = 0$ in the masked
cells. Bitwise, it is what a caller gets by copying $q^n$ into the masked cells of the stage state itself after
`begin_stage` (ADAM's adapter did so before szaghi/adam#134; `selah_test_stepper_options`). Default `.false.`: the
register's masked cells are left untouched. A caller whose masked cells must evolve within the step (ghost values of
an immersed boundary rebuilt from the fluid at every stage) still writes them itself, between `begin_stage` and its
residual.

## Time-dependent residuals

A residual that depends on time (boundary conditions, sources, forcing) is evaluated at
`t + st%stage_time(k)*dt`:

```fortran
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   ! ... the residual of the stage state at time t + st%stage_time(k)*dt into st%rhs_register()
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)
t = t + dt
```

In a Runge–Kutta form `stage_time(k)` is the abscissa $c_k$. A splitting needs a definition: its sub-flows are
explicit Euler steps of one partition, exact only when that partition's residual does not depend on its own variables.
Time is a variable with $t' = 1$, so it advances with the flows of the partition whose residual does **not** depend
on $t$ (those flows stay exact), and it is frozen through the flows of the other one (which stay exact too, their
residual taken at a fixed time). `stage_time(k)` is therefore the sum of the fractions of the *other* partition's flows
before flow $k$. Strang, partition 1 first, flows $(1, \tfrac12), (2, 1), (1, \tfrac12)$: `stage_time` is
$0, \tfrac12, 1$. If partition 2 depends on time, its flow sees $t + \Delta t/2$; if partition 1 does, its flows see
$t$ and $t + \Delta t$. Either way the composition keeps its order (measured: `selah_test_nonautonomous`).

When **both** partitions' residuals depend on time, no sub-flow is exact and no stage time gives the composition's
order: split the problem so that one partition's residual is autonomous.

## Limiters and forward-Euler steps

A strong-stability-preserving (SSP) scheme is a convex combination of forward-Euler steps. In its Shu–Osher form,
stage $i$ is

$$Y_i = \sum_{k<i} \alpha_{ik} \left( Y_k + \Delta t\, \frac{\beta_{ik}}{\alpha_{ik}} F(Y_k) \right),
\qquad \alpha_{ik} \ge 0,\ \beta_{ik} \ge 0,\ \textstyle\sum_k \alpha_{ik} = 1 ,$$

so the residual of stage $k$ enters forward-Euler steps of $\Delta t\, \beta_{ik}/\alpha_{ik}$ from $Y_k$. A
limiter that makes $Y_k + \tau F(Y_k)$ admissible (positive density and pressure, a bound) for $\tau = f_k \Delta t$,
with $f_k$ = `stage_euler_fraction(k)` $= \max_i \beta_{ik}/\alpha_{ik}$, makes it admissible for every
$0 \le \tau \le f_k \Delta t$ ($Y_k + \tau F$ is then a convex combination of $Y_k$ and that state, and the set
is convex). Every later stage and the step are then convex combinations of admissible states (Zhang & Shu 2010;
Gottlieb, Ketcheson & Shu 2011). The limiter's step is $f_k \Delta t$, not $\Delta t$:

```fortran
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   ! ... the residual of the stage state into st%rhs_register()
   tau = st%stage_euler_fraction(k) * dt
   ! ... the limiter: the state of stage k plus tau times its residual admissible (a source term time-centred within
   !     that forward-Euler step at tau/2)
   call st%end_stage(k, bx, dt, q)
enddo
```

| Scheme | $f_k$ | $C$ |
|---|---|---|
| `euler`, `ssp-22`, `ssp-33` | 1 every stage | 1 |
| `ssp-54` | 0.3918, 0.6631, 0.6631, 0.6631, 0.5844 | 1.508 |
| `ssp-104` | 1/6 every stage | 6 |
| `rk4`, `lsrk-33`, `lsrk-54` | 0 (not SSP: no positivity guarantee) | 0 |

Taking $\Delta t$ for every stage is exact for euler, ssp-22 and ssp-33, and over-limits ssp-54 and ssp-104: it
enforces a forward-Euler step longer than the scheme takes. Taking less than $f_k \Delta t$ loses the guarantee
(`selah_test_stage_euler`: the decay $u' = -\lambda u$ at $\lambda \Delta t = 20$, limited with $f_k \Delta t / 2$,
goes negative; with $f_k \Delta t$ it does not). $f_k$ is the same in every storage form of a scheme: the stages are
the same. The time step itself is bounded by $\Delta t \le C\,\Delta t_{FE}$, `ssp_coefficient()` times the
forward-Euler limit of the residual.

## Forms

| Form | Registers | Schemes |
|------|-----------|---------|
| `FORM_BUTCHER` | s + 1 | every Runge–Kutta scheme |
| `FORM_SHU_OSHER_2R` | 2 | `euler`, `ssp-22`, `ssp-33` |
| `FORM_WILLIAMSON_2N` | 2 | `lsrk-33`, `lsrk-54` |
| `FORM_KETCHESON_2R` | 2 | `ssp-104` |
| `FORM_SPLITTING` | 1 | every splitting |

What each register holds, and where each form does its work: [Concepts](./concepts#storage-forms).

## Verification

| What | Case |
|---|---|
| every scheme reaches its order in every form, through the stepper | [A3](/verification/accuracy/stepper-orders) |
| the stage contract, stages on q, masks, caller registers, the error stops | [C3](/verification/contracts/stage-contract) |
| the forms agree with the Butcher form to round-off | [C2](/verification/contracts/forms) |
| `stage_weight` telescopes the stage fluxes to the mass change | [P4](/verification/properties/conservation) |
| `stage_euler_fraction` and `ssp_coefficient` against the literature; a limiter with them keeps positivity | [P2](/verification/properties/ssp-coefficients) |
| the stepper on an A30 | [C5](/verification/contracts/device) |
| the stepper's own cost: none (0.999 of the same kernels called directly) | [Performance](./performance) |
