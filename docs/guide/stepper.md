---
title: Stepper
---

# Stepper

SELAH's L2 layer is the stepper: one time step that the caller drives one stage at a time. The stepper owns the
registers and the coefficients. The caller owns the solution, the time loop, the ghost cells and the residual. The
stepper never calls back.

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
call st%finish_step(bx, dt, q)        ! the Butcher form also returns the step residual: dq=...
```

Every call takes the box (see [Kernels](./kernels)), and the stage calls take an optional mask. A splitting scheme also
needs the components of its two partitions (`partitions=`). For each flow, `stage_partition(k)` tells which partition
to compute the residual of.

## The stage contract

1. Between `begin_stage(k)` and `end_stage(k)`, the stage state holds the state of stage k in the cells of the box. It
   is `q` itself when `stage_state_is_q(k)`, else `stage_state()`. Every form but the Butcher one keeps the stage
   state in `q` (`state_aliases_q()`); the Butcher form does so for its first stage when asked (below).
2. The stepper never writes a ghost cell, of any register or of `q`.
3. The caller may write the ghost cells of the stage state during that window: ghost exchange, AMR seam fill,
   physical boundaries.
4. `rhs_register()` is written by the caller only, in the cells of the box, and is distinct from the stage state.
5. After `finish_step` the registers are scratch.

`stage_weight(k)` is the effective Butcher weight b_k of stage k, in every Runge–Kutta form. A flux register
accumulates the boundary fluxes of each stage with this weight. `stage_time(k)` is the abscissa c_k.

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

During such a stage `stage_state()` is null, so that a caller that does not ask `stage_state_is_q(k)` fails at once
instead of reading a stale register. The caller then writes the ghost cells of `q` at stage 1 (clause 3), which the
plain Butcher form leaves alone. On the A30 the option takes 6% off the integrator's SSP(5,4) step: 2.552 ms against 2.716 ms (see
[Performance](./performance)).

## Forms

| Form | Registers | Schemes | Where the work happens |
|------|-----------|---------|------------------------|
| `FORM_BUTCHER` | s + 1 | every RK scheme | `begin_stage`: the stage state from q and the earlier residuals; `finish_step`: the update (and dq) |
| `FORM_SHU_OSHER_2R` | 2 | `euler`, `ssp-22`, `ssp-33` | `end_stage`: q = a q⁰ + b q + Δt c R |
| `FORM_WILLIAMSON_2N` | 2 | `lsrk-33`, `lsrk-54` | `end_stage`: ΔQ = A ΔQ + Δt R; q += B ΔQ |
| `FORM_KETCHESON_2R` | 2 | `ssp-104` | `end_stage`: Ketcheson's two-register recipe |
| `FORM_SPLITTING` | 1 | every splitting | `end_stage`: the flow's partition, q += c Δt R |

The stepper allocates the registers through FUNDAL: device memory in a device build. A caller that already holds
registers passes them instead (`registers=`, a rank-6 array with the field's bounds and enough slots).

## Verification

All numbers come from `selah_test_stepper`, 248 checks (2026-11-17). The program drives the stepper the way a solver
would: a 1-D finite-volume advection on 4 blocks of 16 cells, whose ghost fill and residual run on the device in a
device build.

**Order.** On a periodic domain, the central residual of a sine has an exact semi-discrete solution, so the error is
the time integrator's alone. Every scheme in every form reaches its order on the last pair of steps:

| Scheme | Butcher form | Other form |
|---|---|---|
| `euler` | 1.007 | 1.007 (two-register) |
| `ssp-22` | 2.000 | 2.000 (two-register) |
| `ssp-33` | 3.000 | 3.000 (two-register) |
| `ssp-54` | 4.001 | — |
| `ssp-104` | 4.000 | 4.000 (Ketcheson) |
| `rk4` | 4.000 | — |
| `lsrk-33` | 3.000 | 3.000 (2N) |
| `lsrk-54` | 4.001 | 4.001 (2N) |

The splittings, on a two-partition oscillator: `strang` 2.000, `yoshida-4` 4.000, `blanes-moan-s6` 4.000,
`blanes-moan-s10` 6.000.

**Stages on q.** With `stages_on_q=.true.` (the `butcher-q` rows of the test), every scheme keeps its order;
`stage_state_is_q(k)` is true exactly for the zero rows of A (and for every stage of the other forms); `stage_state()`
is null during those stages; the solution in the cells is bitwise that of the plain Butcher form. Flagging every
stage on q instead (a mutation) fails that bitwise check for all 7 multi-stage schemes.

**Contract.** Through whole steps:
- ghost cells that the caller never fills keep a sentinel, in `q` and in every register;
- the residual register is never the stage state;
- `state_aliases_q` is false for the Butcher form only.

**Telescoping.** With inflow on the left and outflow on the right, the change in mass over a step equals
−Δt Σₖ stage_weight(k)·(F_right − F_left)ₖ to round-off, in every Runge–Kutta form, the derived weights of the
two-register, 2N and Ketcheson forms included. With uniform weights in place of `stage_weight`, the check fails for
every scheme whose weights are not uniform.

The tests pass with gfortran 16.2 (serial, debug, OpenMP), nvfortran 26.9, ifx 2026.2, and on an NVIDIA A30.
