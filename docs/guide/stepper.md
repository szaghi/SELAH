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
| `initialize(scheme, form, lb, ub [, partitions] [, registers] [, stages_on_q] [, components])` | sets the scheme and form for a field of rank `size(lb)` (1 to 5); allocates `registers_needed(scheme, form)` registers with the bounds `lb:ub` through FUNDAL (device memory in a device build), or takes the caller's | stops if the scheme has no such form, or a splitting comes without `partitions`. `partitions(:,p)` = `(v1, v2)`, the components of partition `p`. `components` is required below rank 5 ([fields of rank 1 to 4](#fields-of-rank-1-to-4)). `registers` (see below), for rank 5 only. Calls `destroy` first, so a stepper can be re-initialised |
| `begin_step(bx, q [, mask, mask_comp])` | opens a step: the two-register and Ketcheson forms keep a copy of $q^n$ | |
| `begin_stage(k, bx, dt, q [, mask, mask_comp])` | opens stage `k`: the Butcher form builds its state into its register (unless `stage_state_is_q(k)`) | `k` from 1 to `stages()`, in order |
| `stage_state()` | pointer to the stage state, with the field's bounds | Butcher form only; **null** on a stage evaluated on `q`, and always null in the other forms (their state is `q`). **Null** outside a stage window (before `begin_stage`, after `end_stage`) |
| `rhs_register()` | pointer to the register the caller writes the residual of the active stage into, with the field's bounds | write the cells of the box; for a splitting, the components of `stage_partition(k)`. **Null** outside a stage window (before `begin_stage`, after `end_stage`) |
| `stage_state_view(u)`, `rhs_register_view(r)` | the same two pointers at the field's rank: `u` and `r` are pointers of rank 1 to 5 | null exactly when `stage_state()` and `rhs_register()` are; the rank must be the field's |
| `end_stage(k, bx, dt, q [, mask, mask_comp])` | closes stage `k`: the two-register, 2N, Ketcheson and splitting forms advance `q` here | |
| `finish_step(bx, dt, q [, dq] [, mask, mask_comp])` | closes the step: the Butcher form computes $q^{n+1}$, and `dq` $= \sum_k b_k F_k$ in the same pass if asked | `dq` with any form other than Butcher stops the program; `dq` is 0 where the mask is inactive |
| `stage_state_is_q(k)` | the state of stage `k` is `q` itself: evaluate its residual on `q` | true for every stage of the forms other than Butcher; for the zero rows of A with `stages_on_q` |
| `stages()` | the stages of a step (the flows of a splitting) | |
| `stage_time(k)` | the abscissa $c_k$ (fraction of the step) | 0 for a splitting |
| `stage_weight(k)` | the effective Butcher weight $b_k$; for a splitting, the fraction of the step of flow `k` | in every Runge–Kutta form, the low-storage ones included |
| `stage_partition(k)` | the partition flow `k` advances (1 or 2) | 0 for a Runge–Kutta form |
| `state_aliases_q()` | the form keeps the stage state in `q` (every form but Butcher) | the per-stage `stage_state_is_q(k)` generalises it |
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
| the stepper on an A30 | [C5](/verification/contracts/device) |
| the stepper's own cost: none (0.999 of the same kernels called directly) | [Performance](./performance) |
