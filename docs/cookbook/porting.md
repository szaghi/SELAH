---
title: Porting from ADAM
---

# Porting from ADAM

SELAH was carved from ADAM's Runge–Kutta layer. Its FNL path computes exactly what ADAM's FNL kernels compute, bit for
bit ([C1](/verification/contracts/adam-equivalence)), so a port can be checked against ADAM's own regression goldens.
This page maps ADAM's calls to SELAH's, in preparation for wiring SELAH into ADAM.

## The step

ADAM's FNL backend (`adam_fnl_rk_kernels`, the version copied verbatim into `src/tests/adam_fnl_rk_kernels.F90`) runs
an SSP step as follows:

```fortran
call rk_initialize_stages_dev(ni, nj, nk, ngc, nv, nb, q_gpu, q_rk_gpu)            ! q copied into every stage buffer
do s=1, nrk
   call rk_compute_stage_dev(ni, nj, nk, ngc, nv, nb, s, dt, alph, q_rk_gpu=q_rk_gpu)  ! stage s in place
   ! ... ghost fill, residual into dq_gpu ...
   call rk_assign_stage_dev(ni, nj, nk, ngc, nv, nb, s, q_gpu=dq_gpu, q_rk_gpu=q_rk_gpu) ! copy the residual in
enddo
call rk_update_q_dev(ni, nj, nk, ngc, nv, nb, nrk, dt, beta, q_rk_gpu=q_rk_gpu, q_gpu=q_gpu)
call compute_rk_ssp_residual_dev(ni, nj, nk, ngc, nv, nb, nrk, beta, q_rk_gpu, dq_gpu)  ! FLUME: rereads the stages
```

The same step with SELAH:

```fortran
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q_gpu), ub=ubound(q_gpu), stages_on_q=.true.)
! ... every step:
call st%begin_step(bx, q_gpu)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q_gpu)
   u => q_gpu ; if (.not. st%stage_state_is_q(k)) u => st%stage_state()
   ! ... ghost fill of u, residual of u into st%rhs_register() ...
   call st%end_stage(k, bx, dt, q_gpu)
enddo
call st%finish_step(bx, dt, q_gpu, dq=dq_gpu)
```

| ADAM | SELAH | What changes |
|---|---|---|
| `q_rk_gpu(..., 1:nrk)`, allocated by the RK object | `registers_needed(scheme, form)` registers, the stepper's or `registers=` | one more slot (the stage state is separate from the residuals), or 2 in a low-storage form |
| `rk_initialize_stages_dev` | nothing | the copy of q into every buffer is gone: each stage state is built from q |
| `rk_compute_stage_dev(s)` (in place) | `begin_stage(k)` | builds the state from q and the earlier residuals; with `stages_on_q`, stage 1 is q itself |
| `rk_assign_stage_dev(s)` | nothing | the residual is written directly into `rhs_register()` |
| `rk_update_q_dev` | `finish_step` | the same sum, the same order, bitwise |
| `compute_rk_ssp_residual_dev` (FLUME) | `finish_step(dq=)` | fused into the update: no reread of the stages |
| `rk_compute_stage_ls_dev` (2-register SSP) | `FORM_SHU_OSHER_2R`, `end_stage(k)` | the same `ark`, `brk`, `crk` coefficients, bitwise |
| `phi_gpu` (IB distance, last component) | `mask=`, `mask_comp=` | the same sign convention (active where < 0) and the same default component |
| the ghost fill into the stage buffers (seams) | the ghost fill of `stage_state()` or `q` | inside the stage window: the [stage contract](/guide/concepts#the-stage-contract) |

The memory traffic per SSP(5,4) step falls from 53 field accesses to 26, and the measured time to 0.52 of FLUME's on an
A30 ([Performance](/guide/performance)).

## The layout and the box

ADAM's device fields are `(b, i, j, k, v)`, allocated `(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv)`: that is
`LAYOUT_V_LAST`, and `box_interior(LAYOUT_V_LAST, nv, ni, nj, nk, ngc, blocks_number)` describes it. The box's block
range is the *active* prefix of the pool (`blocks_number`), so a pool larger than the active set costs nothing.

## Summation order on the CPU

ADAM's CPU backend (`rk_object`) sums the stage terms in a different order from its FNL backend (ADAM#128, D6). SELAH
follows the FNL order everywhere. Ported CPU results therefore move by round-off, within the envelope SELAH measures
against ADAM's CPU loops (4 ε asserted, 1.00 ε observed; [C1](/verification/contracts/adam-equivalence)). The FNL
results do not move.
