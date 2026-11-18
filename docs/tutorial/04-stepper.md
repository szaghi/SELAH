---
title: 4. The stepper
---

# 4. The stepper

Chapter 1 wrote the step by hand. On a block-structured field, the step also has to leave room for what a real solver
does between stages: fill the ghost cells of the stage state (exchange between blocks, AMR seams, physical boundaries)
and evaluate the residual. The **stepper** (`type(stepper)`) owns the registers and the coefficients; you call it
around your own work, one stage at a time. It never calls back.

<<< @/examples/snippets/stepper-step.F90{fortran}

- **`begin_step(bx, q)`** opens the step. Forms that need a copy of $q^n$ take it here.
- **`begin_stage(k, bx, dt, q)`** builds the state of stage `k`. In the Butcher form it goes into the stepper's own
  register, which `stage_state()` returns.
- **Between the two stage calls, the stage is yours:**
  - fill the ghost cells of the stage state;
  - write its residual into `rhs_register()`.
- **`end_stage(k, bx, dt, q)`** closes the stage; the in-place forms (chapter 5) advance `q` here.
- **`finish_step(bx, dt, q)`** closes the step; the Butcher form computes $q^{n+1}$ here.

The setup and the caller's two routines:

<<< @/examples/snippets/stepper-setup.F90{fortran}

<<< @/examples/snippets/stepper-caller.F90{fortran}

- `initialize(scheme, form, lb, ub)` takes the bounds of the field, ghost cells included, and allocates the
  registers with them through FUNDAL: host memory here, device memory in a device build (chapter 9).
- The box describes the field: one component, 16 × 1 × 1 cells per block, one ghost cell, 4 blocks.

## Running it

<<< @/examples/output/stepper.txt{text}

Fourth order through the stepper, as by hand.

::: warning The pointers are valid for one stage
`stage_state()` and `rhs_register()` are pointers into the stepper's registers. Take them after `begin_stage(k)` and
use them before `end_stage(k)`: the next stage uses other registers.
:::

::: tip What you learned
The four calls of a step, and who does what: the stepper builds stage states and updates; the caller fills ghosts and
writes residuals. Reference: [Stepper](/guide/stepper), the [stage contract](/guide/concepts#the-stage-contract).
:::

## A field of any rank

The field so far is rank 5: components, three cell indices and blocks. A 1-D, 2-D or 3-D solver holds fields of rank
1 to 4, and the stepper takes them at their own rank. SELAH views such a field as rank 5 with singleton dimensions,
so the kernels are the same and nothing is copied. The box and the stepper are built from the field's own bounds,
and `components` says whether its first dimension holds components:

<<< @/examples/snippets/ranks-rank1.F90{fortran}

<<< @/examples/snippets/ranks-rank3.F90{fortran}

The step is the one of this chapter, at the field's rank. `stage_state()` and `rhs_register()` keep their rank-5
result. At rank 1 to 4 the pointers come from `stage_state_view(u)` and `rhs_register_view(r)`, with the field's
bounds:

<<< @/examples/snippets/ranks-step.F90{fortran}

The same field declared rank 5 gives the same values, bit for bit:

<<< @/examples/snippets/ranks-rank5.F90{fortran}

<<< @/examples/output/ranks.txt{text}

The rank-1 field is viewed as `(1, i, 1, 1, 1)`, the rank-3 field `w(v, i, j)` as `(v, i, j, 1, 1)`. The whole
mapping is in [Kernels](/guide/kernels#fields-of-rank-1-to-4); the preconditions (a contiguous field, caller registers
for rank 5 only) are in [Stepper](/guide/stepper#fields-of-rank-1-to-4).

Next: [5. Low-storage forms](./05-low-storage).
