---
title: 10. Production
---

# 10. Production

Three things a production solver adds to chapter 4:
- registers from its own memory pool;
- no copy for the first stage;
- the step residual, for convergence monitoring.

<<< @/examples/snippets/production-setup.F90{fortran}

- **Your registers.** `registers=` takes a rank-6 array with the field's bounds and enough slots, numbered from 1
  (`registers_needed` says how many). The stepper checks the bounds and stops on a mismatch. In a device build the pool
  must be device memory. The stepper never frees it, and it does not need it zeroed: no form reads a register before
  writing it.
- **`stages_on_q=.true.`** The state of stage 1 of an explicit scheme is $q^n$ itself. With this option the Butcher
  form skips the copy into its register: `stage_state_is_q(1)` is true, and you evaluate that stage's residual on `q`.
  On the A30 this takes 6% off the integrator's step ([Performance](/guide/performance)).

<<< @/examples/snippets/production-step.F90{fortran}

- **`finish_step(..., dq=dq)`** writes the step residual $dq = \sum_k b_k F_k$ in the same pass as the update, so the
  monitoring costs one write and no extra read of the registers. Only the Butcher form has it: the in-place forms
  never hold all the stage residuals at once.

## Running it

<<< @/examples/output/production.txt{text}

`dq` is exactly the rate of the step, $(q^{n+1} - q^n)/\Delta t$, to round-off.

::: warning Ghost cells of q at stage 1
With `stages_on_q`, stage 1's ghost fill writes the ghost cells of `q` (the stage state *is* `q`). That is allowed:
ghost cells are the caller's (contract clause 3). A plain Butcher form never shows `q` to the ghost fill.
:::

::: tip What you learned
A production loop: your pool, stage 1 on `q`, the fused step residual.
Reference: [Stepper](/guide/stepper), [Cookbook](/cookbook/).
:::

That is the whole API surface a solver needs. The [Cookbook](/cookbook/) has recipes for specific tasks, and
[Verification](/verification/) shows the evidence behind every claim on these pages.
