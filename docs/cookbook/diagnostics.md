---
title: Diagnostics
---

# Diagnostics

[[toc]]

## Monitor convergence with the step residual

The step residual $dq = \sum_k b_k F_k = (q^{n+1} - q^n)/\Delta t$ measures how much a step changes the solution. A
steady-state solver stops when it is small. `finish_step(..., dq=dq)` writes it in the same pass as the update.

<<< @/examples/snippets/cookbook_diagnostics-steady.F90{fortran}

<<< @/examples/output/diagnostics-steady.txt{text}

`dq` belongs to the Butcher form, which holds every stage residual until the update. The in-place forms (two-register,
2N, Ketcheson) never do, and asking them for `dq` stops the program.

## Check the stage contract while you drive the stepper

Two cheap checks catch most mistakes in a new solver loop:
- the stage state is `q` exactly when `stage_state_is_q(k)`, in which case `stage_state()` is null;
- the residual register is never the stage state.

<<< @/examples/snippets/cookbook_diagnostics-contract.F90{fortran}

<<< @/examples/output/diagnostics-contract.txt{text}

::: warning gfortran 16.2: two pointer functions in one associated()
`associated(st%stage_state(), st%rhs_register())` makes gfortran 16.2 stop with an internal compiler error
(`gfc_conv_procedure_call`). Assign the two pointers to local variables first, as above.
:::

## A request the stepper cannot honour

The stepper stops with a message instead of computing something wrong. The errors it detects:
- a scheme asked for in a form it does not have;
- a splitting without its partitions;
- registers that are too few, numbered from anything but 1, or not of the field's bounds;
- `dq` asked of a form other than the Butcher one.

<<< @/examples/snippets/cookbook_diagnostics-wrong-form.F90{fortran}

<<< @/examples/output/diagnostics-wrong-form.txt{text}

Each of these stops is tested by an expected-failure program (`src/tests/selah_test_xfail_*.F90`), which passes only
when it exits with an error.
