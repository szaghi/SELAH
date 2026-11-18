---
title: C3 Stage contract
icon: selah
---

# C3 Stage contract

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C3 | Contracts | C | CPU (gfortran, nvfortran, ifx; R8P and R4P); GPU (A30) | no | `selah_test_stepper`, `selah_test_stepper_options`, `selah_test_kernels`, `selah_test_xfail_*` |

## Case and scope

A solver programs against the [stage contract](/guide/concepts#the-stage-contract). This case checks every clause, and
the options around it, as a solver would use them:

1. Between `begin_stage(k)` and `end_stage(k)` the stage state holds stage $k$: `q` when `stage_state_is_q(k)`,
   else `stage_state()`.
2. The stepper never writes a ghost cell.
3. The caller may write the ghost cells of the stage state.
4. `rhs_register()` is distinct from the stage state, and null outside the window, as `stage_state()` is.
5. After `finish_step` the registers are scratch.

## Setup

| Clause or option | How it is checked | Program |
|---|---|---|
| no ghost cell written (2) | every register and `q` start with a sentinel in their ghost layers; the layers the caller never fills keep it bit for bit through whole steps, every scheme and form | `selah_test_stepper` |
| a kernel writes its box only | every kernel, both layouts, masked and unmasked, on a component range: every cell outside the box (ghosts, pool slots past the active blocks, other components, inactive cells) bitwise untouched | `selah_test_kernels` (316 checks) |
| residual register distinct (4) | `associated(rhs_register(), stage_state)` false at every stage | `selah_test_stepper` |
| no register outside a window (4) | `rhs_register()` and `stage_state()` not associated after `begin_step` and after every `end_stage`, every scheme and form | `selah_test_stepper` |
| `state_aliases_q`, `stage_state_is_q` (1) | true exactly for the aliasing forms, and for the zero rows of A with `stages_on_q`; `stage_state()` null on such a stage | `selah_test_stepper` |
| `stages_on_q` | the solution bitwise that of the plain Butcher form, every scheme | `selah_test_stepper` |
| `registers=` | NaN-filled caller registers give the solution of the own registers bitwise, in the Butcher, two-register, 2N and Ketcheson forms: no form reads a register before writing it | `selah_test_stepper_options` |
| masks | solid blocks of `q` bitwise untouched through whole steps in four forms; `mask_comp` absent = the last component | `selah_test_stepper_options`, `selah_test_kernels` |
| `stage_time`, `finish_step(dq=)` | $c_k$ (0 for a splitting); $dq = (q^{n+1} - q^n)/\Delta t$ to round-off | `selah_test_stepper_options` |
| error stops | `dq` for a non-Butcher form; registers numbered from 0; rk4 in a two-register form: each program must stop with its message | `selah_test_xfail_dq_form`, `_registers`, `_form` |

## Program

```bash
./exe/selah_test_stepper ; ./exe/selah_test_stepper_options ; ./exe/selah_test_kernels
bash scripts/run_tests.sh       # runs the xfail programs and passes them when they exit with an error
```

## Results

All checks pass on every compiler, in both precisions, and on the A30: 248 in `selah_test_stepper`, 23 in
`selah_test_stepper_options`, 316 in `selah_test_kernels`. The three expected failures stop with their messages.

## Findings

- No form reads a register before writing it: the stepper accepts uninitialised registers. Before the documentation
  audit of 2026-11-17 the 2N form's first stage computed $0 \cdot \Delta Q$, and a NaN in a caller's register poisoned
  the step. The 2N stage with $A = 0$ now writes $\Delta Q = \Delta t R$ without reading it.
- `initialize(registers=)` now refuses registers whose slots are not numbered from 1 or whose bounds are not the
  field's. Such registers used to be indexed one slot off without a message.
