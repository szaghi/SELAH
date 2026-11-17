---
title: Concepts
icon: kernels
---

# Concepts

What a solver needs to know to use SELAH correctly: the shape of a field, what a kernel touches, how the work of a
step is split between the stepper and the caller, and what each storage form keeps.

## The field

SELAH works on the fields of block-structured solvers. A field is a rank-5 array: components, the three cell indices
with their ghost cells, and the blocks of a pool. It has one of two layouts:

![The two layouts](/guide/conv-layouts.png)

| Layout | Index order | Contiguous in memory | Used by |
|---|---|---|---|
| `LAYOUT_V_FIRST` | `(v, i, j, k, b)` | the components of a cell | CPU fields (ADAM's CPU backend) |
| `LAYOUT_V_LAST` | `(b, i, j, k, v)` | the same cell of every block | device fields (ADAM's FNL backend): consecutive GPU threads read consecutive blocks |

The cell indices run from `1-ngc` to `n+ngc`, ghost cells included. Components and blocks start at 1. The pool may hold
more blocks than are active.

## The box

A kernel updates a **box**: the interior cells of the active blocks, every component or a range of them. It writes
nothing else.

![The box](/guide/conv-box.png)

- `box_interior(layout, nv, ni, nj, nk, ngc, blocks_number)` builds the box of the first `blocks_number` blocks.
- `bx%lb` holds the array's lower bounds, which the kernels need to address it (ghosts included). `bx%lo:bx%hi` is
  the loop range.
- `box_components(bx, v1, v2)` narrows a box to components `v1..v2`: a partition of a splitting, or the passive
  scalars of a solver.
- A box is passed at every call and never stored, because the active blocks change at every regrid.

Arrays reach the kernels with their own lower bounds: every dummy is declared `q(bx%lb(1):, …, bx%lb(5):)`. An
assumed-shape dummy without bounds would renumber a field with ghost cells from 1 and shift the interior.

## Masks

Every box kernel and every stepper call takes an optional `mask`: an array with the field's bounds and one or more
components. A cell is **active where the mask is negative**, as with a level-set distance that is negative in the
fluid. `mask_comp` picks the component; the default is the last (ADAM's "all solids" summary). Inactive cells are left
untouched, except the step residual `dq`, which is 0 there.

## The stage contract

The stepper and the caller share the work of a step. The stepper builds stage states and updates. The caller fills
ghost cells and evaluates residuals. The **stage contract** is the agreement that makes this safe:

![The stage contract](/guide/conv-contract.png)

1. Between `begin_stage(k)` and `end_stage(k)`, the stage state holds the state of stage $k$ in the cells of the box:
   `q` itself when `stage_state_is_q(k)`, else `stage_state()`.
2. The stepper never writes a ghost cell, of any register or of `q`.
3. The caller may write the ghost cells of the stage state during that window: ghost exchange, AMR seam fill,
   physical boundaries.
4. `rhs_register()` is written by the caller only, in the cells of the box, and is never the stage state.
5. After `finish_step` the registers are scratch.

`stage_weight(k)` is the effective Butcher weight $b_k$ of stage $k$ in every Runge–Kutta form, the weight a flux
register needs ([Tutorial 8](/tutorial/08-flux-register)). Every clause is tested ([C3](/verification/contracts/stage-contract)).

## Storage forms

A scheme can run in several storage forms, which keep different registers and do the work at different calls:

![Register maps](/guide/conv-registers.png)

| Form | Registers | Schemes | The stage state | Where the work happens |
|---|---|---|---|---|
| `FORM_BUTCHER` | s + 1 | every Runge–Kutta scheme | its own register (or q, with `stages_on_q`, for the zero rows of A) | `begin_stage`: the state from q and the earlier residuals; `finish_step`: the update (and dq) |
| `FORM_SHU_OSHER_2R` | 2 | euler, ssp-22, ssp-33 | q | `begin_step`: a copy of $q^n$; `end_stage`: $q = a\,q^n + b\,q + \Delta t\,c\,R$ |
| `FORM_WILLIAMSON_2N` | 2 | lsrk-33, lsrk-54 | q | `end_stage`: $\Delta Q = A\,\Delta Q + \Delta t R$; $q \mathrel{+}= B\,\Delta Q$ |
| `FORM_KETCHESON_2R` | 2 | ssp-104 | q | `begin_step`: a copy; `end_stage`: Ketcheson's recipe |
| `FORM_SPLITTING` | 1 | every splitting | q | `end_stage`: the flow's partition, $q_p \mathrel{+}= c_k\,\Delta t\,R$ |

The low-storage forms give the same solution as the Butcher form, to round-off ([C2](/verification/contracts/forms)).
Only the Butcher form can return the step residual `dq`, because only it holds every stage residual until the update.

## Ownership

| What | Owner |
|---|---|
| the solution `q`, its ghost cells, the time loop, the residual, the mask | the caller |
| the registers | the stepper (allocated through FUNDAL, freed by `destroy`), or the caller (`registers=`: the stepper never frees them) |
| the coefficients | the stepper (device copies of A and b in a device build) |

The stepper has no finalizer: call `destroy`. Its components are public for the tests' sake, but they are not part of
the API, and code should not read them.
