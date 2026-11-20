---
title: Kernels
icon: kernels
---

# Kernels

SELAH's L1 layer holds the kernels: pointwise stage operations over the caller's arrays. They own no memory, call
nothing back, and touch nothing outside the cells they are given.

## The box

A field is a rank-5 array (components, three cell indices with ghost cells, blocks) in one of two layouts:

| Layout | Index order | Contiguous | Used by |
|--------|-------------|------------|---------|
| `LAYOUT_V_FIRST` | `(v, i, j, k, b)` | the components of a cell | ADAM's CPU fields |
| `LAYOUT_V_LAST`  | `(b, i, j, k, v)` | the same cell of every block | ADAM's FNL (device) fields |

A `box` holds the lower bounds of the allocation (`lb`, ghosts included), the range the kernels loop over (`lo:hi`)
and the layout. `box_interior(layout, nv, ni, nj, nk, ngc, blocks_number)` gives the interior cells of the active
blocks. The pool may hold more blocks than are active, and the kernels never touch the slots past `blocks_number`.
`box_components(bx, v1, v2)` restricts a box to the components of one partition of a splitting scheme.

A box is passed at every call, never stored, because the active blocks change at every regrid.

### Fields of rank 1 to 4

A field of lower rank is a rank-5 field with singleton dimensions, so the kernels and the stepper keep one rank-5
implementation. Dimension d of the field is dimension d of the rank-5 view when its first dimension holds the
components, and d + 1 when the field is scalar (`box_dim(bx, d)`):

| Field | As rank 5 `(v, i, j, k, b)` |
|-------|-----------------------------|
| `(i)`, scalar | `(1, i, 1, 1, 1)` |
| `(v, i)` | `(v, i, 1, 1, 1)` |
| `(v, i, j)` | `(v, i, j, 1, 1)` |
| `(v, i, j, k)` | `(v, i, j, k, 1)` |
| `(i, j)`, `(i, j, k)`, `(i, j, k, b)`, scalar | `(1, i, j, 1, 1)`, `(1, i, j, k, 1)`, `(1, i, j, k, b)` |

There are two constructors. Both take the `components` flag explicitly, so the mapping is never guessed:
- `box_interior(lb, ub, ngc, components)` gives the interior of a field allocated `lb:ub`, ghosts included. It drops
  `ngc` cells off each end of every dimension except the components, and keeps all the components.
- `box_range(lb, lo, hi, components)` takes any loop range `lo:hi` of a field whose lower bounds are `lb`.

`lb`, `ub`, `lo` and `hi` have the rank of the field. Fields of rank below 5 use `LAYOUT_V_FIRST`. The kernels take
the rank-5 view: pass a lower-rank field to the stepper, which takes it at its own rank
([Stepper](./stepper#fields-of-rank-1-to-4)).

## Operations

| Generic | Operation |
|---------|-----------|
| `selah_copy`                | y = x |
| `selah_axpy`                | y = y + α x |
| `selah_lincomb`             | y = α x + β y |
| `selah_stage_state`         | u = q + Σₖ Δt coef(k) regs(k): the state of a stage of a Butcher form (n = 0: a copy); with `masked_q`, u = q in the inactive cells, in the same pass |
| `selah_step_update`         | q = q + Δt Σₖ coef(k) regs(k), with the step residual dq = Σₖ coef(k) regs(k) when asked |
| `selah_shu_osher_2r_stage`  | q = a q⁰ + b q + Δt c r: a stage of a two-register Shu–Osher form |
| `selah_williamson_2n_stage` | ΔQ = A ΔQ + Δt r; q = q + B ΔQ: a stage of a Williamson 2N form |

Each name is generic. The box form takes the box first and an optional mask. The flat form works on whole rank-1
arrays with no ghosts and no mask, for particle data or perfectly-matched-layer slabs.

**Mask.** With `mask` present, a cell is active where `mask(…, mask_comp) < 0`. This is ADAM's immersed-boundary
distance, negative in the fluid. `mask_comp` defaults to the last component, ADAM's all-solids summary. Inactive
cells are left untouched, except the step residual, which is 0 there.

**Summation order.** A sum over registers is formed privately and stored once, as ADAM's FNL kernels do. The results
are therefore bitwise those of ADAM's device path (see Verification).

**Rules the kernels do not check.**
- **No aliasing.** An output must not alias an input: `selah_axpy(bx, a, y, y)` is undefined. The kernels are
  written for distinct arrays.
- **Bounds.** A box kernel's arrays must have the box's lower bounds, `bx%lb` (`regs` with its slots on dimension 6,
  from 1). A flat kernel loops over the size of its output, so every array must have at least that size.
- **Contiguity on a device.** A device kernel receives its arrays by address (`deviceptr`), so they must be contiguous
  device memory, the coefficient arrays included. On the host a strided section such as `bt%a(s,:)` works. On a GPU,
  pass a contiguous copy: the stepper keeps the transpose of A for this reason.
- **`mask_comp` without `mask`** is ignored.

**Two shortcuts.** `selah_stage_state` with `n = 0` is a copy, which runs at the speed of a copy. The 2N stage with
`a = 0` (the first stage of every 2N scheme) writes $\Delta Q = \Delta t\, r$ without reading $\Delta Q$: its
previous content, NaN included, cannot reach the step.

**Backends.** One source serves the host (plain loops, or OpenMP with `-fopenmp`) and OpenACC,
through FUNDAL's `fundal.H` macros in the directive style of ADAM's FNL kernels. Each loop body is an `.INC` file,
included unmasked and once per layout masked. OpenMP offload has the directives but no build mode yet
([Installation](./installation#backends)).

## Verification

| What | Case |
|---|---|
| every kernel, cell by cell: values in the box, every other cell bitwise untouched (348 checks) | [C3](/verification/contracts/stage-contract) |
| ADAM's FNL kernels bitwise, ADAM's CPU loops within 4 ε (1.00 ε measured) | [C1](/verification/contracts/adam-equivalence) |
| on an NVIDIA A30 | [C5](/verification/contracts/device) |
| pinned values | [C4](/verification/contracts/golden) |
| timings against ADAM | [Performance](./performance) |
