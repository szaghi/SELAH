---
title: Kernels
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

## Operations

| Generic | Operation |
|---------|-----------|
| `selah_copy`                | y = x |
| `selah_axpy`                | y = y + α x |
| `selah_lincomb`             | y = α x + β y |
| `selah_stage_state`         | u = q + Σₖ Δt coef(k) regs(k): the state of a stage of a Butcher form |
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

**Backends.** One source serves the host (plain loops, or OpenMP with `-fopenmp`), OpenACC and OpenMP offload,
through FUNDAL's `fundal.H` macros in the directive style of ADAM's FNL kernels. Each loop body is an `.INC` file,
included unmasked and once per layout masked. In a device build every array argument must be device resident, the
coefficients included.

## Verification

All numbers come from the test programs (2026-11-17).

**Cell by cell** (`selah_test_kernels`, 270 checks). The fields have 3 components and 4×3×2 interior cells with 2
ghost cells, in a pool of 6 blocks, in both layouts. Every kernel runs on 3 and then 5 active blocks, unmasked and
masked, and on a component range. The checks:

- every cell of the box holds the reference value, within 16 ε;
- every other cell (ghosts, inactive pool slots, components outside the range, masked-out cells) is **bitwise
  untouched**;
- the step residual is exactly 0 on masked cells.

A mutant kernel that writes one ghost layer is caught. The tests pass with gfortran 16.2 (serial, debug, and
OpenMP with 4 threads), nvfortran 26.9 and ifx 2026.2.

**Against ADAM** (`selah_test_kernels_adam`).

| Path | Oracle | Result |
|------|--------|--------|
| FNL (`LAYOUT_V_LAST`) | ADAM's `adam_fnl_rk_kernels.F90`, byte-exact copy | an SSP(5,4) step (stage states, solution) and a two-register SSP(3,3) step, unmasked and masked: **bitwise equal** on gfortran, nvfortran and ifx |
| CPU (`LAYOUT_V_FIRST`) | ADAM's `rk_object` loops, transcribed | first stage and two-register stages bitwise; later stages, dq and solution within **1.00 ε** |

ADAM's CPU loops sum term by term in place, a different order from its FNL kernels. Under FMA contraction (nvfortran)
even a one-term stage can differ: ADAM's in-place update contracts into one rounding, while SELAH's private sum rounds
twice. A mutant that re-associates the stage product is caught by the bitwise checks.

**Device.** The OpenACC build (`tests-nvf-oac`) compiles 116 device kernels. They are not run yet: the device tests
and the benchmarks against ADAM's kernels on the GPU are milestone M3.
