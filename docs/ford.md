---
project: SELAH
summary: Stage-wise Evolution Library for Accelerators and HPC: time-evolution operators for block-structured field solvers
project_url: https://szaghi.github.io/SELAH
author: Stefano Zaghi
email: stefano.zaghi@gmail.com
src_dir: ../src
exclude_dir: ../src/third_party
             ../src/tests
             ../src/bench
include: ../src/third_party/FUNDAL/src/lib
project_github: https://github.com/szaghi/SELAH
docmark: <
display: public
proc_internals: false
coloured_edges: true
graph_maxdepth: 6
---

`use selah` gives the whole public API, and nothing else: the [[selah]] module re-exports every name below from the
module that defines it (its page lists them all, each linked to its definition). The pages of the modules carry every
procedure's signature, arguments and preconditions; this page is the map. Worked examples: the
[Tutorial](/tutorial/), chapters 1 (the kernels by hand), 2 (the catalogue), 3 (a field and a box) and 4 (the stepper).

## Kinds

`R_P` is the working real kind, chosen at build time: R8P by default, R4P with `-D_R_P_IS_R4P` (the `*-r4` modes),
see [Precision](/guide/precision). `I4P`, `I8P`, `R4P` and `R8P` are re-exported from PENF; [[selah_kind_parameters]]
ties SELAH's kind to PENF's.

## L0: coefficients

The coefficients of a scheme, built once, before any field exists.

| Module | What |
|--------|------|
| [[selah_tableaux]] | the coefficient forms (Butcher, Shu–Osher, two-register, Williamson 2N) and the conversions between them |
| [[selah_schemes]] | the Runge–Kutta catalogue: `butcher(scheme)`, `shu_osher`, `shu_osher_2r`, `williamson_2n` |
| [[selah_splittings]] | the splitting schemes, `splitting(scheme[, first])`, and their `splitting_composition` |
| [[selah_storage]] | the storage forms and `registers_needed(scheme, form)`, the field-sized buffers a step needs |

A scheme or a form is named by a character constant; an unknown name returns a form of 0 stages (0 flows).

| Constant | Value | | Constant | Value |
|---|---|---|---|---|
| `SCHEME_EULER` | `'euler'` | | `SCHEME_STRANG` | `'strang'` |
| `SCHEME_SSP_22` | `'ssp-22'` | | `SCHEME_YOSHIDA_4` | `'yoshida-4'` |
| `SCHEME_SSP_33` | `'ssp-33'` | | `SCHEME_BLANES_MOAN_S6` | `'blanes-moan-s6'` |
| `SCHEME_SSP_54` | `'ssp-54'` | | `SCHEME_BLANES_MOAN_S10` | `'blanes-moan-s10'` |
| `SCHEME_SSP_104` | `'ssp-104'` | | `FORM_BUTCHER` | `'butcher'` |
| `SCHEME_RK4` | `'rk4'` | | `FORM_SHU_OSHER_2R` | `'shu-osher-2r'` |
| `SCHEME_LSRK_33` | `'lsrk-33'` | | `FORM_WILLIAMSON_2N` | `'williamson-2n'` |
| `SCHEME_LSRK_54` | `'lsrk-54'` | | `FORM_KETCHESON_2R` | `'ketcheson-2r'` |
| `RK_SCHEMES` | the 8 Runge–Kutta names, `character(7)`, blank-padded: `trim` them | | `FORM_SPLITTING` | `'splitting'` |
| `SPLITTING_SCHEMES` | the 4 splitting names, `character(15)`, blank-padded | | | |

Which scheme in which form, and the SSP coefficients: [Schemes](/guide/schemes).

## L1: kernels

The stage operations, over the part of a block-structured field a [[box]] describes, or over a flat array.

| Module | What |
|--------|------|
| [[selah_box]] | the `box`, its layouts `LAYOUT_V_FIRST` and `LAYOUT_V_LAST`, `box_interior`, `box_components` |
| [[selah_kernels]] | the generic kernels `selah_copy`, `selah_axpy`, `selah_lincomb`, `selah_stage_state`, `selah_step_update`, `selah_shu_osher_2r_stage`, `selah_williamson_2n_stage` |
| [[selah_kernels_box]], [[selah_kernels_flat]] | their specific procedures: the box forms (rank-5 arrays, an optional mask) and the flat forms (rank 1) |

Every output is `intent(inout)`, must not alias an input, and on a device is contiguous device memory: the rules and
the shortcuts are in [Kernels](/guide/kernels).

## L2: stepper

| Module | What |
|--------|------|
| [[selah_stepper]] | the [[stepper]] type: one time step driven stage by stage, the caller filling the ghost cells and writing each residual |

The order of the calls and every precondition: [Stepper](/guide/stepper). The components of `type(stepper)` are
public but are not part of the API.
