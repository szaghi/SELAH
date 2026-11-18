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

### The integrators

Every scheme the API exports, by its constant: pass it to `stepper%initialize` with one of its storage forms (the
field-sized registers each needs in brackets, `registers_needed`). Tableaux, derivations and verification:
[Schemes](/guide/schemes).

| Runge–Kutta | Name | Stages | Order | SSP coefficient | Storage forms | Source |
|---|---|---|---|---|---|---|
| `SCHEME_EULER`   | `'euler'`   | 1  | 1 | 1           | Butcher (2), two-register (2) | forward Euler |
| `SCHEME_SSP_22`  | `'ssp-22'`  | 2  | 2 | 1           | Butcher (3), two-register (2) | Shu & Osher (1988) |
| `SCHEME_SSP_33`  | `'ssp-33'`  | 3  | 3 | 1           | Butcher (4), two-register (2) | Shu & Osher (1988) |
| `SCHEME_SSP_54`  | `'ssp-54'`  | 5  | 4 | 1.508180049 | Butcher (6) | Spiteri & Ruuth (2002) |
| `SCHEME_SSP_104` | `'ssp-104'` | 10 | 4 | 6           | Butcher (11), Ketcheson (2) | Ketcheson (2008) |
| `SCHEME_RK4`     | `'rk4'`     | 4  | 4 | —           | Butcher (5) | Kutta (1901) |
| `SCHEME_LSRK_33` | `'lsrk-33'` | 3  | 3 | —           | Butcher (4), Williamson 2N (2) | Williamson (1980) |
| `SCHEME_LSRK_54` | `'lsrk-54'` | 5  | 4 | —           | Butcher (6), Williamson 2N (2) | Carpenter & Kennedy (1994) |

| Splitting | Name | Flows | Order | Storage form | Source |
|---|---|---|---|---|---|
| `SCHEME_STRANG`          | `'strang'`          | 3  | 2 | splitting (1) | Strang (1968) |
| `SCHEME_YOSHIDA_4`       | `'yoshida-4'`       | 7  | 4 | splitting (1) | Yoshida (1990) |
| `SCHEME_BLANES_MOAN_S6`  | `'blanes-moan-s6'`  | 13 | 4 | splitting (1) | Blanes & Moan (2002), S6 |
| `SCHEME_BLANES_MOAN_S10` | `'blanes-moan-s10'` | 21 | 6 | splitting (1) | Blanes & Moan (2002), S10 |

### Names of forms and lists

A scheme or a form is named by a character constant; an unknown name returns a form of 0 stages (0 flows).

| Constant | Value |
|---|---|
| `FORM_BUTCHER` | `'butcher'` |
| `FORM_SHU_OSHER_2R` | `'shu-osher-2r'` |
| `FORM_WILLIAMSON_2N` | `'williamson-2n'` |
| `FORM_KETCHESON_2R` | `'ketcheson-2r'` |
| `FORM_SPLITTING` | `'splitting'` |
| `RK_SCHEMES` | the 8 Runge–Kutta names, `character(7)`, blank-padded: `trim` them |
| `SPLITTING_SCHEMES` | the 4 splitting names, `character(15)`, blank-padded |

Which scheme in which form, and the SSP coefficients: [Schemes](/guide/schemes).

## L1: kernels

The stage operations, over the part of a block-structured field a [[box]] describes, or over a flat array.

| Module | What |
|--------|------|
| [[selah_box]] | the `box`, its layouts `LAYOUT_V_FIRST` and `LAYOUT_V_LAST`, `box_interior` (rank 5, or a field of rank 1 to 4), `box_range`, `box_components`, `box_dim` (a field of rank k as rank 5) |
| [[selah_kernels]] | the generic kernels `selah_copy`, `selah_axpy`, `selah_lincomb`, `selah_stage_state`, `selah_step_update`, `selah_shu_osher_2r_stage`, `selah_williamson_2n_stage` |
| [[selah_kernels_box]], [[selah_kernels_flat]] | their specific procedures: the box forms (rank-5 arrays, an optional mask) and the flat forms (rank 1) |

Every output is `intent(inout)`, must not alias an input, and on a device is contiguous device memory: the rules and
the shortcuts are in [Kernels](/guide/kernels).

## L2: stepper

| Module | What |
|--------|------|
| [[selah_stepper]] | the [[stepper]] type: one time step driven stage by stage, the caller filling the ghost cells and writing each residual; on fields of rank 1 to 5 (`stage_state_view`, `rhs_register_view` below rank 5) |

The order of the calls and every precondition: [Stepper](/guide/stepper). The components of `type(stepper)` are
public but are not part of the API.
