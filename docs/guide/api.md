---
title: API
---

# API

`use selah` gives everything below. The full generated reference is under [API Reference](/api/).

## Kinds

`R_P` is the working real kind, chosen at build time: R8P by default, R4P with `-D_R_P_IS_R4P` (the `*-r4` modes).
`I4P`, `I8P`, `R4P` and `R8P` are re-exported from PENF.

## L0: coefficients

### Types

| Type | Components |
|------|------------|
| `butcher_tableau`    | `stages`, `order`, `a(s,s)`, `b(s)`, `c(s)`, `b_hat(s)` (embedded weights, unallocated when none), `ssp_coefficient` |
| `shu_osher_form`     | `stages`, `alpha(s+1,s)`, `beta(s+1,s)` |
| `shu_osher_2r_form`  | `stages`, `a(s)`, `b(s)`, `c(s)` |
| `williamson_2n_form` | `stages`, `a(s)`, `b(s)` |
| `splitting_composition` | `flows`, `order`, `partition(flows)` (1 or 2), `coefficient(flows)` |

### Catalogue

| Function | Returns |
|----------|---------|
| `butcher(scheme)`       | the Butcher tableau of any RK scheme |
| `shu_osher(scheme)`     | the Shu–Osher form of an SSP scheme |
| `shu_osher_2r(scheme)`  | the two-register form of `euler`, `ssp-22`, `ssp-33` |
| `williamson_2n(scheme)` | the 2N form of `lsrk-33`, `lsrk-54` |
| `splitting(scheme[, first])` | a splitting scheme, `first` the partition of the first flow (default 1) |

An unknown name returns a form of 0 stages (0 flows). The names are the constants `SCHEME_*`; the lists of all of
them are `RK_SCHEMES` and `SPLITTING_SCHEMES`.

### Conversions and queries

| Procedure | Purpose |
|-----------|---------|
| `butcher_from_shu_osher(so[, order])`    | Shu–Osher → Butcher (the summation order of ADAM's `ssp_tableau`) |
| `butcher_from_shu_osher_2r(r[, order])`  | two-register → Butcher |
| `butcher_from_williamson_2n(w[, order])` | 2N → Butcher |
| `shu_osher_from_2r(r)`                   | two-register → Shu–Osher |
| `is_explicit(bt)`                        | the tableau is strictly lower triangular |
| `ssp_coefficient(so)`                    | $\min \alpha_{ik}/\beta_{ik}$, 0 if the form is not SSP |
| `palindrome(first, half, order)`         | a symmetric composition from its half sequence |
| `registers_needed(scheme, form)`         | field-sized buffers besides the solution, -1 if the scheme has no such form |

## L1: kernels

| Name | Purpose |
|------|---------|
| `box`, `LAYOUT_V_FIRST`, `LAYOUT_V_LAST` | the part of a rank-5 field a kernel updates, and its layout |
| `box_interior(layout, nv, ni, nj, nk, ngc, blocks_number)` | the interior cells of the active blocks |
| `box_components(bx, v1, v2)`, `component_dim(bx)` | a component range (a partition); the component dimension |
| `selah_copy([bx,] x, y [, mask, mask_comp])` | y = x |
| `selah_axpy([bx,] alpha, x, y [, mask, mask_comp])` | y = y + alpha x |
| `selah_lincomb([bx,] alpha, x, beta, y [, mask, mask_comp])` | y = alpha x + beta y |
| `selah_stage_state([bx,] n, dt, coef, q, regs, u [, mask, mask_comp])` | u = q + sum dt coef(k) regs(k) |
| `selah_step_update([bx,] n, dt, coef, regs, q [, dq] [, mask, mask_comp])` | q = q + dt sum coef(k) regs(k) |
| `selah_shu_osher_2r_stage([bx,] dt, a, b, c, q0, r, q [, mask, mask_comp])` | q = a q0 + b q + dt c r |
| `selah_williamson_2n_stage([bx,] dt, a, b, r, dqr, q [, mask, mask_comp])` | dqr = a dqr + dt r; q = q + b dqr |

With the box the arrays are rank 5 (`regs` rank 6) with the box's lower bounds; without it, rank 1 (`regs` rank 2).
See [Kernels](./kernels).

## Example: an SSP(3,3) step with the kernels

```fortran
program selah_kernels_example
use selah
implicit none
integer(I4P), parameter :: NV=2, NI=4, NJ=4, NK=1, NGC=2, NB=3
real(R_P),    parameter :: DT=0.1_R_P
real(R_P), allocatable  :: q(:,:,:,:,:), u(:,:,:,:,:), regs(:,:,:,:,:,:)
type(butcher_tableau)   :: bt
type(box)               :: bx
integer(I4P)            :: s

bt = butcher(SCHEME_SSP_33)
bx = box_interior(LAYOUT_V_FIRST, NV, NI, NJ, NK, NGC, 2_I4P) ! 2 active blocks of a pool of 3
allocate(q(NV, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, NB), source=1._R_P)
allocate(u, mold=q) ; u = 0._R_P
allocate(regs(NV, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, NB, bt%stages), source=0._R_P)
do s=1, bt%stages
   call selah_stage_state(bx, s-1, DT, bt%a(s,:), q, regs, u) ! the state of stage s
   regs(:,:,:,:,:,s) = -u                                     ! its residual, u' = -u: the caller's job
enddo
call selah_step_update(bx, bt%stages, DT, bt%b, regs, q)       ! q = q + dt sum_s b_s F_s
print '(A,F12.10)', 'interior cell:  ', q(1,1,1,1,1)   ! 1 - h + h^2/2 - h^3/6 at h = 0.1
print '(A,F12.10)', 'ghost cell:     ', q(1,0,1,1,1)   ! untouched
print '(A,F12.10)', 'inactive block: ', q(1,1,1,1,3)   ! untouched
endprogram selah_kernels_example
```

```text
interior cell:  0.9048333333
ghost cell:     1.0000000000
inactive block: 1.0000000000
```

## Example: the coefficients

```fortran
program selah_example
use selah
implicit none
type(butcher_tableau)       :: bt
type(splitting_composition) :: sc

bt = butcher(SCHEME_SSP_54)
print '(A,I0,A,I0,A,F6.3)', 'ssp-54: ', bt%stages, ' stages, order ', bt%order, ', SSP coefficient ', &
                            bt%ssp_coefficient
print '(A,I0,A,I0)', 'ssp-104 buffers: Butcher ', registers_needed(SCHEME_SSP_104, FORM_BUTCHER), &
                     ', two-register ', registers_needed(SCHEME_SSP_104, FORM_KETCHESON_2R)
sc = splitting(SCHEME_BLANES_MOAN_S6)
print '(A,I0,A,I0)', 'blanes-moan-s6: ', sc%flows, ' flows, order ', sc%order
endprogram selah_example
```

```text
ssp-54: 5 stages, order 4, SSP coefficient  1.508
ssp-104 buffers: Butcher 11, two-register 2
blanes-moan-s6: 13 flows, order 4
```
