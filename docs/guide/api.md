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

## Example

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
