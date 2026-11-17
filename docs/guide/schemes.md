---
title: Schemes
---

# Schemes

SELAH's L0 layer holds the coefficients: every scheme in the form its source publishes, the other forms derived from
it, and every value verified by a test. Nothing in this layer touches a field. The kernels (L1) and the stepper (L2)
read these coefficients.

## Runge–Kutta schemes

| Name | Stages | Order | SSP coefficient | Primary form | Other forms | Source |
|------|--------|-------|-----------------|--------------|-------------|--------|
| `euler`   | 1  | 1 | 1           | Butcher       | Shu–Osher, two-register | forward Euler |
| `ssp-22`  | 2  | 2 | 1           | Butcher       | Shu–Osher, two-register | Shu & Osher (1988) |
| `ssp-33`  | 3  | 3 | 1           | Butcher       | Shu–Osher, two-register | Shu & Osher (1988) |
| `ssp-54`  | 5  | 4 | 1.508180049 | Shu–Osher     | Butcher | Spiteri & Ruuth (2002); 15 digits of Isherwood, Grant & Gottlieb (2018) |
| `ssp-104` | 10 | 4 | 6           | Shu–Osher     | Butcher | Ketcheson (2008) |
| `rk4`     | 4  | 4 | —           | Butcher       | — | Kutta (1901) |
| `lsrk-33` | 3  | 3 | —           | Williamson 2N | Butcher | Williamson (1980) |
| `lsrk-54` | 5  | 4 | —           | Williamson 2N | Butcher | Carpenter & Kennedy (1994), NASA TM-109112 |

The SSP coefficient is computed from the Shu–Osher form (`ssp_coefficient`), never typed in.

### The coefficient forms

For an autonomous system $u' = F(u)$, with $Y_i$ the stages:

- **Butcher**: $Y_i = u^n + \Delta t \sum_j a_{ij} F(Y_j)$, $u^{n+1} = u^n + \Delta t \sum_j b_j F(Y_j)$, $c = A\mathbf 1$.
  The form the order conditions and the reflux weights $b_j$ are stated in.
- **Shu–Osher**: $Y_1 = u^n$, $Y_i = \sum_{k<i} \left(\alpha_{ik} Y_k + \Delta t\, \beta_{ik} F(Y_k)\right)$,
  $u^{n+1} = Y_{s+1}$. With $\alpha_{ik}, \beta_{ik} \ge 0$ each stage is a convex combination of forward-Euler
  steps, and the SSP coefficient is $\min \alpha_{ik}/\beta_{ik}$.
- **Two-register Shu–Osher**: $Y_{i+1} = a_i Y_1 + b_i Y_i + c_i \Delta t F(Y_i)$: only $u^n$ and the current
  stage are kept. These are ADAM's `runge-kutta-1/2/3` (`ark`, `brk`, `crk`).
- **Williamson 2N**: $\Delta Q_i = A_i \Delta Q_{i-1} + \Delta t F(Y_i)$, $Y_{i+1} = Y_i + B_i \Delta Q_i$.

`runge-kutta-2` and `runge-kutta-3` of ADAM are `ssp-22` and `ssp-33` in two-register form: the same methods, with
less storage and no access to the per-stage residuals.

### Storage

`registers_needed(scheme, form)` counts the field-sized buffers a step needs besides the solution, the residual
buffer included:

| Form | Buffers | Schemes |
|------|---------|---------|
| `FORM_BUTCHER`       | s + 1 (stage state + s residuals) | every RK scheme |
| `FORM_SHU_OSHER_2R`  | 2 | `euler`, `ssp-22`, `ssp-33` |
| `FORM_WILLIAMSON_2N` | 2 | `lsrk-33`, `lsrk-54` |
| `FORM_KETCHESON_2R`  | 2 | `ssp-104` |
| `FORM_SPLITTING`     | 1 | every splitting scheme |

`ssp-104` is the reason the low-storage forms matter. In Butcher form it needs 11 buffers, but 2 in Ketcheson's. Its
effective SSP coefficient per stage is 0.6, against 0.30 for `ssp-54`.

## Splitting schemes

For $u' = F_1(u) + F_2(u)$, a step is the product of the sub-flows of the two partitions, each over a fraction of
$\Delta t$. Every scheme is a palindrome, stored as its half sequence, so its order is even:

| Name | Flows | Order | Source |
|------|-------|-------|--------|
| `strang`          | 3  | 2 | Strang (1968) |
| `yoshida-4`       | 7  | 4 | Yoshida (1990), the triple jump |
| `blanes-moan-s6`  | 13 | 4 | Blanes & Moan (2002), S6 |
| `blanes-moan-s10` | 21 | 6 | Blanes & Moan (2002), S10 |

## Verification

All numbers come from the test programs (gfortran 16.2.0, 2026-11-17); the tests also pass with nvfortran 26.9
and ifx 2026.2.

**Order conditions.** The eight conditions up to order 4 are evaluated in quadruple precision on the stored
coefficients (R8P with nvfortran, which has no quad kind), so the residual is the tableau's. The largest residual up
to the declared order:

| `euler` | `ssp-22` | `ssp-33` | `ssp-54` | `ssp-104` | `rk4` | `lsrk-33` | `lsrk-54` |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 5.6e-17 | 5.1e-16 | 6.0e-17 | 5.6e-17 | 1.2e-17 | 8.7e-17 |

For the schemes of order below 4, some condition of the next order fails: the declared order is exact.

**Observed order.** The orders measured on $y' = -y$ and $y' = -y^2$, between 64 and 128 steps on $[0, 1]$, all lie
within 0.02 of the declared order. For the splittings, the system is $q' = p$, $p' = -V'(q)$ with exact sub-flows,
measured with either partition first, on the harmonic oscillator and on the pendulum:

| Scheme | Observed order, finest pair |
|--------|-----------------------------|
| `strang`          | 2.000 |
| `yoshida-4`       | 3.999–4.000 |
| `blanes-moan-s6`  | 4.000 |
| `blanes-moan-s10` | 6.002–6.006 |

**Published values reproduced.** The SSP coefficients (1, 1, 1, 1.508, 6), the weights of `ssp-104` (all 1/10) and
its abscissas, and the abscissas of Carpenter & Kennedy.

**Equivalence with ADAM.** The `euler`, `ssp-22`, `ssp-33` and `ssp-54` tableaux are bitwise ADAM's `ssp_tableau`,
and the two-register forms are bitwise ADAM's `ark`, `brk`, `crk`. The test uses ADAM's code verbatim as its
oracle (`src/tests/selah_adam_reference.F90`), so wiring SELAH's coefficients into ADAM changes no golden output.
