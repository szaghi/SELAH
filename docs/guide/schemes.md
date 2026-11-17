---
title: Schemes
icon: tableaux
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

<SelahIcon name="ssp" :size="40" label="strong stability preserving" float />

The SSP coefficient is computed from the Shu–Osher form (`ssp_coefficient`), never typed in. It is the minimum of
$\alpha_{ik}/\beta_{ik}$ over the coefficients SELAH stores: another Shu–Osher representation of the same tableau can
give a larger value. For the catalogue the stored forms are the optimal ones, so the value is the scheme's SSP
coefficient ([P2](/verification/properties/ssp-coefficients)). For a tableau built by hand, 0 means "not computed" as
well as "not SSP".

The names are the constants `SCHEME_EULER`, `SCHEME_SSP_22`, …, `SCHEME_LSRK_54`; `RK_SCHEMES` lists them all. The
entries of `RK_SCHEMES` are fixed-length strings padded with blanks: `trim` them before passing them on. An unknown
name gives a tableau of 0 stages.

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

<SelahIcon name="low-storage" :size="40" label="low storage" float />

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

### Every tableau

Printed from the library itself (`docs/examples/src/tableaux.F90`, regenerated with the documentation):

<<< @/examples/output/tableaux.txt{text}

`butcher_tableau%b_hat`, the embedded weights of a pair, is reserved: no scheme of the catalogue has one, and it stays
unallocated.

## Splitting schemes

<SelahIcon name="splitting" :size="40" label="splitting" float />

For $u' = F_1(u) + F_2(u)$, a step is the product of the sub-flows of the two partitions, each over a fraction of
$\Delta t$. Every scheme is a palindrome, stored as its half sequence, so its order is even:

| Name | Flows | Order | Source |
|------|-------|-------|--------|
| `strang`          | 3  | 2 | Strang (1968) |
| `yoshida-4`       | 7  | 4 | Yoshida (1990), the triple jump |
| `blanes-moan-s6`  | 13 | 4 | Blanes & Moan (2002), S6 |
| `blanes-moan-s10` | 21 | 6 | Blanes & Moan (2002), S10 |

The flows apply left to right, flow 1 first, and `splitting(name, first=2)` starts with partition 2. Through the
stepper a flow is one forward-Euler step of its partition, $y_p \mathrel{+}= c_k\,\Delta t\,R_p$. That is the exact
flow of partition $p$ only when $R_p$ does not depend on partition $p$ itself, as with the positions and momenta of a
separable Hamiltonian. Otherwise the composition is first order, whatever its declared order
([A4](/verification/accuracy/splittings), [Tutorial 6](/tutorial/06-splittings)).

## Verification

| What | Case |
|---|---|
| the order conditions, in quadruple precision | [A1](/verification/accuracy/order-conditions) |
| the observed orders, scalar problems, through the stepper, splittings | [A2](/verification/accuracy/scalar-orders), [A3](/verification/accuracy/stepper-orders), [A4](/verification/accuracy/splittings) |
| the stability regions | [A5](/verification/accuracy/stability) |
| the SSP coefficients and strong stability | [P2](/verification/properties/ssp-coefficients), [P1](/verification/properties/ssp-tvd) |
| ADAM's coefficients, bitwise | [C1](/verification/contracts/adam-equivalence) |
