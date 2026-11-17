---
title: 5. Low-storage forms
---

# 5. Low-storage forms

The Butcher form of an $s$-stage scheme keeps $s+1$ field-sized registers. On a GPU, memory is the scarce resource,
and several schemes can be run with two:

| Form | Schemes | Registers | The update, at the end of stage k |
|---|---|---|---|
| `FORM_SHU_OSHER_2R` | `euler`, `ssp-22`, `ssp-33` | 2 | $q = a_k q^n + b_k q + \Delta t\, c_k R$ |
| `FORM_WILLIAMSON_2N` | `lsrk-33`, `lsrk-54` | 2 | $\Delta Q = A_k \Delta Q + \Delta t R$; $q = q + B_k \Delta Q$ |
| `FORM_KETCHESON_2R` | `ssp-104` | 2 | Ketcheson's recipe: Euler steps of $\Delta t/6$, a combination at stage 5 |

In these forms the stage state is `q` itself, advanced in place at every `end_stage`. The caller asks which array holds
the stage state with `stage_state_is_q(k)`:

<<< @/examples/snippets/low_storage-step.F90{fortran}

## Running it

<<< @/examples/output/low_storage.txt{text}

Each scheme gives the same solution in both forms, to round-off (the forms add the same terms in a different order),
with a third of the registers or less.

::: tip Which form, then?
On a GPU, take the low-storage form when the scheme has one: less memory and fewer field passes per step
([Performance](/guide/performance)). Keep the Butcher form when you need the step residual `dq` (chapter 10), which
only it can compute, or a scheme that has no other form (`ssp-54`, `rk4`).
:::

Next: [6. Splittings](./06-splittings).
