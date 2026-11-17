---
title: 6. Splittings
---

# 6. Splittings

A splitting scheme advances a system $y' = f_1(y) + f_2(y)$ by composing *flows*, each moving one **partition** of the
components by a fraction of the step. Hamiltonian systems split into positions and momenta, and the compositions of
SELAH's catalogue (Strang, Yoshida, Blanes–Moan) are symmetric and symplectic.

The stepper drives a splitting like a Runge–Kutta scheme, one flow per stage. You declare the partitions as component
ranges, and `stage_partition(k)` tells you which partition flow `k` advances, so you compute only that half of the
residual.

<<< @/examples/snippets/splittings-setup.F90{fortran}

<<< @/examples/snippets/splittings-step.F90{fortran}

- `partitions(:,p)` holds `(v1, v2)`, the components of partition `p`: here `q` alone and `p` alone.
- `end_stage` advances the partition's components by $c_k\,\Delta t\,R$: a forward-Euler sub-flow.

## Running it

<<< @/examples/output/splittings.txt{text}

::: warning When a splitting reaches its order
The sub-flow $y_p \mathrel{+}= c_k \Delta t\, R_p$ is the exact flow of partition $p$ only when $R_p$ does not depend on
partition $p$ itself: $\dot q = p$ moves $q$ with $p$ frozen, and $\dot p = -q$ moves $p$ with $q$ frozen. This holds
for separable Hamiltonians, kinetic and potential. For other splittings the forward-Euler sub-flow is first order and
caps the composition at first order.
:::

::: tip What you learned
Partitions, flows, `stage_partition`, and when the declared order holds.
Reference: [Schemes: splittings](/guide/schemes#splittings), [A4](/verification/accuracy/splittings).
:::

Next: [7. Masks](./07-masks).
