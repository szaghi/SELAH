---
title: 8. The flux register
---

# 8. The flux register

On an AMR grid, a coarse cell next to a fine level is updated with the coarse flux through the shared face. The fine
level computes its own, more accurate fluxes. To keep the scheme conservative, the solver accumulates the difference
over the step in a **flux register** and corrects the coarse cell afterwards (refluxing).

With a Runge–Kutta step, "the flux over the step" is a weighted sum of the stage fluxes, and the weights are the
scheme's $b_k$:

$$
q^{n+1} - q^n = -\frac{\Delta t}{\Delta x} \sum_k b_k \left(F_{\text{right}} - F_{\text{left}}\right)_k .
$$

`stage_weight(k)` returns $b_k$ in every Runge–Kutta form, including the low-storage ones, whose coefficients are not
$b_k$ but imply them.

<<< @/examples/snippets/flux_register-register.F90{fortran}

The program advects a profile with inflow on the left and outflow on the right. For each stage it adds the stage's net
boundary flux, once with `stage_weight(k)` and once with uniform weights $1/s$, then compares both sums with the
change of the mass over the step.

## Running it

<<< @/examples/output/flux_register.txt{text}

With `stage_weight` the identity holds to round-off. With uniform weights it misses by the size of the time error,
which is how a refluxed AMR solver silently stops conserving.

::: tip What you learned
`stage_weight(k)` is the weight of stage `k` in the step's fluxes, in every form.
Reference: [Stepper](/guide/stepper), [P4 Conservation](/verification/properties/conservation),
[Cookbook: AMR](/cookbook/amr).
:::

Next: [9. On the GPU](./09-gpu).
