---
title: 1. A step by hand
---

# 1. A step by hand

Before any field, any box or any stepper: what does SELAH compute? An explicit Runge–Kutta step of
$y' = f(y)$ with $s$ stages, written in Butcher's form, is two operations repeated:

$$
U_s = y^n + \Delta t \sum_{j<s} a_{sj} F_j, \qquad F_s = f(U_s), \qquad
y^{n+1} = y^n + \Delta t \sum_{s} b_s F_s .
$$

SELAH provides the first and the third as kernels. The second, the residual, is always yours. This chapter writes a
whole step of SSP(5,4) with them, on three decoupled equations $y_i' = -\lambda_i y_i$.

<<< @/examples/snippets/by_hand.F90{fortran}

- `butcher(SCHEME_SSP_54)` returns the scheme's coefficients as a `butcher_tableau`: `stages`, `order`, the matrix `a`,
  the weights `b` and the abscissae `c`.
- `selah_stage_state(n, dt, coef, q, regs, u)` writes `u = q + dt Σ_{k≤n} coef(k) regs(:,k)`. For stage `s`, `n = s-1`
  and `coef` is row `s` of A. With `n = 0` it is a copy.
- The loop body `regs(:,s) = -LAMBDA * u` is the caller's residual: SELAH never calls back into your code.
- `selah_step_update(n, dt, coef, regs, y)` writes `y = y + dt Σ_k coef(k) regs(:,k)`.

These are the **flat** forms of the kernels: they work on rank-1 arrays of any length (particles, a reduced model, a
whole field flattened). Chapter 3 introduces their **box** forms for block-structured fields.

## Running it

<<< @/examples/output/by_hand.txt{text}

Halving the step divides the error by 16: fourth order, as SSP(5,4) promises.

::: warning Contiguous coefficients on a device
The program copies row `s` of A into `a_row` before the call. On the host a strided section `bt%a(s,:)` would work too,
but a device kernel receives its arrays by address (`deviceptr`), so on a GPU the coefficients must be contiguous. The
stepper keeps the transpose of A for this reason.
:::

::: tip What you learned
A Runge–Kutta step is stage states, residuals and an update; SELAH computes the first and the last, you compute the
residual. Reference: [Kernels](/guide/kernels), [Schemes](/guide/schemes).
:::

Next: [2. The catalogue](./02-catalogue).
