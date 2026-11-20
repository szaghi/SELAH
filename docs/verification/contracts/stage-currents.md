---
title: C9 Stage currents
icon: selah
---

# C9 Stage currents

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C9 | Contracts | C | CPU (gfortran, nvfortran, ifx; R8P and R4P); GPU | no | `selah_test_stage_current` |

## Case and scope

`lookahead(k)` and `stage_current(k)` make the weighted sums of a field's stage currents the charge-conserving current
of the particles' displacement to each stage ([guide](/guide/particles#charge-conserving-stage-currents)), so that
Gauss's law holds at every stage and step. The case checks it on a stencil-free model, every Runge–Kutta scheme.

## Setup

| Item | Value |
|---|---|
| Model | a particle $x' = v$, $v' = -x$ (a companion), a field $D' = -J$ (one cell), the deposit $H(x^n \to X) = (\sin X - \sin x^n)/\Delta t$: a nonlinear stand-in for Esirkepov's current; the conserved charge is $D + \sin x$ |
| Stage $k$ | the particle's residual; `lookahead(k)` $\to X$; $H$; `stage_current(k)` $\to J_k$; the field's residual $-J_k$ |
| Steps | 10 of $\Delta t = 0.3$, $x_0 = 0.7$, $v_0 = 0.4$, $D_0 = 0.2$ |
| Admissibility | the catalogue; ssp-33 with $a_{32} = 0$; strang |

## Program

```bash
./exe/selah_test_stage_current
```

## Results

| Check | Result | Bound |
|---|---|---|
| $\lvert D + \sin x - (D + \sin x)^0\rvert$ in every stage state and after every step | at most 3.3e-16 (ssp-104), 0 for euler | 64 ε |
| `lookahead(k)` against the state of stage $k+1$ | bitwise, every scheme | exact |
| `lookahead(s)` against the step's result | to round-off | 8 ε |
| admissibility | every catalogue scheme yes; $a_{32} = 0$ no; strang no | exact |

Four mutants each fail it: the row $k$ in place of $k+1$, the look-ahead without the current stage's residual, the
current not divided by $a_{k+1,k}$, one term of the sum dropped.

## Findings

- The construction keeps the discrete Gauss law to round-off at every stage, in every scheme of the catalogue, with a
  nonlinear deposit: it needs nothing of the deposit but that it be a function of the end points.
