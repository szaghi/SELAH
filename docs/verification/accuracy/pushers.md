---
title: A8 Particle pushers
icon: order-conditions
---

# A8 Particle pushers

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A8 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P and R4P, the volume values in R8P); GPU | no | `selah_test_pushers` |

## Case and scope

`selah_push` advances particles by one staggered leapfrog step with one of four pushers ([guide](/guide/particles)).
The case checks each against its defining implicit equation and against the properties the literature gives it (Boris
1970; Vay 2008; Higuera & Cary 2017; Qin et al. 2013), with the reference values computed at 40 digits by the
research of `notes/companions-pushers-pml.md`.

## Setup

| Item | Value |
|---|---|
| Units | $q/m = 1$, $c = 1$ (relativistic), one particle per call through device arrays |
| Boris | pure $B = \omega\hat z$, $\omega\Delta t$ = 0.1, 1, 3: the turning angle of one step; $\lvert v\rvert$ over 1000 steps in a generic $B$; 16 steps from $v_E$ ($E = 0.5\hat y$, $B = 2\hat z$); both layouts, 7 steps in generic fields |
| Defining equations | $E = (0.3, -0.2, 0.15)$, $B = (0.1, 0.4, 1.3)$, $u = (0.7, 1.9, -0.4)$, $\Delta t = 0.8$ |
| Force-free | $\gamma_0 = 10$ along $\hat x$, $B = \hat z$, $E = -v\times B$, $qB\Delta t/m = 0.1$, one step |
| Relativistic drift | $E = 0.5\hat y$, $B = \hat z$, $qB\Delta t/m = 1$, one step from $u_E = \gamma_E v_E$ |
| Volume | the momentum map's Jacobian by centred differences ($h = 10^{-5}$), the fields and $u$ of the defining equations, $q\Delta t/m = 1$ |
| Non-relativistic limit | $c = 10^6$: each relativistic pusher against Boris |

## Program

```bash
./exe/selah_test_pushers
```

## Results

| Check | Result | Bound |
|---|---|---|
| Boris's angle at $\omega\Delta t$ = 0.1, 1, 3 | $2\arctan(\omega\Delta t/2)$: 0.0999168, 0.9272952, 1.9655874 | 8 ε |
| Boris's speed, 1000 steps; $v_E$, 16 steps; layouts | kept; kept; bitwise | $10^3$ ε; 8 ε; exact |
| defining equations, relative residual | 2.2e-16, 5.4e-17, 2.2e-16, 2.2e-16 (Boris, rel. Boris, Vay, HC) | 64 ε |
| force-free, $\lvert u^1 - u^0\rvert/\lvert u^0\rvert$ | rel. Boris 1.24e-7 (reference 1.237e-7); Vay 0; HC 0 | 5 % / 64 ε |
| relativistic drift from $u_E$ | rel. Boris 1.8e-2; Vay 1.9e-16; HC 2.2e-16 | moved / 64 ε |
| $\det - 1$ | Boris 1.5e-11, rel. Boris 6.6e-12, HC 3.9e-12 (the differences' floor); Vay −1.220e-4 (reference −1.22e-4) | $10^{-8}$; $\pm 2\cdot10^{-6}$ |
| $c = 10^6$ | every relativistic pusher = Boris | $10^{-10}$ |

## Findings

- The Fortran kernels reproduce the research's 40-digit reference values (the force-free error of relativistic
  Boris, the volume defect of Vay) to the digits printed: two independent implementations agree.
- Higuera–Cary is the only one of the four with all three properties (volume, force-free, drift) in the
  relativistic regime.
