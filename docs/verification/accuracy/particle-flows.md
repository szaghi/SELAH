---
title: A7 Particle sub-flows
icon: order-conditions
---

# A7 Particle sub-flows

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A7 | Accuracy | P | CPU (gfortran, nvfortran, ifx; R8P, the oscillator in R4P); GPU | no | `selah_test_particle_flows` |

## Case and scope

A particle array `p(6, np)` stepped by the stepper's splitting form, the velocities as partition 1 (the kick) and the
positions as partition 2 (the drift): `strang` is kick–drift–kick leapfrog with SELAH's fractions
([guide](/guide/composition#particles-kick-and-drift)). The case checks it against KDK's own discrete solution, the
orders and long-time energy of every splitting on a Kepler orbit, and a gas–sink system coupled through
[A6](./additive-splitting)'s additive splitting.

## Setup

| Item | Value |
|---|---|
| Oscillator | $x'' = -x$, $x_0 = 1$, $v_0 = 0$, $\omega h$ = 0.1 and 0.5, 100 periods; KDK's map has trace $2 - (\omega h)^2$ and determinant 1, so $x_n = \cos(n\theta)$, $\cos\theta = 1 - (\omega h)^2/2$ |
| Kepler | $GM = 1$, $e = 0.5$, from the perihelion ($r = 0.5$, $v = \sqrt 3$), period $2\pi$; orders from one orbit at 4 halvings (from 128 steps, 64 for order 6); energy over 50 orbits of 256 steps, its largest error over every step of orbits 1–10 and 41–50 |
| Gas and sink | $m' = -\kappa x^2 m$ (the gas, rk4 through the stepper, the sink frozen), $x'' = -(1 + m) x$ (the sink, one KDK step over the substep, the gas frozen), $\kappa = 2$, $t = 2$; reference rk4 on the whole system, 40 000 steps |
| Bounds | $\lvert x_n - \cos n\theta \rvert \le 10^4\,\varepsilon$; phase error per period $\pi(\omega h)^2/12$ to 5 %; orders within 0.2 (Kepler), 0.15 (gas–sink); energy error of orbits 41–50 at most 1.2 times that of orbits 1–10 |

## Program

```bash
./exe/selah_test_particle_flows
```

## Results

| Check | Value |
|---|---|
| $\max_n \lvert x_n - \cos n\theta \rvert$, $\omega h$ = 0.1 / 0.5 | 2.9e-13 / 7.2e-14 |
| phase error per period, $\omega h$ = 0.1 / 0.5 (against $\pi(\omega h)^2/12$) | 2.621e-3 (2.618e-3) / 6.736e-2 (6.545e-2) |
| Kepler orders, strang / yoshida-4 / blanes-moan-s6 / blanes-moan-s10 | 2.00 / 4.00 / 4.00 / 5.99 |
| Kepler energy error, orbits 1–10 and 41–50 (equal) | 8.2e-4 / 3.8e-6 / 2.9e-8 / 2.2e-11 |
| gas and sink, lie / strang orders | 0.99 / 2.00 |

## Findings

- **SELAH's strang on a particle array is KDK leapfrog**, to round-off: its discrete solution, not only its order.
- **The energy error is bounded with a fixed step.** Sampled once per orbit, the energy error seemed to grow 4-fold
  over 50 orbits; that sample follows the precessing perihelion. The largest error over every step of an orbit is the
  same in orbit 1 and orbit 50 for every scheme.
- **Gas and sinks compose at the expected orders** through the additive splitting with the sink's KDK step as the
  flow.
