---
title: P3 Splitting energy
---

# P3 Splitting energy

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P3 | Properties | P | CPU (gfortran, nvfortran, ifx; R8P) | yes, and a cast | `selah_vv_splitting_energy` |

## Case and scope

The splittings of SELAH's catalogue are symmetric compositions of exact flows of a separable Hamiltonian: they are
symplectic. A symplectic method conserves a modified energy close to the true one, so its energy error oscillates
with a bounded amplitude forever. A Runge–Kutta scheme is not symplectic. On the oscillator its amplification factor
has modulus below 1, and the energy decays at every step. The case checks both behaviours, which decide what to use
for long-time Hamiltonian dynamics.

## Setup

| Item | Value |
|---|---|
| Problem | $q' = p$, $p' = -q$, $H = (q^2 + p^2)/2$, from $(1, 0)$ |
| Step | $\Delta t = 0.1$, 62 832 steps (1000 periods) |
| Schemes | strang, yoshida-4, blanes-moan-s6 (splittings, partitions $q$ and $p$); rk4, ssp-54 (Butcher form, stages on q) |
| Bound | over every step, $\max|H - H_0|$ in the last tenth of the run versus the first tenth: at most 2 for the splittings, at least 5 for the Runge–Kutta schemes |
| Precision | R8P builds (the high-order splittings' errors, about $10^{-8}$, are below R4P round-off over 62 832 steps) |

## Program

```bash
./exe/selah_vv_splitting_energy
SELAH_VV_DIR=exe/vv-data ./exe/selah_vv_splitting_energy   # splitting-energy-*.dat, splitting-phase-*.dat
```

## Results

![Energy error](/verification/splitting-energy.png)

*The largest energy error in each of 400 windows of the run (every step counted).*

<!--@include: ../_tables/splitting-energy.md-->

<div class="selah-casts">
<figure class="selah-cast" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-light.png" aria-label="Phase portrait, Strang against rk4"><source src="/verification/casts/splitting-phase-light.webm" type="video/webm"><source src="/verification/casts/splitting-phase-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-dark.png" aria-label="Phase portrait, Strang against rk4"><source src="/verification/casts/splitting-phase-dark.webm" type="video/webm"><source src="/verification/casts/splitting-phase-dark.mp4" type="video/mp4"></video>
<figcaption>Δt = 0.6, 400 steps: the Strang splitting stays on a closed curve near the exact circle; rk4 spirals inwards.</figcaption>
</figure>
</div>

## Findings

- The splittings' energy error is flat: the same amplitude in the last tenth of the run as in the first, at
  $1.3\cdot10^{-3}$ (Strang), $3.8\cdot10^{-6}$ (Yoshida) and $1.4\cdot10^{-8}$ (Blanes–Moan S6).
- rk4 and ssp-54 drift linearly: ten times the error after ten times the time. Over long times a fourth-order
  Runge–Kutta scheme ends up less accurate in energy than the second-order Strang splitting would be at a finer step.
- The property needs exact sub-flows. With a sub-flow that is only first order ([A4](../accuracy/splittings)), the
  composition is neither of its order nor symplectic.
