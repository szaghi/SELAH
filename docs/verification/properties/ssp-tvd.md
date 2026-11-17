---
title: P1 Strong stability (TVD)
---

# P1 Strong stability (TVD)

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P1 | Properties | P | CPU (gfortran, nvfortran, ifx; R8P and R4P) | yes, and a cast | `selah_vv_ssp_tvd` |

## Case and scope

Strong stability is the reason SSP schemes exist. If a spatial operator keeps a convex functional bounded under a
forward-Euler step up to $\Delta t_{FE}$ (total variation, positivity, a maximum principle), an SSP scheme with
coefficient $C$ keeps it bounded up to $C\,\Delta t_{FE}$ (Gottlieb, Shu & Tadmor, 2001). The case checks that SELAH's
schemes, run through the stepper, deliver that property, and measures where it actually breaks.

What the case catches:
- **a scheme or coefficient that is not SSP as claimed**: the total variation grows below $\nu = C$;
- **a wrong SSP coefficient**: the bound is checked at the value `butcher(scheme)` reports.

## Setup

| Item | Value |
|---|---|
| Operator | first-order upwind for $u_t + u_x = 0$, periodic, 100 cells: forward Euler is TVD for $\nu \le 1$ |
| Data | a square wave ($u = 1$ on $(1/4, 1/2]$, 0 elsewhere) |
| Integrator | every scheme through the stepper, Butcher form with stages on q, 60 steps per CFL |
| Sweep | $\nu = 0.05\,C, 0.10\,C, \ldots, 2\,C$ (in units of 1 for the schemes that are not SSP) |
| Bound | for every SSP scheme and every $\nu \le C$, no step increases $TV$ by more than $1000\,\varepsilon\,TV_0$ |
| Reported | the first swept $\nu$ at which the total variation grows |

## Program

```bash
./exe/selah_vv_ssp_tvd
SELAH_VV_DIR=exe/vv-data ./exe/selah_vv_ssp_tvd     # ssp-tvd-sweep.dat, ssp-tvd-summary.dat, ssp-tvd-frames-*.dat
```

## Results

![TVD sweep](/verification/ssp-tvd.png)

*The largest one-step growth of the total variation over 60 steps, against the CFL. Left: the SSP schemes, with the
CFL in units of their coefficient; the dashed line is $\nu = C$. Right: the schemes without an SSP coefficient.*

<!--@include: ../_tables/ssp-tvd.md-->

<div class="selah-casts">
<figure class="selah-cast selah-cast-wide" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-light.png" aria-label="SSP(5,4) at the SSP limit and above its threshold"><source src="/verification/casts/ssp-tvd-light.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-dark.png" aria-label="SSP(5,4) at the SSP limit and above its threshold"><source src="/verification/casts/ssp-tvd-dark.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-dark.mp4" type="video/mp4"></video>
<figcaption>SSP(5,4), 60 steps. At ν = C (left) the total variation (bottom) only falls. At 1.15 times the measured threshold (right) it overshoots in the first steps, the new extrema of a non-TVD step.</figcaption>
</figure>
</div>

## Findings

- Every SSP scheme is TVD up to its coefficient, as the theorem promises.
- On this *linear* operator euler, ssp-22, ssp-33 and ssp-104 lose the property at the first swept CFL above $C$,
  1.05 C (the sweep's resolution): their bound is sharp here. ssp-54 holds up to 1.25 C. $C$ is the bound valid for
  *every* operator that is TVD under forward Euler; a particular operator can tolerate more. **Do not size a step on a
  measured threshold.**
- rk4 and lsrk-33 lose the TVD property just above $\nu = 1$, lsrk-54 at 1.7: without an SSP coefficient there is
  no guarantee at any $\nu$ for a nonlinear operator.
