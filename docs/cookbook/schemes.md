---
title: Schemes
---

# Schemes

[[toc]]

## Choose a scheme and its storage form

Pick the order you need, then the form with the fewest registers. Use the Butcher form only when you need what only it
gives (the step residual `dq`), or when the scheme has no other form.

<<< @/examples/snippets/cookbook_schemes-choose.F90{fortran}

<<< @/examples/output/schemes-choose.txt{text}

`ssp-54` and `rk4` have only the Butcher form. Every other scheme runs in two registers.

## The largest strong-stability-preserving time step

If your spatial operator is TVD (or positivity-preserving, or bounded in any convex functional) under a forward-Euler
step up to $\Delta t_{FE}$, an SSP scheme with coefficient $C$ keeps that property up to $\Delta t = C\,\Delta t_{FE}$
(Gottlieb, Shu & Tadmor 2001). The cost per unit time is the number of stages divided by $C$.

<<< @/examples/snippets/cookbook_schemes-ssp-dt.F90{fortran}

<<< @/examples/output/schemes-ssp-dt.txt{text}

- `C / stages` is the effective coefficient: the fraction of a forward-Euler step gained per residual evaluation.
- `ssp-104` (Ketcheson's SSPRK(10,4)) has the best: 0.6. It needs half the residuals of `ssp-54` per unit time, in two
  registers instead of six.

The bound is verified on a TVD advection in [P1](/verification/properties/ssp-tvd).

::: warning Not SSP is not unstable
`rk4` and the `lsrk-*` schemes have no SSP coefficient: they are not guaranteed to preserve TVD-like properties at any
step. They are still stable for smooth problems up to their linear stability limit ([A5](/verification/accuracy/stability)).
:::
