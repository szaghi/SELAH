---
title: Properties
---

# Properties

The properties a solver relies on beyond accuracy: strong stability, the size of the SSP step, the long-time behaviour
of a Hamiltonian, conservation through a flux register.

| Case | Type | What it verifies |
|---|---|---|
| [P1 Strong stability (TVD)](./ssp-tvd) | physics | every SSP scheme keeps a TVD operator TVD up to $\nu = C$ |
| [P2 SSP coefficients](./ssp-coefficients) | physics | the coefficients computed from the Shu–Osher forms match the published ones |
| [P3 Splitting energy](./splitting-energy) | physics | the energy error of a splitting stays bounded over 1000 periods; a Runge–Kutta scheme's drifts |
| [P4 Conservation](./conservation) | physics | the stage fluxes weighted by `stage_weight` telescope to the mass change, in every form |
