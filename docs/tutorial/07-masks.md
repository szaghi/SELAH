---
title: 7. Masks
---

# 7. Masks

Immersed-boundary solvers keep solid cells in the field, and the integrator must leave them alone. Every kernel, and
every stepper call, takes an optional **mask**: an array with the field's bounds and one or more components. A cell is
active where the mask is negative, as with ADAM's level-set distance $\phi$, which is negative in the fluid.

<<< @/examples/snippets/masks-mask.F90{fortran}

<<< @/examples/snippets/masks-step.F90{fortran}

- Pass the same mask to all four calls. `begin_step` copies $q^n$ only where it is active, and `finish_step` updates
  only there.
- With several mask components, `mask_comp=` picks one; the default is the last (ADAM's "all solids" summary).
- The residual can be computed everywhere: whatever it holds in solid cells, the masked update ignores.

## Running it

<<< @/examples/snippets/masks-result.F90{fortran}

<<< @/examples/output/masks.txt{text}

Block 3 kept every value bit for bit.

::: tip What you learned
Masks: the sign convention, the component, and that they apply to every stepper call.
Reference: [Concepts: masks](/guide/concepts#masks).
:::

Next: [8. The flux register](./08-flux-register).
