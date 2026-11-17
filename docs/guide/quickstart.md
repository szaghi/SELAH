---
title: Quick start
---

# Quick start

One step of a solver with SELAH: the stepper builds the stage states and the update; you fill the ghost cells and
write the residual. This is the loop of [Tutorial 4](/tutorial/04-stepper), on a periodic 1-D advection.

<<< @/examples/snippets/stepper-setup.F90{fortran}

<<< @/examples/snippets/stepper-step.F90{fortran}

- `box_interior` describes the field: layout, components, cells per block, ghost cells, active blocks.
- `initialize` takes the scheme, the storage form and the field's bounds, and allocates the registers.
- Between `begin_stage(k)` and `end_stage(k)` the stage is yours: the ghost cells of the stage state, and its residual
  in `rhs_register()`.

The whole program, and its output (fourth order, as SSP(5,4) promises):

<<< @/examples/output/stepper.txt{text}

## Build and test

```bash
fobis fetch && fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis build --mode tests-gnu-r4 && for t in exe-r4/selah_*; do $t > /dev/null || echo "FAIL $t"; done
bash scripts/run_vv.sh          # the verification report: data, figures, tables
```

Next: the [Concepts](./concepts), the [Tutorial](/tutorial/), or the [Cookbook](/cookbook/).
