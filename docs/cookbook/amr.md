---
title: AMR
---

# AMR

[[toc]]

## Reflux with stage_weight

At a face where two sides compute different fluxes (a coarse and a fine level, two solvers coupled), conservation needs
one side corrected with the other's flux, summed over the step with the stage weights. Accumulate the difference inside
each stage window, and apply it once after `finish_step`.

<<< @/examples/snippets/cookbook_amr-reflux.F90{fortran}

<<< @/examples/output/amr-reflux.txt{text}

- The register is zeroed at the start of each step and accumulated with `stage_weight(k)`, in every Runge–Kutta form
  (the low-storage forms included: their weights are derived).
- The correction replaces block 1's flux through the shared face with block 2's: $q_N \mathrel{-}= \frac{\Delta t}{\Delta x}\sum_k b_k (F_2 - F_1)_k$.
- Accumulate from the stage state, in the window between `begin_stage` and `end_stage`, after the ghost fill: the flux
  of stage `k` is the flux of the state of stage `k`.

## Fill the ghost cells from the stage state

The ghost cells of a stage are a function of the neighbours' **stage state**. Filling them from the neighbours' $q^n$
(a stale exchange) is a classic bug: the program still runs, and still converges, but at first order.

<<< @/examples/snippets/cookbook_amr-ghosts.F90{fortran}

<<< @/examples/output/amr-ghosts.txt{text}

In a multi-block or MPI solver, the exchange happens inside the stage window, on the arrays `stage_state()` (or `q`,
when `stage_state_is_q(k)`) of every block.
