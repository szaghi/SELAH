---
title: About
---

# About SELAH

SELAH (**S**tage-wise **E**volution **L**ibrary for **A**ccelerators and **H**PC) provides time-evolution operators
for block-structured field solvers on CPUs and GPUs. It is carved out of the Runge–Kutta machinery of
[ADAM](https://github.com/szaghi/adam) and is the kernel-level successor of
[FOODIE](https://github.com/szaghi/FOODIE). FOODIE's operator-overloaded, polymorphic integrands are elegant, but
they cost too much where the cells are swept.

*Selah* is the mark that appears 71 times in the Psalms and 3 times in Habakkuk, usually read as a pause: the music
stops and the singer breathes. SELAH's stage boundaries are such pauses. The caller's orchestrator acts there (ghost
exchange, AMR seam fill, reflux) before the integrator goes on.

## Status

Milestone **M0**: the build, test and documentation scaffolding. The design (layers L0–L3, principles, milestones,
risks) is in `notes/PLAN.md` of the repository.

## Design principles

1. No polymorphism, procedure pointers or callbacks in or around a kernel.
2. The caller owns the time loop and the residual evaluation; the library never calls back.
3. Kernels act on plain contiguous arrays with explicit lower bounds.
4. Box kernels update the interior only; the integrator never writes a ghost cell.
5. One working precision per build (PENF's `R_P`: R8P by default, R4P with `-D_R_P_IS_R4P`).
6. One kernel source for host OpenMP, OpenACC and OpenMP offload.
7. The floating-point summation order of every kernel is specified, so golden tests can be bitwise.
