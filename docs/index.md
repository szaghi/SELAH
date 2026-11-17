---
layout: home

hero:
  name: SELAH
  text: Stage-wise Evolution Library for Accelerators and HPC
  tagline: The time integrators of ADAM, carved out as stage-wise operators, on every CPU and GPU. Time-evolution operators for block-structured field solvers, in modern Fortran.
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/SELAH

features:
  - title: Verified tableaux
    details: Butcher, Shu–Osher, low-storage and splitting coefficients, with their order conditions checked in quadruple precision.
  - title: Kernels over your arrays
    details: Interior-only stage and update operations on the caller's arrays, one source for host OpenMP, OpenACC and OpenMP offload.
  - title: A stepper you drive
    details: begin_stage, end_stage, finish_step, under a documented stage contract that an AMR orchestrator can rely on.
---

::: warning Pre-release
SELAH is at milestone M2: the coefficients (L0) and the kernels (L1, verified on the host) are in; the device tests
and the stepper are next. The design is in `notes/PLAN.md` of the repository.
:::
