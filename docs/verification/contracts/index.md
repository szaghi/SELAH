---
title: Contracts
icon: selah
---

# Contracts

The guarantees that are not physics: identities with ADAM and between the storage forms, the stage contract a solver
programs against, values pinned between releases, the device, and the suite's own power to see defects.

| Case | Type | What it verifies |
|---|---|---|
| [C1 Equivalence with ADAM](./adam-equivalence) | contract | ADAM's coefficients and FNL kernels, bitwise; ADAM's CPU loops within 4 ε |
| [C2 Storage forms](./forms) | contract | every low-storage form computes its Butcher form's solution to round-off |
| [C3 Stage contract](./stage-contract) | contract | ghost cells, aliasing, stages on q, masks, caller registers, the error stops |
| [C4 Golden regression](./golden) | contract | 1312 values, written once, compared on every compiler |
| [C5 Device](./device) | contract | kernels and stepper on an NVIDIA A30, against the host and ADAM |
| [C6 Mutation testing](./mutation) | contract | 28 deliberate defects, each caught by the suite |
| [C7 Precision](./precision) | contract | what runs, and what is asserted, in single precision |
| [C8 Companion states](./companion) | contract | a companion and its field are the coupled Runge–Kutta step, bitwise; active extent, layouts, rebuild; particles pushed on it |
| [C9 Stage currents](./stage-currents) | contract | Gauss's law at every stage and step with the stage currents; the look-ahead state |
