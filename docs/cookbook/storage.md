---
title: Storage
---

# Storage

[[toc]]

## Registers from your own pool

A solver that manages its own memory (a pool, a device arena) passes registers to the stepper instead of letting it
allocate them. Size the pool once for the largest need: `registers_needed(scheme, form)` for each scheme the solver may
use.

<<< @/examples/snippets/cookbook_storage-pool.F90{fortran}

<<< @/examples/output/storage-pool.txt{text}

The rules, checked by `initialize` (it stops with a message otherwise):
- the array is rank 6, with the field's bounds on dimensions 1–5 (`lb`, `ub`);
- the slots are numbered from 1 on dimension 6, and there are at least `registers_needed` of them.

The rules the stepper cannot check:
- **Memory space:** in a device build the pool must be device memory (FUNDAL's `dev_alloc`).
- **Ownership:** the stepper never frees it; `destroy` only forgets it.
- **No sharing between live steppers:** two steppers in use at the same time need two pools, because the registers are
  scratch between `begin_step` and `finish_step`.
- **No zeroing needed:** no form reads a register before writing it.

## The working precision

SELAH computes in one real kind per build, `R_P`: double precision by default. The `*-r4` build modes compile SELAH and
PENF with `-D_R_P_IS_R4P` and give single precision. Quadruple precision (`-D_R_P_IS_R16P`) is rejected at compile time.
A build in which SELAH and PENF disagree on the kind does not compile.

<<< @/examples/snippets/cookbook_storage-precision.F90{fortran}

<<< @/examples/output/storage-precision.txt{text}

See [Precision](/guide/precision) for what single precision does to the coefficients and the verification.
