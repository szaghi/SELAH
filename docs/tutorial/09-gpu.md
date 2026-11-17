---
title: 9. On the GPU
---

# 9. On the GPU

SELAH runs on the GPU through [FUNDAL](https://github.com/szaghi/FUNDAL):
- every kernel is written once, with FUNDAL's portable directives;
- the stepper allocates its registers with `dev_alloc`.

In an OpenACC build both live in device memory. The caller does the same with its own field, and nothing crosses the
bus during a step.

<<< @/examples/snippets/gpu-field.F90{fortran}

The caller's ghost fill and residual become device loops, written with the same macros SELAH uses
(`fundal.H`: `DEVICEVAR` for OpenACC, `DEVICEPTR`/`OMPLOOP` for OpenMP; on the host they become plain loops):

<<< @/examples/snippets/gpu-kernels.F90{fortran}

The step itself is chapter 4's, unchanged, called with the device array `q_d`.

## Building it for the GPU

```bash
fobis build --mode static-nvf     # or your own nvfortran build of SELAH with -acc
nvfortran -cpp -acc -gpu=cc80 -DDEV_OAC -DCOMPILER_NVF -DFUNDAL_NO_MPI \
          -I <selah mods> -I src/third_party/FUNDAL/src/lib gpu.F90 libselah.a -o gpu
```

## Running it

The documentation is built on the host, so this is the host build's output; the device build prints the same error.
The device runs are verified separately on an NVIDIA A30 ([C5](/verification/contracts/device)).

<<< @/examples/output/gpu.txt{text}

::: warning Device pointers stay on the device
On a GPU, `q_d`, `stage_state()` and `rhs_register()` are device addresses. Read and write them only inside device
loops or through FUNDAL's copies (`dev_memcpy_to_device`, `dev_memcpy_from_device`), never element by element from the
host. Coefficient arrays passed to kernels must be device arrays too, and contiguous.
:::

::: tip What you learned
The same stepper on device memory; the caller's loops on the device; one copy back at the end.
Reference: [Installation: backends](/guide/installation#backends), [Cookbook: performance](/cookbook/performance).
:::

Next: [10. Production](./10-production).
