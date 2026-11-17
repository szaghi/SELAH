---
title: C5 Device
---

# C5 Device

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C5 | Contracts | C | GPU (NVIDIA A30, nvfortran 26.9 OpenACC, `tests-nvf-oac --varset iac`) | no | `selah_test_device_kernels`, `selah_test_device_adam`, `selah_test_stepper`, `selah_test_stepper_options`, `selah_vv_stability` |

## Case and scope

On a GPU the kernels run as device loops on FUNDAL device memory, and the stepper's registers live there too. The
case runs the device-capable programs on an NVIDIA A30 (CNR-IAC medooza, one GPU of a `dev-gn` node). In a device
build `dev_init(require_device=.true.)` forbids a silent fallback to the host.

## Setup

| Program | On the device |
|---|---|
| `selah_test_device_kernels` (57 checks) | every box kernel, both layouts, masked and unmasked, the stage state with no register, the step update without dq, the 2N stage with $A = 0$, a component range, an explicit mask component, every flat kernel: each cell of the box within 16 ε of the host reference (the device may contract to FMA), every other cell bitwise untouched |
| `selah_test_device_adam` (18 checks) | SELAH's kernels against ADAM's FNL kernels, both on the device: bitwise |
| `selah_test_stepper` (248 checks) | the whole stepper with the caller's ghost fill and residual as device loops: orders, contract, telescoping, stages on q |
| `selah_test_stepper_options` (23 checks) | `stage_time`, NaN-filled caller registers (device memory), `dq`, masks |
| `selah_vv_stability`, the L0 programs | the coefficient checks, built by nvfortran with OpenACC |

The programs that pass host arrays to the kernels or the stepper are excluded from the device build. Their device
twins are listed above.

## Program

```bash
fobis build --mode tests-nvf-oac --varset iac && bash scripts/run_tests.sh   # on a GPU node
sbatch --export=ALL,SELAH_SRC=$PWD scripts/medooza/selah-tests.slurm           # the medooza job script
```

## Results

| Job (2026-11-17) | Tree | Outcome |
|---|---|---|
| 1368142 | the audit fixes | 14/14 |
| 1368177 | the verification programs | 15/15 |
| the final job of the night | the documented tree | see the morning brief of the night log |

## Findings

- Every device check passes. The kernels agree with ADAM's FNL kernels bitwise on the device too.
- The documentation's examples run on the host ([Tutorial 9](/tutorial/09-gpu)). The device behaviour is verified
  here.
