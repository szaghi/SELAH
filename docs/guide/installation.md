---
title: Installation
icon: cpu-gpu
---

# Installation

## Requirements

- **Compiler:** Fortran 2018. The reference is gfortran 16.2; nvfortran 26.9 (also for OpenACC) and ifx 2026.2 are
  tested too.
- **Build tool:** [FoBiS](https://github.com/szaghi/FoBiS) 3.8 or later (`fobis`), or [fpm](https://fpm.fortran-lang.org).
- **Dependencies:** [PENF](https://github.com/szaghi/PENF) (kinds) and [FUNDAL](https://github.com/szaghi/FUNDAL)
  (device memory and the portable directives), both fetched by `fobis fetch` and both compiled into SELAH's library.

## Build with FoBiS

```bash
git clone https://github.com/szaghi/SELAH && cd SELAH
fobis fetch                      # PENF and FUNDAL into src/third_party/
fobis build                      # the default mode: static library, gfortran, lib/gnu/static/libselah.a
fobis build --mode tests-gnu     # the library and the tests, into exe/
bash scripts/run_tests.sh        # run them
fobis build --lmodes             # every mode
```

## Modes

| Modes | What they build |
|---|---|
| `static-gnu`, `shared-gnu`, `static-nvf`, `shared-nvf`, `static-intel`, `shared-intel` | the library, optimised |
| `tests-gnu`, `tests-gnu-debug`, `tests-gnu-omp`, `tests-gnu-coverage` | the tests with gfortran: optimised, every runtime check, OpenMP, coverage |
| `tests-nvf`, `tests-nvf-debug`, `tests-intel`, `tests-intel-debug` | the tests with nvfortran and ifx, on the host |
| `tests-nvf-oac` | the tests with nvfortran and OpenACC, on the GPU |
| `*-r4` (`static-gnu-r4`, `tests-gnu-r4`, `tests-nvf-oac-r4`, …) | the same in single precision, into their own directories (`libselah_r4.a`, `exe-r4/`) |
| `bench-nvf-oac`, `bench-stepper-nvf-oac` (and their `-gnu` build checks) | the benchmarks |

The OpenACC modes take the compute capability from a varset: `--varset cc80` (A100, A30), `cc89`, `cc90`, or `iac`
(medooza: cc80 with CUDA 12.9).

## Backends

| Backend | How | Status |
|---|---|---|
| host, serial | any mode | tested with gfortran, nvfortran, ifx |
| host, OpenMP | `tests-gnu-omp` (`-fopenmp`): the kernels' loops become `parallel do` | tested, 4 threads |
| GPU, OpenACC | `tests-nvf-oac` (`-acc`, `DEV_OAC`, `COMPILER_NVF`) | tested on an NVIDIA A30 ([C5](/verification/contracts/device)) |
| GPU, OpenMP offload | the kernels carry the directives (FUNDAL's `DEV_OMP` macros) | **not built or tested**: no fobos mode defines `DEV_OMP` |

In a device build every array a kernel receives must be device memory, the coefficient arrays included, and
contiguous. The stepper allocates its own registers through FUNDAL, so they are device memory.

## Use from another FoBiS project

```ini
[dependencies]
deps_dir = src/third_party
PENF     = https://github.com/szaghi/PENF :: branch=master    # SELAH's own dependencies: fobis fetch is not transitive
FUNDAL   = https://github.com/szaghi/FUNDAL :: branch=main
SELAH    = https://github.com/szaghi/SELAH :: branch=main
```

Build your program with `-DFUNDAL_NO_MPI` (unless it uses FUNDAL's MPI handler) and with FUNDAL's `src/lib` among the
include directories: the kernels include `fundal.H`.

## Use from fpm

```toml
[dependencies]
SELAH = { git = "https://github.com/szaghi/SELAH", branch = "main" }
```

fpm fetches PENF and FUNDAL from SELAH's own manifest. Build with `-DFUNDAL_NO_MPI` unless your program uses FUNDAL's
MPI handler.
