<div align="center">

# SELAH

#### *the time integrators of ADAM, carved out as stage-wise operators, on every CPU and GPU*

### **S**tage-wise **E**volution **L**ibrary for **A**ccelerators and **H**PC

> Time-evolution operators for block-structured field solvers, in modern Fortran. Extracted from the Runge–Kutta
> machinery of ADAM and successor of FOODIE: verified scheme tableaux, bandwidth-bound kernels over the caller's
> arrays, and a stepper the caller drives one stage at a time, with no dynamic dispatch where the cells are swept.

[![CI](https://github.com/szaghi/SELAH/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/SELAH/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/SELAH/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/SELAH/actions/workflows/matrix.yml)
[![Fortran](https://img.shields.io/badge/Fortran-modern-734f96?logo=fortran&logoColor=white)](https://fortran-lang.org/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue)](LICENSE.gpl3.md)
[![License: BSD-2](https://img.shields.io/badge/license-BSD--2--Clause-blue)](LICENSE.bsd-2.md)
[![License: BSD-3](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.bsd-3.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE.mit.md)

</div>

---

> [!WARNING]
> **Pre-release (M0).** The repository holds the build, test and documentation scaffolding only. The design and the
> milestones are in [notes/PLAN.md](notes/PLAN.md).

## The name

*Selah* is the mark that appears 71 times in the Psalms and 3 times in Habakkuk, usually read as a pause: the music
stops and the singer breathes. SELAH's stage boundaries are such pauses. The caller's orchestrator acts there (ghost
exchange, AMR seam fill, reflux) before the integrator goes on.

## Design in one paragraph

Four layers (details in [notes/PLAN.md](notes/PLAN.md)):

| Layer | What it holds |
|---|---|
| **L0** tableaux | Butcher, Shu–Osher, low-storage and splitting coefficients, verified in quadruple precision |
| **L1** kernels | stage and update operations over the caller's arrays: interior-only boxes and flat arrays, one source for host OpenMP, OpenACC and OpenMP offload |
| **L2** stepper | `begin_stage` / `end_stage` / `finish_step`, driven by the caller, under a documented stage contract |
| **L3** adapters | in the consumer (ADAM), not in SELAH |

## Build

```bash
fobis fetch                      # PENF
fobis build --mode tests-gnu     # the library and the tests
bash scripts/run_tests.sh
fobis build                      # the static library, lib/gnu/static/libselah.a
fobis build --lmodes             # every mode (gfortran, nvfortran, ifx; R8P and R4P; OpenACC)
```

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@gmail.com
> *Chief Yak Shaver & Keeper of the Time Step*: believes every problem is one more library away from being solved,
> and keeps proving it

**[Claude](https://claude.ai)** (Anthropic)
> *Tireless Rubber Duck & Professional Pedant*: AI pair programmer; checks the order conditions so humans don't have
> to, and has an opinion on every comma

## License

SELAH is available under your choice of license:

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](LICENSE.mit.md) |

© [Stefano Zaghi](https://github.com/szaghi)
