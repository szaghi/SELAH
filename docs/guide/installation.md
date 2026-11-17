---
title: Installation
---

# Installation

## Requirements

- a Fortran 2018 compiler: gfortran (16.2 is the reference), nvfortran (26.9; OpenACC modes), ifx;
- [FoBiS](https://github.com/szaghi/FoBiS) 3.8 or later (`fobis`), or [fpm](https://fpm.fortran-lang.org);
- [PENF](https://github.com/szaghi/PENF), fetched by `fobis fetch`.

## Build with FoBiS

```bash
git clone https://github.com/szaghi/SELAH && cd SELAH
fobis fetch                      # PENF, compiled from its sources by each mode
fobis build                      # the default mode: static library, gfortran, lib/gnu/static/libselah.a
fobis build --mode tests-gnu     # the library and the tests, into exe/
bash scripts/run_tests.sh        # run them
fobis build --lmodes             # every mode
```

The single-precision modes (`*-r4`) build PENF and SELAH with `R_P = R4P` into their own directories
(`libselah_r4.a`, `exe-r4/`). The OpenACC modes take the compute capability from a varset:
`fobis build --mode tests-nvf-oac --varset cc80`.

## Use from another FoBiS project

```ini
[dependencies]
SELAH = https://github.com/szaghi/SELAH :: branch=main
```

## Use from fpm

```toml
[dependencies]
SELAH = { git = "https://github.com/szaghi/SELAH", branch = "main" }
```
