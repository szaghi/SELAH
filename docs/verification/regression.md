---
title: Regression
---

# Regression

What is pinned between versions of SELAH, how it is checked, and the rules for changing it.

| What | Where | Checked by | Rule |
|---|---|---|---|
| the values of every scheme, form and kernel | `src/tests/golden/selah_golden_v1.txt`, `selah_golden_r4_v1.txt` | `selah_test_golden`, every build | never edited; a deliberate change of results gets a new file (v2) and a new default path ([C4](./contracts/golden)) |
| ADAM's kernels, the oracle of equivalence | `src/tests/adam_fnl_rk_kernels.F90` (byte-exact, sha256 `a693ca35…`), `src/bench/adam_flume_rk_residual.F90` (verbatim, sha256 of its ADAM source `974e532a…`) | `selah_test_kernels_adam`, `selah_test_device_adam` | replaced only by a new verbatim copy of a later ADAM, with the version in the header ([C1](./contracts/adam-equivalence)) |
| the documentation's examples and their outputs | `docs/examples/snippets/`, `docs/examples/output/` | `scripts/docs_examples.sh`; the Docs examples workflow fails when the regenerated files differ | regenerated, never edited |
| the verification tables and figures | `docs/verification/_tables/`, `docs/public/verification/`, `docs/public/guide/` | `scripts/run_vv.sh` regenerates them from the programs | regenerated, never edited; the tables of a report run are byte-identical to the committed ones |
| the suite's power to see defects | `scripts/mutation_test.py` | 28 of 28 mutants killed ([C6](./contracts/mutation)) | a new feature gets its mutants |
| the A30 measurements | `notes/m3-bench.md`, `notes/m5-bench.md`, `notes/m5-bench-data/` | the benchmarks, by hand | new measurements get new files and job numbers |

## Continuous integration

On every push the workflows run:
- the test suite in the fortran-toolchains images: gfortran, nvfortran and oneAPI (`matrix.yml`);
- the fpm build and tests (`fpm.yml`);
- the installation check (`install.yml`);
- the documentation's examples (`docs-examples.yml`).

The documentation is built and published from the CI's coverage data (`docs.yml`). The GPU runs are not in CI. They
are made on medooza with `scripts/medooza/selah-tests.slurm`, and their job numbers are recorded on the pages that use
them.
