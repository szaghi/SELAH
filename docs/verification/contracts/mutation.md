---
title: C6 Mutation testing
icon: selah
---

# C6 Mutation testing

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C6 | Contracts | C | CPU (gfortran 16.2, `tests-gnu`) | no | `scripts/mutation_test.py` |

## Case and scope

A test suite can pass because the code is right, or because the tests cannot see the defect. Mutation testing asks the
second question directly. It plants one defect at a time in a scratch copy of the tree, rebuilds, runs the whole suite,
and the defect (the *mutant*) is *killed* when at least one program fails. A surviving mutant is a defect the suite is
blind to.

## Setup

| Mutant | What it changes |
|---|---|
| ssp22-weight, ssp33-weight, rk4-weight | a Butcher weight in the fifth digit |
| ssp54-shu-osher-beta, ssp104-alpha | a Shu–Osher coefficient (converted at run time) |
| ssp33-two-register-c | a two-register coefficient |
| lsrk33-A, lsrk54-B | a 2N coefficient (lsrk54-B in the 13th digit) |
| yoshida-w1, blanes-moan-s6-a1 | a splitting coefficient |
| ssp-coefficient-max | the SSP coefficient as the maximum instead of the minimum ratio |
| axpy-sign, lincomb-beta, stage-state-dt, step-update-dq-zero, shu-osher-2r-swap | a kernel body |
| williamson-a0-reads | the 2N stage with $A = 0$ reads its register again |
| mask-sign | the mask's sign convention inverted in one kernel |
| box-ghost-bound | the box's allocation lower bound off by one |
| stage-state-copy-n1 | the copy shortcut taken for one register instead of none |
| registers-butcher | one register too few for the Butcher form |
| rhs-slot | the residual register is the stage state's |
| stage-weight-c | `stage_weight` returns $c_k$ instead of $b_k$ |
| ketcheson-euler-step, ketcheson-combination | Ketcheson's recipe: the Euler step, the mid-step combination |
| stages-on-q-all | every stage treated as evaluated on q |
| two-register-no-qn | the two-register form forgets to keep $q^n$ |
| splitting-subflow | the splitting sub-flow scaled by 1.0001 |

## Program

```bash
python3 scripts/mutation_test.py [--only NAME ...]     # about half an hour for the 28 mutants
```

## Results

**28 of 28 mutants killed** (2026-11-17; the run record is `notes/mutation-run-2026-11-17.txt`).

## Findings

- The first run killed 27. The survivor, `ssp104-alpha`, changed the coefficient of $u^n$ in ssp-104's last Shu–Osher
  stage. Nothing in the library reads it: the Butcher conversion skips the row of $u^n$, and the SSP coefficient skips
  $\beta = 0$. So a wrong value in the public `shu_osher()` result was invisible. `selah_test_schemes` now checks that
  every Shu–Osher form is consistent (each row of $\alpha$ sums to 1), and the mutant is killed.
- Most mutants are caught by several programs at once: the golden file, the orders and the contract checks overlap.
