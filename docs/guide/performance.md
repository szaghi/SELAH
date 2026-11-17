---
title: Performance
---

# Performance

SELAH's kernels, and whole steps through the [stepper](./stepper), against ADAM's FNL kernels, the code they were carved from, on an NVIDIA A30 (medooza, driver 575.51.03,
nvfortran 26.9, OpenACC `-gpu=cc80,cuda12.9`). The programs are `src/bench/selah_bench_kernels.F90` (kernels, M3) and
`src/bench/selah_bench_stepper.F90` (whole steps, M5). The method and every job are in `notes/m3-bench.md` and
`notes/m5-bench.md` of the repository.

## Setup

- **Field:** 5 components, 16³ interior cells with 3 ghost cells, 256 blocks, in ADAM's device layout; 109 MB per field.
- **Timing:** every launch is synchronous: `!$acc wait` before `system_clock`.
  - First, an unreported warm-up pass of every operation.
  - Then 3 passes, alternating which side is timed first.
  - Each operation: the median of 7 trials of 20 launches.
- **Bandwidth reference:** a STREAM triad over arrays the size of one field, **816 GB/s**.

## Results

| Operation | SELAH/ADAM | SELAH [GB/s] | of triad |
|---|---|---|---|
| stage state (stage 5 of SSP(5,4)) | 0.969 | 536 | 0.66 |
| step update (5 registers) | 0.995 | 504 | 0.62 |
| copy of a residual into a register | 0.995 | 513 | 0.63 |
| two-register stage | 1.001 | 763 | 0.93 |
| stage state, masked | 1.002–1.014 | 491 | 0.60 |
| step update fused with the step residual | — | 550 | 0.67 |
| **whole SSP(5,4) step** | **0.795** | 427 | 0.52 |

**Per kernel, SELAH matches ADAM**, to within 3%. **On a whole step SELAH is 20% faster.** It builds each stage state
directly from the solution: 37 field accesses per step, against ADAM's 47, which first copy the solution into every
stage buffer (M3 counted 43: the in-place stage also reads its buffer). The time ratio, 0.795, sits on the access
ratio, 0.787. This row still copies each residual into its register, as ADAM does; the stepper (below) does not.

**Neither code reaches 80% of the triad.** The kernels stream 6 to 8 arrays over the interior of each block, a gap
every 256 contiguous values; the triad streams 3 contiguous arrays.

## Whole steps through the stepper

The stepper driven as a solver drives it, the residual left in its register (the caller's work, the same on both sides).
Field accesses per step over the interior, counted from the code before the runs; job 1368100.

| Step | Accesses | Time [ms] | GB/s |
|---|---|---|---|
| ADAM SSP(5,4) | 47 | 4.569 | 432 |
| ADAM FLUME SSP(5,4) + its step residual kernel | 53 | 4.948 | 449 |
| SELAH stepper SSP(5,4), Butcher | 27 | 2.688 | 421 |
| SELAH stepper SSP(5,4), Butcher, step residual fused | 28 | 2.714 | 433 |
| the same, `stages_on_q=.true.` (stage 1 evaluated on q; job 1368116) | 26 | 2.552 | 427 |
| SELAH stepper SSP(10,4), Butcher | 77 | 6.399 | 505 |
| SELAH stepper SSP(10,4), Ketcheson's two registers | 39 | 2.322 | 705 |

- **With the step residual, a SELAH step costs 0.548 of ADAM FLUME's** (accesses: 0.528); without it, 0.588 of
  ADAM's (0.574). Three passes of the field are gone: the copy of the solution into every stage buffer, the copy of
  every residual into its register (the caller writes it in place, into `rhs_register()`), and the rereading of every
  stage for the step residual (fused into the update).
- **The stepper costs nothing**: the same kernels called without it take 0.999 of the time.
- **Or no copy at all.** With `stages_on_q=.true.` the caller evaluates stage 1 on `q` and the stepper skips the
  copy: 0.940 of the step, 0.515 of ADAM FLUME's (see [Stepper](./stepper#stages-on-q)).
- **The state of the first stage is a copy.** The summing kernel with nothing to sum ran at 0.57 of a copy's speed;
  `selah_stage_state` now hands that case to the copy kernel, which made the SSP(5,4) step 4.3% faster.

**SSP(10,4) in two registers.** Its kernels stream 3 or 4 arrays and reach 0.86 of the triad, so it costs 0.363 of the
Butcher form of the same scheme, below its access ratio (0.506). At equal simulated time, with the time step at the
SSP limit (coefficient 6 against SSP(5,4)'s 1.508), against SSP(5,4):

| | SSP(10,4) Ketcheson / SSP(5,4) Butcher |
|---|---|
| integrator time | 0.217 |
| residual evaluations | 0.503 |
| stored fields (with q) | 3 / 7 |

A solver's ratio lies between the first two rows, by the share of its step the residual takes. Where accuracy rather
than the SSP limit sets the step, this comparison does not apply; it has not been measured.
