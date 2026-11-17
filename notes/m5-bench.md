# M5: whole steps through the stepper, on the A30

*2026-11-17, night m5-fused. The program is `src/bench/selah_bench_stepper.F90`, mode `bench-stepper-nvf-oac --varset
iac`, run by `scripts/medooza/selah-tests.slurm`.*

## Setup

- **Machine, build, field, timing:** as in [m3-bench.md](m3-bench.md).
  - Machine: medooza `gnode01`, one A30.
  - Build: nvfortran 26.9, `-gpu=cc80,cuda12.9`.
  - Field: 5 components, 16³ interior cells with 3 ghosts, 256 blocks, layout V_LAST, 109 MB per field.
  - Timing: an unreported warm-up pass, then 3 passes with the order alternated; for each operation, the median of 7
    trials of 20 launches, `!$acc wait` before `system_clock`.
- **Operations:** whole steps of the integrator alone.
  - SELAH is driven through the stepper as a solver drives it: `begin_step`, then `begin_stage`/`end_stage` for each
    stage, then `finish_step`.
  - The residual is the caller's work and is the same on both sides, so it is left in the register.
  - The registers are the caller's (`registers=`), initialised on the device.
- **ADAM side:**
  - the FNL kernels of `src/tests/adam_fnl_rk_kernels.F90`;
  - FLUME's step residual kernel `compute_rk_ssp_residual_dev`, copied verbatim into
    `src/bench/adam_flume_rk_residual.F90` (`adam_flume_fnl_kernels.F90:578-614`, sha256 `974e532a…`).

## The traffic model

These are field accesses per step over the interior, each field read or written once per updated value. They were
counted from the code before the runs.

| Path | Accesses | Breakdown |
|---|---|---|
| ADAM, SSP(5,4) | 47 | q copied into the 5 stage buffers: 10. In-place stage s, reading its buffer and the s−1 earlier ones, then writing it: 2+3+4+5+6. Five assigns: 10. Update: 7 |
| ADAM FLUME, SSP(5,4) with the step residual | 53 | 47, plus the residual kernel, which rereads the 5 stages and writes dq (ADAM#128 D4) |
| SELAH stepper, Butcher SSP(5,4) | 27 | Stage states from q: 2+3+4+5+6. The residual is written into the register in place, with no assign. Update: 7 |
| SELAH stepper, Butcher SSP(5,4), dq fused | 28 | 27, plus 1 for dq written by the update |
| SELAH, stage 1 evaluated on q | 26 | 28 without the stage-1 copy (not in the stepper's API, see below) |
| SELAH stepper, Butcher SSP(10,4) | 77 | 2+…+11, plus the update, 12 |
| SELAH stepper, Ketcheson SSP(10,4) (two registers) | 39 | begin_step copy 2; 9 axpy at 3 each; stage 5, two lincomb at 3 each; stage 10, 4 |

**Two corrections to earlier counts:**
- **M3's count of ADAM's step was wrong**: it is 47 accesses, not 43, because the in-place stage also reads its own
  buffer. With 47, M3's measured 0.795 sits on the access ratio 37/47 = 0.787. The "faster than the counts predict"
  in M3 was a miscount.
- **PLAN's "49 → 28" was built on the same 43** (43 + 6). With the corrected count it is 53 → 28.

**Criterion, set before the runs:** for each pair, the measured time ratio is within ±10% of the access ratio.

## Results

Job 1368100 is on the final code: the empty-sum stage state is a copy (finding 3). The table shows pass 2; the
passes agree to within 0.004 on every ratio.

| Operation | Accesses | Time [ms] | GB/s |
|---|---|---|---|
| ADAM, SSP(5,4) step | 47 | 4.569 | 432 |
| ADAM FLUME, step + residual kernel | 53 | 4.948 | 449 |
| SELAH stepper, SSP(5,4) Butcher | 27 | 2.688 | 421 |
| SELAH stepper, SSP(5,4) Butcher + dq fused | 28 | 2.714 | 433 |
| SELAH, SSP(5,4) + dq, stage 1 on q (kernels) | 26 | 2.545 | 429 |
| SELAH stepper, SSP(10,4) Butcher | 77 | 6.399 | 505 |
| SELAH stepper, SSP(10,4) Ketcheson | 39 | 2.322 | 705 |
| SELAH, SSP(5,4) + dq (kernels only, the control) | 28 | 2.710 | 433 |

| Pair | Access ratio | Time ratio | Time / access | Criterion |
|---|---|---|---|---|
| SELAH + dq / ADAM FLUME | 0.528 | **0.548** | 1.038 | met |
| SELAH / ADAM | 0.574 | **0.588** | 1.024 | met |
| Kernels only / stepper (+ dq) | 1.000 | 0.999 | 0.999 | — |
| Stage 1 on q / SELAH + dq (kernels) | 0.929 | 0.939 | 1.011 | met |
| Ketcheson / Butcher, SSP(10,4) | 0.506 | **0.363** | 0.716 | **missed** |

**Equal simulated time.** The time step is SSP-limited, Δt ∝ C, with C = 1.508 for SSP(5,4) and 6 for SSP(10,4).
Comparing SSP(10,4) Ketcheson with SSP(5,4) Butcher per unit simulated time:
- the integrator costs **0.217** as much (accesses predict 0.363);
- the residual is evaluated **0.503** times as often (10/6 against 5/1.508, exact);
- registers: 2 instead of 6, so 3 fields with q instead of 7.

A solver's ratio lies between 0.217 and 0.503, depending on the share of the step its residual takes. When the residual
dominates the step, as a WENO residual does, the gain approaches 2× per unit time. This holds only where the stable step is the SSP limit.
When accuracy limits the step, the comparison is not this one, and it was not measured.

## Findings

1. **The fused paths deliver.**
   - Through the stepper, SSP(5,4) with the step residual costs **0.548** of ADAM FLUME's step with its residual
     kernel.
   - Without dq it costs 0.588 of ADAM's step.
   - In both cases the time follows the traffic to within 4%.

   Three things are gone: the copy of q into every stage buffer, the assign of each residual, and the rereading of
   the stages for dq.
2. **The stepper is free.** Driving the same sequence of kernels without the stepper changes the time by 0.1%.
3. **The empty-sum stage state is now a copy.**
   - Before the fix (job 1368098), the state of stage 1 (u = q + 0) ran at 291 GB/s and a plain copy at 513 GB/s, on
     the same two fields.
   - `stage_state_box` and `stage_state_flat` now hand n = 0 to the copy kernel. That made the SSP(5,4) step 4.3%
     faster (2.836 → 2.714 ms), and SELAH + dq against ADAM FLUME went from 0.573 to 0.548.
   - ADAM's first stage is the same empty sum, run in place after the copy.
4. **The access count misses Ketcheson by 28%, and the miss has a cause.**
   - Ketcheson's kernels stream 3 or 4 arrays and run at 705 GB/s, 0.86 of the triad.
   - The Butcher stage states of SSP(10,4) stream up to 12 arrays and run at 505 GB/s.
   - Counting accesses ignores the bandwidth a kernel shape reaches. M3 had already shown this: the two-register stage
     reached 0.93 of the triad, the stage states 0.60–0.66.
   - The miss is in Ketcheson's favour, and the criterion was not retuned after the run.
5. **Stage 1 evaluated on q would save another 6.1% of the step**, for 7.1% of the accesses. The state of stage 1 is
   q itself, so the copy into the stage-state register is avoidable. The stepper's API does not allow it yet: see the
   proposal in the night log.

## Jobs

| Job | Content | Outcome |
|---|---|---|
| 1368096 | bench, 7 operations | criterion met for the two ADAM pairs (1.085, 1.071), missed for Ketcheson (0.703); the stage-1 row mixed two effects |
| 1368097 | + the kernels-only control | the stepper costs nothing; stage 1 on q saved more than its accesses |
| 1368098 | + the empty-sum stage state against a copy | 1.759: the cause of the extra saving |
| 1368099 | tests-nvf-oac, the copy fix | 10/10 |
| 1368100 | bench, the copy fix | the tables above |
