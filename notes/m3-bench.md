# M3: SELAH's kernels against ADAM's on the A30

*2026-11-17, night m3-device. The program is `src/bench/selah_bench_kernels.F90`, mode `bench-nvf-oac --varset iac`,
run by `scripts/medooza/selah-tests.slurm`.*

## Setup

- **Machine:** medooza, `gnode01`, one NVIDIA A30 (24 GB, driver 575.51.03, 1440 MHz).
- **Build:** nvfortran 26.9, OpenACC, `-gpu=cc80,cuda12.9`, release profile.
- **Field:** a FLUME-like field in ADAM's device layout `(b, i, j, k, v)`: 5 components, 16³ interior cells with
  3 ghost cells, 256 blocks. That is 109 MB per field.
- **ADAM side:** ADAM's `adam_fnl_rk_kernels.F90`, the byte-exact copy in `src/tests/` (the version at ff879596).
- **Timing:** every launch is synchronous: `!$acc wait` before the clock is read, `system_clock` and never `cpu_time`.
  - First, one unreported pass of every operation.
  - Then 3 passes; SELAH is timed first on passes 1 and 3, ADAM first on pass 2.
  - Each operation: 3 warm-up launches, then the median of 7 trials of 20 launches.
- **Bandwidth:** each field read or written counts once per updated value, over the interior cells only. The reference
  is a STREAM triad `a = b + s c` over arrays the size of one field.

## Results (job 1368089)

STREAM triad: **816 GB/s**.

| Operation | Fields | SELAH [ms] | ADAM [ms] | SELAH/ADAM | SELAH [GB/s] | of triad |
|---|---|---|---|---|---|---|
| stage state, stage 5 of 5 | 6 | 0.469 | 0.484 | 0.969 | 536 | 0.66 |
| step update, 5 registers | 7 | 0.583 | 0.585 | 0.995 | 504 | 0.62 |
| copy residual to register | 2 | 0.164 | 0.164 | 0.995 | 513 | 0.63 |
| two-register stage | 4 | 0.220 | 0.220 | 1.001 | 763 | 0.93 |
| stage state, masked (half the blocks solid) | 7 | 0.598 | 0.597 | 1.002 | 491 | 0.60 |
| step update + dq (fused) | 8 | 0.610 | — | — | 550 | 0.67 |
| **whole SSP(5,4) step** (integrator kernels only) | 37 / 43 | 3.637 | 4.576 | **0.795** | 427 | 0.52 |

The table shows pass 2, timed with ADAM first. Passes 1 and 3 agree with it to within 0.012 on every ratio; the
largest spread is on the masked stage state (1.002–1.014).

The whole-step row is the comparison that matters:
- **SELAH:** five stage states built from `q`, five residuals copied into registers, one update. That is 37 field
  accesses.
- **ADAM:** `initialize_stages` copies `q` into all five buffers, then five in-place stages, five assigns and one
  update. That is 43 field accesses.

SELAH is 20% faster, more than the 14% the access counts predict.

## Findings

1. **Parity holds.** Every SELAH kernel is within 0.97–1.015 of ADAM's on the same work.
2. **SELAH's design is faster where the work differs.** SELAH builds each stage state directly from `q`, so ADAM's
   copy of `q` into every stage buffer disappears.
3. **A benchmark artifact, now fixed.** The first operations timed after the allocations ran slow for one pass. Job
   1368086, with one pass and no warm-up pass, put the stage state at 1.31× ADAM. Job 1368088 showed pass 1 unstable
   and passes 2–3 identical. The unreported warm-up pass removes it: in job 1368089 all three passes agree.
4. **The plan's "≥ 80% of the triad" target is not met, by SELAH or by ADAM.** The kernels reach 0.60–0.67 of the
   triad; the two-register stage reaches 0.93. Two properties of the access pattern explain the gap:
   - each kernel streams 6 to 8 arrays at once, where the triad streams 3;
   - an interior-only box covers 16³ of 22³ cells per block, so each contiguous run of 256 blocks is followed by a
     jump.

   Closing it would take kernel-shape work, such as tiling or vector-length tuning, that ADAM's kernels don't do either.
   This is a decision for the owner, not a defect of the extraction.

## Jobs

| Job | Content | Outcome |
|---|---|---|
| 1368084 | tests-nvf-oac, first device tests | device kernels pass; the smoke test read device memory on the host (fixed) |
| 1368085 | tests-nvf-oac | 9/9; device parity with ADAM bitwise |
| 1368086 | bench, one pass | stage state 1.31× ADAM (a first-timed artifact, see finding 3) |
| 1368088 | bench, 3 passes, order alternated | pass 1 unstable, passes 2–3 steady |
| 1368089 | bench, warm-up pass + 3 passes | the table above |

Not measured: the WSL RTX (cc89). That GPU is reserved for ADAM's tests; running the bench there with
`--varset wsl` is the owner's call.
