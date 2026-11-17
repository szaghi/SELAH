---
title: Performance
---

# Performance

SELAH's kernels against ADAM's FNL kernels, the code they were carved from, on an NVIDIA A30 (medooza, driver 575.51.03,
nvfortran 26.9, OpenACC `-gpu=cc80,cuda12.9`). The program is `src/bench/selah_bench_kernels.F90`. The method and every
job are in `notes/m3-bench.md` of the repository.

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
directly from the solution: 37 field accesses per step, against ADAM's 43, which first copy the solution into every
stage buffer. The step residual comes in the same pass as the update, where ADAM needs a separate kernel.

**Neither code reaches 80% of the triad.** The kernels stream 6 to 8 arrays over the interior of each block, a gap
every 256 contiguous values; the triad streams 3 contiguous arrays.
