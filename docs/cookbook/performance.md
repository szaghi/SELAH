---
title: Performance
---

# Performance

[[toc]]

## Time a device kernel

A GPU launch returns before the kernel finishes. Time it with `system_clock` after an explicit wait, never with
`cpu_time` (CPU time of the host thread). Warm up first, take the median of several trials, and alternate the order
of the things you compare. The pattern of SELAH's benchmarks (`src/bench/selah_bench_stepper.F90`):

```fortran
do k=1, WARM                       ! warm-up launches, not timed
   call launch(op)
enddo
!$acc wait
do k=1, TRIALS
   call system_clock(c0, rate)
   do l=1, LAUNCHES
      call launch(op)
   enddo
   !$acc wait                      ! the kernels have finished
   call system_clock(c1)
   t(k) = real(c1 - c0, R8P) / real(rate, R8P) / LAUNCHES
enddo
ms = median(t) * 1.e3_R8P
```

The benchmarks also run one unreported pass of every operation before the timed passes: on the A30 the first
operations timed after the allocations ran up to 30% slow (notes/m3-bench.md).

## The paths that save memory traffic

An explicit integrator does almost no arithmetic per value: its speed is memory traffic. SELAH's design removes field
passes that ADAM's step makes:

| Path | How | Saving per SSP(5,4) step (field accesses) |
|---|---|---|
| stage states built from q | `begin_stage` sums from q and the residuals; no copy of q into the stage buffers | 10 |
| residuals written in place | the caller writes `rhs_register()`; no assign copy | 10 |
| step residual fused | `finish_step(dq=)` writes dq in the update's pass | 5 |
| stage 1 on q | `stages_on_q=.true.` | 2 |
| low-storage forms | 2 registers; few-stream kernels that reach a higher fraction of the bandwidth | memory: s − 1 registers |

On an A30 the stepper's SSP(5,4) step with dq takes 0.52 of ADAM FLUME's, and the stepper itself costs nothing
([Performance](/guide/performance)).
