!< SELAH, L1: the kernels, one generic name per operation over a box (rank 5) or a flat array (rank 1).
module selah_kernels
!< SELAH, L1: the kernels, one generic name per operation over a box (rank 5) or a flat array (rank 1).
!<
!<| Generic | Operation | Box (selah_kernels_box) | Flat (selah_kernels_flat) |
!<|---------|-----------|-------------------------|---------------------------|
!<| `selah_copy`                | y = x                                   | masked | unmasked |
!<| `selah_axpy`                | y = y + alpha x                         | masked | unmasked |
!<| `selah_lincomb`             | y = alpha x + beta y                    | masked | unmasked |
!<| `selah_stage_state`         | u = q + sum_k dt coef(k) regs(k)        | masked [, inactive u = q] | unmasked |
!<| `selah_step_update`         | q = q + dt sum_k coef(k) regs(k) [, dq] | masked, dq = 0 inactive | unmasked |
!<| `selah_shu_osher_2r_stage`  | q = a q0 + b q + dt c r                 | masked | unmasked |
!<| `selah_williamson_2n_stage` | dqr = a dqr + dt r; q = q + b dqr       | masked | unmasked |
!<
!< The box form takes the box first; the flat form has no box and no mask.
use selah_kernels_box
use selah_kernels_flat
implicit none
private
public :: selah_copy, selah_axpy, selah_lincomb, selah_stage_state, selah_step_update
public :: selah_shu_osher_2r_stage, selah_williamson_2n_stage

interface selah_copy !< y = x.
   module procedure copy_box, copy_flat
endinterface selah_copy

interface selah_axpy !< y = y + alpha x.
   module procedure axpy_box, axpy_flat
endinterface selah_axpy

interface selah_lincomb !< y = alpha x + beta y.
   module procedure lincomb_box, lincomb_flat
endinterface selah_lincomb

interface selah_stage_state !< Stage state: u = q + sum_k dt coef(k) regs(k), k = 1..n (a copy for n = 0).
   module procedure stage_state_box, stage_state_flat
endinterface selah_stage_state

interface selah_step_update !< Step update: q = q + dt sum_k coef(k) regs(k), and the step residual dq if present.
   module procedure step_update_box, step_update_flat
endinterface selah_step_update

interface selah_shu_osher_2r_stage !< Two-register Shu–Osher stage: q = a q0 + b q + dt c r.
   module procedure shu_osher_2r_stage_box, shu_osher_2r_stage_flat
endinterface selah_shu_osher_2r_stage

interface selah_williamson_2n_stage !< Williamson 2N stage: dqr = a dqr + dt r; q = q + b dqr.
   module procedure williamson_2n_stage_box, williamson_2n_stage_flat
endinterface selah_williamson_2n_stage
endmodule selah_kernels
