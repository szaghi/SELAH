!< SELAH benchmark (M5): whole steps through the stepper against ADAM's FNL step, and SSP(10,4) in two registers.

#include "fundal.H"

program selah_bench_stepper
!< SELAH benchmark (M5): whole steps through the stepper against ADAM's FNL step, and SSP(10,4) in two registers.
!<
!< The field is the one of selah_bench_kernels: layout V_LAST, 5 components, 16^3 interior cells with 3 ghosts, NB blocks
!< (default 256, ~109 MB per field). Each operation is one whole step of the integrator alone: the stepper is driven as
!< a solver drives it (begin_step, begin_stage/end_stage per stage, finish_step), with the residual left in the
!< register, since the residual is the caller's work and is the same on both sides. The registers are the caller's
!< (`registers=`), initialised on the device, so that every value stays finite over the launches.
!<
!< Field accesses per step over the interior (each field read or written once per updated value):
!<
!<| Operation | Accesses |
!<|---|---|
!<| ADAM SSP(5,4): q copied into the 5 stage buffers, 5 in-place stages, 5 assigns, the update | 47 |
!<| ADAM FLUME FNL: the same + `compute_rk_ssp_residual_dev` (rereads the 5 stages, writes dq) | 53 |
!<| SELAH stepper, Butcher SSP(5,4): stage states from q, residuals written in place, the update | 27 |
!<| SELAH stepper, Butcher SSP(5,4), the step residual dq fused in the update | 28 |
!<| SELAH kernels, Butcher SSP(5,4) + dq, stage 1 evaluated on q (no stage-1 copy; not in the stepper's API) | 26 |
!<| SELAH stepper, Butcher SSP(10,4) | 77 |
!<| SELAH stepper, Ketcheson SSP(10,4), two registers | 39 |
!<| SELAH kernels, Butcher SSP(5,4) + dq, the stepper's sequence without the stepper (the control of the stage-1 row) | 28 |
!<| The state of stage 1 alone (no register summed), and a plain copy doing the same work | 2, 2 |
!<
!< Timing as in selah_bench_kernels: an unreported pass of every operation, then PASSES passes with the order of each
!< pair alternated; each operation 3 warm-up launches and the median of 7 trials of 20 launches, `!$acc wait` before
!< the clock is read (system_clock, never cpu_time).
!<
!< Usage: selah_bench_stepper [NB [PASSES]]
use fundal, only : dev_alloc, dev_assign_to_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use penf, only : I4P, I8P, R8P
use selah
use adam_fnl_rk_kernels
use adam_flume_rk_residual
implicit none
integer(I4P), parameter :: NV=5, NI=16, NJ=16, NK=16, NGC=3, WARM=3, TRIALS=7, LAUNCHES=20, NOPS=10
real(R8P),    parameter :: DT=1.e-3_R8P
character(*), parameter :: OP_NAME(NOPS) = [character(48) :: &
                           'ADAM SSP(5,4) step',                               &
                           'ADAM FLUME SSP(5,4) step + residual kernel',       &
                           'SELAH stepper SSP(5,4) Butcher',                   &
                           'SELAH stepper SSP(5,4) Butcher + dq fused',        &
                           'SELAH SSP(5,4) + dq, stage 1 on q (kernels)',      &
                           'SELAH stepper SSP(10,4) Butcher',                  &
                           'SELAH stepper SSP(10,4) Ketcheson 2 registers',    &
                           'SELAH SSP(5,4) + dq (kernels, the control of 5)',  &
                           'stage state of stage 1 (no register): u = q + 0',  &
                           'copy u = q']
integer(I4P), parameter :: OP_FIELDS(NOPS) = [47, 53, 27, 28, 26, 77, 39, 28, 2, 2]
integer(I4P)            :: nb
real(R8P), pointer      :: q_d(:,:,:,:,:)=>null(), r_d(:,:,:,:,:)=>null(), dq_d(:,:,:,:,:)=>null()
real(R8P), pointer      :: u_d(:,:,:,:,:)=>null()
real(R8P), pointer      :: q_rk_d(:,:,:,:,:,:)=>null(), reg54_d(:,:,:,:,:,:)=>null(), reg104_d(:,:,:,:,:,:)=>null()
real(R8P), pointer      :: reg2_d(:,:,:,:,:,:)=>null()
real(R8P), pointer      :: alph_d(:,:)=>null(), beta_d(:)=>null(), at_d(:,:)=>null()
type(butcher_tableau)   :: bt54, bt104
type(box)               :: bx
type(stepper)           :: st54, st104, stk
real(R8P)               :: cells, t(NOPS), ratio
character(32)           :: arg
integer(I4P)            :: ierr, passes, pass, op

nb = 256 ; passes = 3
if (command_argument_count() > 0) then
   call get_command_argument(1, arg) ; read(arg, *) nb
endif
if (command_argument_count() > 1) then
   call get_command_argument(2, arg) ; read(arg, *) passes
endif
#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
#endif
bt54 = butcher(SCHEME_SSP_54) ; bt104 = butcher(SCHEME_SSP_104)
bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, nb)
cells = real(NI, R8P) * NJ * NK * nb * NV ! updated values per field
call field(q_d, 1._R8P) ; call field(r_d, 0.5_R8P) ; call field(dq_d, 0._R8P) ; call field(u_d, 0._R8P)
call registers(q_rk_d, bt54%stages) ; call registers(reg54_d, bt54%stages + 1_I4P)
call registers(reg104_d, bt104%stages + 1_I4P) ; call registers(reg2_d, 2_I4P)
call dev_assign_to_device(dst=alph_d, src=bt54%a) ; call dev_assign_to_device(dst=beta_d, src=bt54%b)
call dev_assign_to_device(dst=at_d, src=transpose(bt54%a)) ! column s = row s of A, contiguous
call st54%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q_d), ub=ubound(q_d), registers=reg54_d)
call st104%initialize(SCHEME_SSP_104, FORM_BUTCHER, lb=lbound(q_d), ub=ubound(q_d), registers=reg104_d)
call stk%initialize(SCHEME_SSP_104, FORM_KETCHESON_2R, lb=lbound(q_d), ub=ubound(q_d), registers=reg2_d)

print '(A,I0,A,F7.1,A)', 'blocks ', nb, ', field ', 8._R8P*NV*(NI+2*NGC)*(NJ+2*NGC)*(NK+2*NGC)*nb/1.e6_R8P, ' MB'
print '(A,F6.3,A,F6.3)', 'SSP coefficients: SSP(5,4) ', bt54%ssp_coefficient, ', SSP(10,4) ', bt104%ssp_coefficient
do op=1, NOPS ! an unreported pass of every operation first (see selah_bench_kernels)
   t(op) = timed(op)
enddo
do pass=1, passes
   print '(A,I0,A)', 'pass ', pass, merge(' (SELAH timed first)', ' (ADAM timed first) ', mod(pass, 2) == 1)
   if (mod(pass, 2) == 1) then
      do op=NOPS, 1, -1
         t(op) = timed(op)
      enddo
   else
      do op=1, NOPS
         t(op) = timed(op)
      enddo
   endif
   print '(A)', '| operation | accesses | time [ms] | GB/s |'
   print '(A)', '|---|---|---|---|'
   do op=1, NOPS
      print '(A,I0,A,F9.3,A,F8.1,A)', '| '//trim(OP_NAME(op))//' | ', OP_FIELDS(op), ' | ', t(op), ' | ', &
            real(OP_FIELDS(op), R8P) * 8._R8P * cells / (t(op) * 1.e-3_R8P) / 1.e9_R8P, ' |'
   enddo
   print '(A)', '| pair | access ratio | time ratio | time / access |'
   print '(A)', '|---|---|---|---|'
   call pair_out('SELAH + dq / ADAM FLUME', 4, 2)
   call pair_out('SELAH / ADAM', 3, 1)
   call pair_out('SELAH kernels / SELAH stepper, + dq', 8, 4)
   call pair_out('SELAH stage 1 on q / SELAH + dq, kernels', 5, 8)
   call pair_out('Ketcheson / Butcher SSP(10,4)', 7, 6)
   call pair_out('stage 1 stage state / copy', 9, 10)
   ! equal simulated time, dt proportional to the SSP coefficient: steps per unit time scale as 1/C
   ratio = (t(7) / bt104%ssp_coefficient) / (t(3) / bt54%ssp_coefficient)
   print '(A,F6.3,A,F6.3,A,F6.3)', 'equal simulated time, SSP(10,4) Ketcheson / SSP(5,4) Butcher: integrator ', ratio, &
         ' (accesses ', (39._R8P / bt104%ssp_coefficient) / (27._R8P / bt54%ssp_coefficient), &
         '), residual evaluations ', (10._R8P / bt104%ssp_coefficient) / (5._R8P / bt54%ssp_coefficient)
enddo

contains
   subroutine pair_out(name, a, b)
   !< One row: access ratio, time ratio, and their quotient (1 when the time follows the accesses).
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: a, b
   real(R8P)                :: ra, rt

   ra = real(OP_FIELDS(a), R8P) / OP_FIELDS(b) ; rt = t(a) / t(b)
   print '(A,F6.3,A,F6.3,A,F6.3,A)', '| '//name//' | ', ra, ' | ', rt, ' | ', rt / ra, ' |'
   endsubroutine pair_out

   subroutine field(a_d, v)
   !< A device field, initialised to v.
   real(R8P), pointer, intent(inout) :: a_d(:,:,:,:,:)
   real(R8P),          intent(in)    :: v
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NV], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P], init_value=v, ierr=ierr)
   if (ierr /= 0) error stop 'selah_bench_stepper: dev_alloc failed'
   endsubroutine field

   subroutine registers(a_d, n)
   !< n device registers with the field's bounds, initialised to 0.25.
   real(R8P), pointer, intent(inout) :: a_d(:,:,:,:,:,:)
   integer(I4P),       intent(in)    :: n
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NV, n], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P, 1_I4P], init_value=0.25_R8P, ierr=ierr)
   if (ierr /= 0) error stop 'selah_bench_stepper: dev_alloc failed'
   endsubroutine registers

   subroutine drive(st, with_dq)
   !< One step through the stepper, as a solver drives it; the caller's residual is left in the register.
   type(stepper), intent(inout) :: st
   logical,       intent(in)    :: with_dq
   integer(I4P)                 :: k

   call st%begin_step(bx, q_d)
   do k=1, st%stages()
      call st%begin_stage(k, bx, DT, q_d)
      call st%end_stage(k, bx, DT, q_d)
   enddo
   if (with_dq) then
      call st%finish_step(bx, DT, q_d, dq=dq_d)
   else
      call st%finish_step(bx, DT, q_d)
   endif
   endsubroutine drive

   subroutine launch(op)
   !< One step of operation op.
   integer(I4P), intent(in) :: op
   integer(I4P)             :: s

   select case(op)
   case(1, 2) ! ADAM: q copied into every stage, the stages in place, the residuals assigned, the update
      call rk_initialize_stages_dev(NI, NJ, NK, NGC, NV, nb, q_d, q_rk_d)
      do s=1, 5
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, nb, s, DT, alph_d, q_rk_gpu=q_rk_d)
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, nb, s, q_gpu=r_d, q_rk_gpu=q_rk_d)
      enddo
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, DT, beta_d, q_rk_gpu=q_rk_d, q_gpu=q_d)
      if (op == 2) call compute_rk_ssp_residual_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, beta_d, q_rk_d, dq_d)
   case(3) ; call drive(st54, .false.)
   case(4) ; call drive(st54, .true.)
   case(5) ! stage 1 evaluated on q itself: stage states 2..5 from q, then the update with dq
      do s=2, 5
         call selah_stage_state(bx, s-1, DT, at_d(:,s), q_d, reg54_d(:,:,:,:,:,2:), u_d)
      enddo
      call selah_step_update(bx, 5_I4P, DT, beta_d, reg54_d(:,:,:,:,:,2:), q_d, dq=dq_d)
   case(6) ; call drive(st104, .false.)
   case(7) ; call drive(stk, .false.)
   case(9) ; call selah_stage_state(bx, 0_I4P, DT, at_d(:,1), q_d, reg54_d(:,:,:,:,:,2:), u_d)
   case(10); call selah_copy(bx, q_d, u_d)
   case(8) ! op 4 without the stepper: stage states 1..5 from q, then the update with dq
      do s=1, 5
         call selah_stage_state(bx, s-1, DT, at_d(:,s), q_d, reg54_d(:,:,:,:,:,2:), u_d)
      enddo
      call selah_step_update(bx, 5_I4P, DT, beta_d, reg54_d(:,:,:,:,:,2:), q_d, dq=dq_d)
   endselect
   endsubroutine launch

   function timed(op) result(ms)
   !< Median over TRIALS of the time per launch [ms], after WARM launches.
   integer(I4P), intent(in) :: op
   real(R8P)                :: ms
   real(R8P)                :: tt(TRIALS)
   integer(I8P)             :: c0, c1, rate
   integer(I4P)             :: k, l

   do k=1, WARM
      call launch(op)
   enddo
   !$acc wait
   do k=1, TRIALS
      call system_clock(c0, rate)
      do l=1, LAUNCHES
         call launch(op)
      enddo
      !$acc wait
      call system_clock(c1)
      tt(k) = real(c1 - c0, R8P) / real(rate, R8P) / LAUNCHES * 1.e3_R8P
   enddo
   ms = median(tt)
   endfunction timed

   pure function median(x) result(m)
   real(R8P), intent(in) :: x(:)
   real(R8P)             :: m, s(size(x)), tmp
   integer(I4P)          :: i, j

   s = x
   do i=2, size(s)
      tmp = s(i) ; j = i - 1
      do while (j >= 1)
         if (s(j) <= tmp) exit
         s(j+1) = s(j) ; j = j - 1
      enddo
      s(j+1) = tmp
   enddo
   m = s((size(s)+1)/2)
   endfunction median
endprogram selah_bench_stepper
