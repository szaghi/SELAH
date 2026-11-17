!< SELAH benchmark: the kernels against ADAM's FNL kernels on the device, and the bandwidth they reach.

#include "fundal.H"

program selah_bench_kernels
!< SELAH benchmark: the kernels against ADAM's FNL kernels on the device, and the bandwidth they reach.
!<
!< A FLUME-like device field (layout V_LAST): 5 components, 16^3 interior cells with 3 ghosts, NB blocks (default 256,
!< ~109 MB per field). Each operation of SELAH is timed against ADAM's kernel doing the same work (src/tests/
!< adam_fnl_rk_kernels.F90, verbatim): the stage state of stage 5 of SSP(5,4) (4 registers), the step update
!< (5 registers), the copy of a residual into a register, a two-register stage, and the masked stage state. A launch is
!< synchronous (`!$acc wait` before the clock is read: system_clock, never cpu_time); 3 warm-up launches, then the
!< median of 7 trials of 20 launches. The effective bandwidth counts the fields read and written over the interior
!< cells only, and is compared with a STREAM triad (a = b + s c) over the same number of bytes.
!<
!< An unreported pass of every operation comes first (the first operations timed after the allocations are slow for a
!< while on the A30). The table is then printed PASSES times (default 3), SELAH timed first on odd passes and ADAM first
!< on even ones, so that an effect of the order shows up as a ratio that changes from pass to pass.
!<
!< Usage: selah_bench_kernels [NB [PASSES]]
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_assign_to_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use penf, only : I4P, I8P, R8P
use selah
use adam_fnl_rk_kernels
implicit none
integer(I4P), parameter :: NV=5, NI=16, NJ=16, NK=16, NGC=3, NS=1, WARM=3, TRIALS=7, LAUNCHES=20
real(R8P),    parameter :: DT=1.e-3_R8P
integer(I4P)            :: nb
real(R8P), pointer      :: q_d(:,:,:,:,:)=>null(), u_d(:,:,:,:,:)=>null(), r_d(:,:,:,:,:)=>null()
real(R8P), pointer      :: phi_d(:,:,:,:,:)=>null(), q_rk_d(:,:,:,:,:,:)=>null(), regs_d(:,:,:,:,:,:)=>null()
real(R8P), pointer      :: alph_d(:,:)=>null(), beta_d(:)=>null(), row_d(:)=>null(), at_d(:,:)=>null()
real(R8P), pointer      :: ta(:)=>null(), tb(:)=>null(), tc(:)=>null()
type(butcher_tableau)   :: bt
type(box)               :: bx
real(R8P), allocatable  :: row(:)
real(R8P)               :: cells, triad_gbs, t_selah, t_adam
character(32)           :: arg
integer(I4P)            :: ierr, passes, pass

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
bt = butcher(SCHEME_SSP_54)
bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, nb)
cells = real(NI, R8P) * NJ * NK * nb * NV ! updated values per field
call field(q_d, 1._R8P) ; call field(u_d, 0._R8P) ; call field(r_d, 0.5_R8P) ; call mask_field(phi_d)
call dev_alloc(fptr_dev=q_rk_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NV, bt%stages], &
               lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P, 1_I4P], init_value=0.25_R8P, ierr=ierr)
call dev_alloc(fptr_dev=regs_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NV, bt%stages], &
               lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P, 1_I4P], init_value=0.25_R8P, ierr=ierr)
call dev_assign_to_device(dst=alph_d, src=bt%a) ; call dev_assign_to_device(dst=beta_d, src=bt%b)
row = bt%a(5,:) ; call dev_assign_to_device(dst=row_d, src=row)
call dev_assign_to_device(dst=at_d, src=transpose(bt%a)) ! column s = row s of A, contiguous
triad_gbs = triad()

print '(A,I0,A,F7.1,A)', 'blocks ', nb, ', field ', 8._R8P*NV*(NI+2*NGC)*(NJ+2*NGC)*(NK+2*NGC)*nb/1.e6_R8P, ' MB'
print '(A,F8.1,A)', 'STREAM triad: ', triad_gbs, ' GB/s'
! an unreported pass of every operation first: on the A30 the first operations timed after the allocations ran up to
! 30% slow for a pass (job 1368086 put stage state at 1.31 x ADAM, job 1368088 at 0.97 once warm), whatever the order
do ierr=1, 13
   t_selah = timed(ierr)
enddo
do pass=1, passes
   print '(A,I0,A)', 'pass ', pass, merge(' (SELAH timed first)', ' (ADAM timed first) ', mod(pass, 2) == 1)
   print '(A)', '| operation | fields | SELAH [ms] | ADAM [ms] | SELAH/ADAM | SELAH [GB/s] | of triad |'
   print '(A)', '|---|---|---|---|---|---|---|'
   call pair(1, 'stage state, stage 5 of 5', 6)
   call pair(3, 'step update, 5 registers', 7)
   call pair(5, 'copy residual to register', 2)
   call pair(7, 'two-register stage', 4)
   call pair(9, 'stage state, masked', 7)
   t_selah = timed(11) ; call row_out('step update + dq (fused)', 8, t_selah, -1._R8P)
   call pair(12, 'whole SSP(5,4) step (SELAH 37, ADAM 43 accesses)', 37)
enddo

contains
   subroutine pair(op, name, fields)
   !< Time SELAH (op) and ADAM (op+1), in the order of the pass, and print the row.
   integer(I4P), intent(in) :: op
   character(*), intent(in) :: name
   integer,      intent(in) :: fields

   if (mod(pass, 2) == 1) then
      t_selah = timed(op) ; t_adam = timed(op+1)
   else
      t_adam = timed(op+1) ; t_selah = timed(op)
   endif
   call row_out(name, fields, t_selah, t_adam)
   endsubroutine pair

   subroutine field(a_d, v)
   !< A device field of the pool, initialised to v.
   real(R8P), pointer, intent(inout) :: a_d(:,:,:,:,:)
   real(R8P),          intent(in)    :: v
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NV], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P], init_value=v, ierr=ierr)
   if (ierr /= 0) error stop 'selah_bench_kernels: dev_alloc failed'
   endsubroutine field

   subroutine mask_field(m_d)
   !< A device mask (one component): a checkerboard of fluid (-1) and solid (+1) blocks.
   real(R8P), pointer, intent(inout) :: m_d(:,:,:,:,:)
   real(R8P), allocatable            :: m(:,:,:,:,:)
   integer(I4P)                      :: b, ierr

   allocate(m(1:nb, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NS))
   do b=1, nb
      m(b,:,:,:,:) = merge(-1._R8P, 1._R8P, mod(b, 2) == 1)
   enddo
   call dev_alloc(fptr_dev=m_d, ubounds=[nb, NI+NGC, NJ+NGC, NK+NGC, NS], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=m_d, src=m)
   endsubroutine mask_field

   subroutine launch(op)
   !< One launch of operation op (SELAH: 1, 3, 5, 7, 9, 11, 12; ADAM: the next number).
   integer(I4P), intent(in) :: op
   integer(I4P)             :: s

   select case(op)
   case(1) ; call selah_stage_state(bx, 4_I4P, DT, row_d, q_d, regs_d, u_d)
   case(2) ; call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, DT, alph_d, q_rk_gpu=q_rk_d)
   case(3) ; call selah_step_update(bx, 5_I4P, DT, beta_d, regs_d, q_d)
   case(4) ; call rk_update_q_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, DT, beta_d, q_rk_gpu=q_rk_d, q_gpu=q_d)
   case(5) ; call selah_copy(bx, r_d, regs_d(:,:,:,:,:,1))
   case(6) ; call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, nb, 1_I4P, q_gpu=r_d, q_rk_gpu=q_rk_d)
   case(7) ; call selah_shu_osher_2r_stage(bx, DT, 0.75_R8P, 0.25_R8P, 0.25_R8P, u_d, r_d, q_d)
   case(8) ; call rk_compute_stage_ls_dev(NI, NJ, NK, NGC, NV, nb, DT, 0.75_R8P, 0.25_R8P, 0.25_R8P, &
                                          q_n_gpu=q_rk_d(:,:,:,:,:,1), dq_gpu=r_d, q_rk_gpu=q_d)
   case(9) ; call selah_stage_state(bx, 4_I4P, DT, row_d, q_d, regs_d, u_d, mask=phi_d)
   case(10); call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, DT, alph_d, phi_gpu=phi_d, q_rk_gpu=q_rk_d)
   case(11); call selah_step_update(bx, 5_I4P, DT, beta_d, regs_d, q_d, dq=u_d)
   case(12) ! a whole SSP(5,4) step, SELAH: the stage states from q, the residuals copied into the registers, the update
      do s=1, 5
         call selah_stage_state(bx, s-1, DT, at_d(:,s), q_d, regs_d, u_d)
         call selah_copy(bx, r_d, regs_d(:,:,:,:,:,s))
      enddo
      call selah_step_update(bx, 5_I4P, DT, beta_d, regs_d, q_d)
   case(13) ! the same step, ADAM: q copied into every stage, the stages in place, the residuals assigned, the update
      call rk_initialize_stages_dev(NI, NJ, NK, NGC, NV, nb, q_d, q_rk_d)
      do s=1, 5
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, nb, s, DT, alph_d, q_rk_gpu=q_rk_d)
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, nb, s, q_gpu=r_d, q_rk_gpu=q_rk_d)
      enddo
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, nb, 5_I4P, DT, beta_d, q_rk_gpu=q_rk_d, q_gpu=q_d)
   endselect
   endsubroutine launch

   function timed(op) result(ms)
   !< Median over TRIALS of the time per launch [ms], after WARM launches.
   integer(I4P), intent(in) :: op
   real(R8P)                :: ms
   real(R8P)                :: t(TRIALS)
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
      t(k) = real(c1 - c0, R8P) / real(rate, R8P) / LAUNCHES * 1.e3_R8P
   enddo
   ms = median(t)
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

   function triad() result(gbs)
   !< STREAM triad over 3 arrays of the size of one field: the bandwidth reference [GB/s].
   real(R8P)    :: gbs
   real(R8P)    :: t(TRIALS), s
   integer(I8P) :: n, i, c0, c1, rate
   integer(I4P) :: k, l, ierr

   n = int(NV, I8P) * (NI+2*NGC) * (NJ+2*NGC) * (NK+2*NGC) * nb
   call dev_alloc(fptr_dev=ta, ubounds=[int(n, I4P)], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=tb, ubounds=[int(n, I4P)], init_value=1._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=tc, ubounds=[int(n, I4P)], init_value=2._R8P, ierr=ierr)
   s = 3._R8P
   do k=1, TRIALS
      call system_clock(c0, rate)
      do l=1, LAUNCHES
         !$acc parallel loop independent DEVICEVAR(ta, tb, tc)
         do i=1, n
            ta(i) = tb(i) + s * tc(i)
         enddo
      enddo
      !$acc wait
      call system_clock(c1)
      t(k) = real(c1 - c0, R8P) / real(rate, R8P) / LAUNCHES
   enddo
   gbs = 3._R8P * 8._R8P * real(n, R8P) / median(t) / 1.e9_R8P
   call dev_free(ta) ; call dev_free(tb) ; call dev_free(tc)
   endfunction triad

   subroutine row_out(name, fields, ms_selah, ms_adam)
   !< One table row; fields = arrays read plus written per updated value.
   character(*), intent(in) :: name
   integer,      intent(in) :: fields
   real(R8P),    intent(in) :: ms_selah, ms_adam
   real(R8P)                :: gbs

   gbs = real(fields, R8P) * 8._R8P * cells / (ms_selah * 1.e-3_R8P) / 1.e9_R8P
   if (ms_adam > 0._R8P) then
      print '(A,I0,A,F9.3,A,F9.3,A,F6.3,A,F8.1,A,F5.2,A)', '| '//name//' | ', fields, ' | ', ms_selah, ' | ', ms_adam, &
            ' | ', ms_selah/ms_adam, ' | ', gbs, ' | ', gbs/triad_gbs, ' |'
   else
      print '(A,I0,A,F9.3,A,F8.1,A,F5.2,A)', '| '//name//' | ', fields, ' | ', ms_selah, ' | — | — | ', gbs, ' | ', &
            gbs/triad_gbs, ' |'
   endif
   endsubroutine row_out
endprogram selah_bench_kernels
