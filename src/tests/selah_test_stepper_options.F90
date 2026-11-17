!< SELAH test: the stepper's optional arguments and queries, on the host or the device.

#include "fundal.H"

program selah_test_stepper_options
!< SELAH test: the stepper's optional arguments and queries, on the host or the device.
!<
!< The caller is the periodic 1-D advection of selah_test_stepper (NB blocks of NI cells, one ghost cell, central
!< residual), here driven through the arguments that test does not use:
!<
!<+ **stage_time.** stage_time(k) is the abscissa c_k of the tableau in every Runge–Kutta form, 0 for a splitting.
!<+ **registers=.** The caller's registers, filled with NaN, give the solution of the stepper's own (zeroed) registers
!<  bit for bit in every form of lsrk-54, ssp-54, ssp-33 and ssp-104: no form reads a register before writing it (the
!<  2N form's first stage, A_1 = 0, included).
!<+ **finish_step(dq=).** The step residual is (q^{n+1} - q^n)/dt to round-off.
!<+ **mask=.** With half of the blocks solid (mask > 0), the solid cells of q keep their values bit for bit through
!<  whole steps, in every form, and the fluid cells stay finite.
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: NI=16, NB=4, NGC=1                 !< Cells per block, blocks, ghosts.
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P),    parameter :: DX = 1._R_P / real(NI*NB, R_P)     !< Cell size.
integer(I4P), parameter :: LB(5) = [1, 0, 0, 0, 1], UB(5) = [1, NI+1, 2, 2, NB] !< Field bounds.
integer(I4P)            :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
call stage_times
call caller_registers(SCHEME_LSRK_54, FORM_WILLIAMSON_2N)
call caller_registers(SCHEME_LSRK_54, FORM_BUTCHER)
call caller_registers(SCHEME_SSP_54, FORM_BUTCHER)
call caller_registers(SCHEME_SSP_33, FORM_SHU_OSHER_2R)
call caller_registers(SCHEME_SSP_104, FORM_KETCHESON_2R)
do s=1, size(RK_SCHEMES)
   call step_residual(trim(RK_SCHEMES(s)))
enddo
call masked(SCHEME_SSP_54, FORM_BUTCHER)
call masked(SCHEME_SSP_33, FORM_SHU_OSHER_2R)
call masked(SCHEME_LSRK_54, FORM_WILLIAMSON_2N)
call masked(SCHEME_SSP_104, FORM_KETCHESON_2R)
if (failures > 0) error stop 'selah_test_stepper_options: failed'
print '(A)', 'selah_test_stepper_options: ok'

contains
   subroutine stage_times
   !< stage_time(k) = c_k in every Runge–Kutta form; 0 for a splitting.
   type(stepper)         :: st
   type(butcher_tableau) :: bt
   integer(I4P)          :: k
   logical               :: ok

   ok = .true.
   do s=1, size(RK_SCHEMES)
      bt = butcher(trim(RK_SCHEMES(s)))
      call st%initialize(trim(RK_SCHEMES(s)), FORM_BUTCHER, lb=LB, ub=UB)
      do k=1, st%stages()
         if (st%stage_time(k) /= bt%c(k)) ok = .false.
      enddo
      call st%destroy
   enddo
   call check(ok, 'stage_time(k) = c_k for every scheme')
   call st%initialize(SCHEME_STRANG, FORM_SPLITTING, lb=LB, ub=UB, partitions=reshape([1, 1, 1, 1], [2, 2]))
   ok = .true.
   do k=1, st%stages()
      if (st%stage_time(k) /= 0._R_P) ok = .false.
   enddo
   call st%destroy
   call check(ok, 'stage_time(k) = 0 for a splitting')
   endsubroutine stage_times

   subroutine caller_registers(name, form)
   !< The caller's NaN-filled registers against the stepper's own: the same solution, bit for bit.
   character(*), intent(in) :: name, form
   real(R_P), pointer       :: regs_d(:,:,:,:,:,:)=>null()
   real(R_P), allocatable   :: q_own(:,:,:,:,:), q_caller(:,:,:,:,:)
   integer(I4P)             :: ierr

   call dev_alloc(fptr_dev=regs_d, ubounds=[UB, registers_needed(name, form)], lbounds=[LB, 1_I4P], &
                  init_value=ieee_value(1._R_P, ieee_quiet_nan), ierr=ierr)
   call advect(name, form, 8_I4P, q_own)
   call advect(name, form, 8_I4P, q_caller, regs_d=regs_d)
   call check(all(q_own(:,1:NI,1,1,:) == q_caller(:,1:NI,1,1,:)), &
              name//' '//form//': NaN-filled caller registers give the solution of the own registers, bitwise')
   call dev_free(regs_d)
   endsubroutine caller_registers

   subroutine step_residual(name)
   !< finish_step(dq=) in the Butcher form: dq = (q^{n+1} - q^n)/dt to round-off, at every step.
   character(*), intent(in) :: name
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q_d(:,:,:,:,:)=>null(), dq_d(:,:,:,:,:)=>null()
   real(R_P), allocatable   :: q0(:,:,:,:,:), q1(:,:,:,:,:), dq(:,:,:,:,:)
   real(R_P)                :: dt, worst
   integer(I4P)             :: step, ierr

   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, NGC, NB)
   call field(q_d)
   call dev_alloc(fptr_dev=dq_d, ubounds=UB, lbounds=LB, init_value=0._R_P, ierr=ierr)
   allocate(q0(1:1,0:NI+1,0:2,0:2,1:NB), q1(1:1,0:NI+1,0:2,0:2,1:NB), dq(1:1,0:NI+1,0:2,0:2,1:NB))
   call st%initialize(name, FORM_BUTCHER, lb=LB, ub=UB)
   dt = 0.25_R_P / 8._R_P
   worst = 0._R_P
   do step=1, 4
      call dev_memcpy_from_device(dst=q0, src=q_d)
      call one_step(st, bx, dt, q_d, dq_d=dq_d)
      call dev_memcpy_from_device(dst=q1, src=q_d)
      call dev_memcpy_from_device(dst=dq, src=dq_d)
      worst = max(worst, maxval(abs((q1(:,1:NI,1,1,:) - q0(:,1:NI,1,1,:)) - dt * dq(:,1:NI,1,1,:))))
   enddo
   call check(worst <= 8._R_P * epsilon(1._R_P), name//' butcher: dq = (q^{n+1} - q^n)/dt to round-off')
   call st%destroy
   call dev_free(q_d) ; call dev_free(dq_d)
   endsubroutine step_residual

   subroutine masked(name, form)
   !< Blocks 2 and 4 solid: their cells of q keep their values bit for bit; the fluid cells stay finite.
   character(*), intent(in) :: name, form
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q_d(:,:,:,:,:)=>null(), m_d(:,:,:,:,:)=>null()
   real(R_P), allocatable   :: q0(:,:,:,:,:), q1(:,:,:,:,:), m(:,:,:,:,:)
   integer(I4P)             :: step, b, ierr

   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, NGC, NB)
   call field(q_d)
   allocate(m(1:1,0:NI+1,0:2,0:2,1:NB), q0(1:1,0:NI+1,0:2,0:2,1:NB), q1(1:1,0:NI+1,0:2,0:2,1:NB))
   do b=1, NB
      m(:,:,:,:,b) = merge(1._R_P, -1._R_P, mod(b, 2) == 0)
   enddo
   call dev_alloc(fptr_dev=m_d, ubounds=UB, lbounds=LB, ierr=ierr)
   call dev_memcpy_to_device(dst=m_d, src=m)
   call dev_memcpy_from_device(dst=q0, src=q_d)
   call st%initialize(name, form, lb=LB, ub=UB)
   do step=1, 4
      call one_step(st, bx, 0.25_R_P / 8._R_P, q_d, mask_d=m_d)
   enddo
   call dev_memcpy_from_device(dst=q1, src=q_d)
   call check(all(q1(:,1:NI,1,1,2:NB:2) == q0(:,1:NI,1,1,2:NB:2)), name//' '//form//': solid cells of q untouched')
   call check(all(ieee_is_finite(q1(:,1:NI,1,1,1:NB:2))), name//' '//form//': fluid cells finite')
   call st%destroy
   call dev_free(q_d) ; call dev_free(m_d)
   endsubroutine masked

   subroutine advect(name, form, nsteps, q_final, regs_d)
   !< nsteps periodic steps to t = 0.25; the final solution on the host.
   character(*),           intent(in)           :: name, form
   integer(I4P),           intent(in)           :: nsteps
   real(R_P), allocatable, intent(out)          :: q_final(:,:,:,:,:)
   real(R_P), pointer,     intent(in), optional :: regs_d(:,:,:,:,:,:)
   type(stepper)                                :: st
   type(box)                                    :: bx
   real(R_P), pointer                           :: q_d(:,:,:,:,:)=>null()
   integer(I4P)                                 :: step

   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, NGC, NB)
   call field(q_d)
   call st%initialize(name, form, lb=LB, ub=UB, registers=regs_d)
   do step=1, nsteps
      call one_step(st, bx, 0.25_R_P / real(nsteps, R_P), q_d)
   enddo
   allocate(q_final(1:1,0:NI+1,0:2,0:2,1:NB))
   call dev_memcpy_from_device(dst=q_final, src=q_d)
   call st%destroy
   call dev_free(q_d)
   endsubroutine advect

   subroutine one_step(st, bx, dt, q_d, mask_d, dq_d)
   !< One step as a solver drives it: ghost fill and residual of each stage state, optional mask and dq.
   type(stepper),      intent(inout)        :: st
   type(box),          intent(in)           :: bx
   real(R_P),          intent(in)           :: dt
   real(R_P), pointer, intent(inout)        :: q_d(:,:,:,:,:)
   real(R_P), pointer, intent(in), optional :: mask_d(:,:,:,:,:), dq_d(:,:,:,:,:)
   real(R_P), pointer                       :: u_d(:,:,:,:,:), r_d(:,:,:,:,:)
   integer(I4P)                             :: k

   if (present(mask_d)) then
      call st%begin_step(bx, q_d, mask=mask_d)
   else
      call st%begin_step(bx, q_d)
   endif
   do k=1, st%stages()
      if (present(mask_d)) then
         call st%begin_stage(k, bx, dt, q_d, mask=mask_d)
      else
         call st%begin_stage(k, bx, dt, q_d)
      endif
      if (st%stage_state_is_q(k)) then
         u_d => q_d
      else
         u_d => st%stage_state()
      endif
      r_d => st%rhs_register()
      call fill_ghosts(u_d)
      call residual(u_d, r_d)
      if (present(mask_d)) then
         call st%end_stage(k, bx, dt, q_d, mask=mask_d)
      else
         call st%end_stage(k, bx, dt, q_d)
      endif
   enddo
   if (present(dq_d)) then
      call st%finish_step(bx, dt, q_d, dq=dq_d)
   elseif (present(mask_d)) then
      call st%finish_step(bx, dt, q_d, mask=mask_d)
   else
      call st%finish_step(bx, dt, q_d)
   endif
   endsubroutine one_step

   subroutine field(q_d)
   !< A device field holding sin(2 pi x) in the cells, 0 in the ghost cells.
   real(R_P), pointer, intent(inout) :: q_d(:,:,:,:,:)
   real(R_P)                         :: q(1:1,0:NI+1,0:2,0:2,1:NB), x
   integer(I4P)                      :: b, i, ierr

   q = 0._R_P
   do b=1, NB
      do i=1, NI
         x = (real((b-1)*NI + i, R_P) - 0.5_R_P) * DX
         q(1,i,1,1,b) = sin(2._R_P*PI*x)
      enddo
   enddo
   call dev_alloc(fptr_dev=q_d, ubounds=UB, lbounds=LB, ierr=ierr)
   call dev_memcpy_to_device(dst=q_d, src=q)
   endsubroutine field

   subroutine fill_ghosts(u)
   !< Periodic ghost fill along i.
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   integer(I4P)             :: b, bl, br

   !$acc parallel loop independent DEVICEVAR(u) private(bl, br)
   !$omp OMPLOOP DEVICEPTR(u) private(bl, br)
   do b=1, NB
      bl = b - 1 ; if (bl < 1) bl = NB
      br = b + 1 ; if (br > NB) br = 1
      u(1,0,1,1,b) = u(1,NI,1,1,bl)
      u(1,NI+1,1,1,b) = u(1,1,1,1,br)
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r)
   !< Flux form, central face fluxes (u_L + u_R)/2.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   integer(I4P)             :: b, i

   !$acc parallel loop independent collapse(2) DEVICEVAR(u, r)
   !$omp OMPLOOP collapse(2) DEVICEPTR(u, r)
   do b=1, NB
      do i=1, NI
         r(1,i,1,1,b) = -(0.5_R_P*(u(1,i,1,1,b) + u(1,i+1,1,1,b)) - 0.5_R_P*(u(1,i-1,1,1,b) + u(1,i,1,1,b))) / DX
      enddo
   enddo
   endsubroutine residual
endprogram selah_test_stepper_options
