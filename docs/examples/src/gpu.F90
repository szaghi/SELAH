!run gpu gpu
#include "fundal.H"

program gpu
!< Tutorial, chapter 9: chapter 4 on device memory. The field is a FUNDAL array (device memory in an OpenACC build,
!< host memory otherwise), the caller's loops are FUNDAL's portable directives, and the stepper's registers and
!< kernels live on the device with it: nothing crosses the bus during a step.
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
implicit none
integer(I4P), parameter :: NI=16, NB=4, NSTEPS=16
real(R_P),    parameter :: PI=acos(-1._R_P), T_END=0.25_R_P
real(R_P), pointer      :: q_d(:,:,:,:,:)=>null(), u_d(:,:,:,:,:), r_d(:,:,:,:,:)
real(R_P), allocatable  :: q(:,:,:,:,:)
type(stepper)           :: st
type(box)               :: bx
real(R_P)               :: dt, dx, c, err
integer(I4P)            :: step, k, b, i, ierr

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
#endif
dx = 1._R_P / (NI * NB) ; dt = T_END / NSTEPS
!region field
allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
q = 0._R_P
do b=1, NB
   do i=1, NI
      q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
   enddo
enddo
call dev_alloc(fptr_dev=q_d, ubounds=ubound(q), lbounds=lbound(q), ierr=ierr)   ! the field on the device
call dev_memcpy_to_device(dst=q_d, src=q)
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))      ! registers through FUNDAL: device
!endregion field
do step=1, NSTEPS
   call st%begin_step(bx, q_d)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, q_d)
      u_d => st%stage_state() ; r_d => st%rhs_register()
      call fill_ghosts(u_d)
      call residual(u_d, r_d, dx)
      call st%end_stage(k, bx, dt, q_d)
   enddo
   call st%finish_step(bx, dt, q_d)
enddo
call dev_memcpy_from_device(dst=q, src=q_d)                                    ! back to the host, once
c = sin(2._R_P*PI*dx) / (2._R_P*PI*dx)
err = 0._R_P
do b=1, NB
   do i=1, NI
      err = max(err, abs(q(1,i,1,1,b) - sin(2._R_P*PI*(((b-1)*NI + i - 0.5_R_P)*dx - c*T_END))))
   enddo
enddo
print '(A,ES10.3)', 'error at t = 0.25: ', err
call st%destroy
call dev_free(q_d)

contains
   !region kernels
   subroutine fill_ghosts(u)
   !< Periodic ghost cells, one device thread per block.
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   integer(I4P)             :: b

   !$acc parallel loop independent DEVICEVAR(u)
   !$omp OMPLOOP DEVICEPTR(u)
   do b=1, NB
      u(1,0,1,1,b)    = u(1,NI,1,1,1 + modulo(b-2, NB))
      u(1,NI+1,1,1,b) = u(1,1,1,1,1 + modulo(b, NB))
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r, dx)
   !< The central residual, one device thread per cell.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   integer(I4P)             :: b, i

   !$acc parallel loop independent collapse(2) DEVICEVAR(u, r)
   !$omp OMPLOOP collapse(2) DEVICEPTR(u, r)
   do b=1, NB
      do i=1, NI
         r(1,i,1,1,b) = -(u(1,i+1,1,1,b) - u(1,i-1,1,1,b)) / (2._R_P*dx)
      enddo
   enddo
   endsubroutine residual
   !endregion kernels
endprogram gpu
