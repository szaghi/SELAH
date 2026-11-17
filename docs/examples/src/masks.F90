!run masks masks
program masks
!< Tutorial, chapter 7: a mask freezes solid cells; the two-register form of ssp-33 on the periodic advection.
use selah
implicit none
integer(I4P), parameter        :: NI=16, NB=4, NSTEPS=8
real(R_P),    parameter        :: PI=acos(-1._R_P)
real(R_P), allocatable, target :: q(:,:,:,:,:)
real(R_P), allocatable         :: q0(:,:,:,:,:), mask(:,:,:,:,:)
type(stepper)                  :: st
type(box)                      :: bx
real(R_P), pointer             :: r(:,:,:,:,:)
real(R_P)                      :: dt, dx
integer(I4P)                   :: step, k, b, i

dx = 1._R_P / (NI * NB) ; dt = 0.25_R_P / NSTEPS
allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
do b=1, NB
   do i=1, NI
      q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
   enddo
enddo
q0 = q
!region mask
! the mask has the field's bounds and one component; a cell is active (fluid) where mask < 0
allocate(mask(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
mask = -1._R_P
mask(:,:,:,:,3) = 1._R_P                                          ! block 3 is solid
!endregion mask
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=lbound(q), ub=ubound(q))
do step=1, NSTEPS
   !region step
   call st%begin_step(bx, q, mask=mask)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, q, mask=mask)
      r => st%rhs_register()
      call fill_ghosts(q)
      call residual(q, r, dx)                                     ! computed everywhere: the stepper masks the update
      call st%end_stage(k, bx, dt, q, mask=mask)
   enddo
   call st%finish_step(bx, dt, q, mask=mask)
   !endregion step
enddo
!region result
do b=1, NB
   print '(A,I0,A,I2,A,I0,A)', 'block ', b, ': ', count(q(1,1:NI,1,1,b) /= q0(1,1:NI,1,1,b)), ' of ', NI, &
         ' cells changed'//trim(merge(' (solid)', '        ', b == 3))
enddo
!endregion result
call st%destroy

contains
   subroutine fill_ghosts(u)
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   integer(I4P)             :: b
   do b=1, NB
      u(1,0,1,1,b)    = u(1,NI,1,1,1 + modulo(b-2, NB))
      u(1,NI+1,1,1,b) = u(1,1,1,1,1 + modulo(b, NB))
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r, dx)
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   integer(I4P)             :: b, i
   do b=1, NB
      do i=1, NI
         r(1,i,1,1,b) = -(u(1,i+1,1,1,b) - u(1,i-1,1,1,b)) / (2._R_P*dx)
      enddo
   enddo
   endsubroutine residual
endprogram masks
