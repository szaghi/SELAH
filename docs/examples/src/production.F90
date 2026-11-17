!run production production
program production
!< Tutorial, chapter 10: a production loop. The registers come from the caller's pool, the first stage is evaluated on
!< q itself (stages_on_q), and the step residual dq is fused in the update and monitored.
use selah
implicit none
integer(I4P), parameter        :: NI=16, NB=4, NSTEPS=6
real(R_P),    parameter        :: PI=acos(-1._R_P)
real(R_P), allocatable, target :: q(:,:,:,:,:), dq(:,:,:,:,:), q_old(:,:,:,:,:)
real(R_P), pointer             :: pool(:,:,:,:,:,:)=>null(), u(:,:,:,:,:), r(:,:,:,:,:)
type(stepper)                  :: st
type(box)                      :: bx
real(R_P)                      :: dt, dx
integer(I4P)                   :: step, k, b, i

dx = 1._R_P / (NI * NB) ; dt = 0.5_R_P * dx
allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB), dq(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
do b=1, NB
   do i=1, NI
      q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
   enddo
enddo
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
!region setup
! the caller's pool: the field's bounds, the slots numbered from 1, as many as the form needs
allocate(pool(1:1, 0:NI+1, 0:2, 0:2, 1:NB, 1:registers_needed(SCHEME_SSP_54, FORM_BUTCHER)))
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q), registers=pool, stages_on_q=.true.)
!endregion setup
print '(A)', 'step   max |dq|    max |q_new - q_old - dt dq|'
do step=1, NSTEPS
   q_old = q
   !region step
   call st%begin_step(bx, q)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, q)
      if (st%stage_state_is_q(k)) then                    ! stage 1: its state is q, no copy was made
         u => q
      else
         u => st%stage_state()
      endif
      r => st%rhs_register()
      call fill_ghosts(u)
      call residual(u, r, dx)
      call st%end_stage(k, bx, dt, q)
   enddo
   call st%finish_step(bx, dt, q, dq=dq)                  ! the update writes dq = sum_k b_k F_k in the same pass
   !endregion step
   print '(I4,ES12.3,ES14.2)', step, maxval(abs(dq(1,1:NI,1,1,:))), &
         maxval(abs(q(1,1:NI,1,1,:) - q_old(1,1:NI,1,1,:) - dt * dq(1,1:NI,1,1,:)))
enddo
call st%destroy
deallocate(pool)

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
endprogram production
