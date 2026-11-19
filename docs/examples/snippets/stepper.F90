program stepper_tour
!< Tutorial, chapter 4: the stepper drives a 1-D periodic advection; the caller owns ghosts and residual.
use selah
use gideon, only : gideon_fd_d1_c2   ! the space operator: GIDEON's centred first derivative
implicit none
integer(I4P), parameter :: NI=16, NB=4                       ! 4 blocks of 16 cells, 1 ghost cell
real(R_P),    parameter :: PI=acos(-1._R_P), T_END=0.25_R_P
real(R_P)               :: err(3)
integer(I4P)            :: k

do k=1, 3
   err(k) = advect(SCHEME_SSP_54, 8 * 2**(k-1))
   print '(A,I3,A,ES10.3)', 'steps ', 8 * 2**(k-1), ': error ', err(k)
enddo
print '(A,2F6.2)', 'observed orders:', log(err(1:2) / err(2:3)) / log(2._R_P)

contains
   function advect(scheme, nsteps) result(err)
   !< nsteps steps of u_t + u_x = 0 to T_END; the error against the exact semi-discrete solution.
   character(*), intent(in) :: scheme
   integer(I4P), intent(in) :: nsteps
   real(R_P)                :: err
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), allocatable   :: q(:,:,:,:,:)              ! the solution: the caller's array
   real(R_P), pointer       :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                :: dt, dx, c
   integer(I4P)             :: step, k, b, i

   dx = 1._R_P / (NI * NB) ; dt = T_END / nsteps
   allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))              ! (v,i,j,k,b): 1 component, NI x 1 x 1 cells, 1 ghost
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
   call st%initialize(scheme, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
   do b=1, NB
      do i=1, NI
         q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
      enddo
   enddo
   do step=1, nsteps
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)               ! the stepper builds the state of stage k
         u => st%stage_state()                           ! Butcher form: its own register
         r => st%rhs_register()                          ! where the residual of stage k goes
         call fill_ghosts(u)                             ! the caller: ghost cells of the stage state
         call residual(u, r, dx)                         ! the caller: r = -du/dx
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)                     ! q^{n+1}
   enddo
   c = sin(2._R_P*PI*dx) / (2._R_P*PI*dx)                ! the semi-discrete phase speed of the central residual
   err = 0._R_P
   do b=1, NB
      do i=1, NI
         err = max(err, abs(q(1,i,1,1,b) - sin(2._R_P*PI*(((b-1)*NI + i - 0.5_R_P)*dx - c*T_END))))
      enddo
   enddo
   call st%destroy
   endfunction advect

   subroutine fill_ghosts(u)
   !< Periodic ghost cells along i: each block's ghosts from its neighbours.
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   integer(I4P)             :: b
   do b=1, NB
      u(1,0,1,1,b)    = u(1,NI,1,1,1 + modulo(b-2, NB))
      u(1,NI+1,1,1,b) = u(1,1,1,1,1 + modulo(b, NB))
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r, dx)
   !< Central flux form: r_i = -((u_i + u_{i+1})/2 - (u_{i-1} + u_i)/2)/dx.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   integer(I4P)             :: b, i
   real(R_P)                :: d
   do b=1, NB
      do i=1, NI
         call gideon_fd_d1_c2([u(1,i-1,1,1,b), u(1,i,1,1,b), u(1,i+1,1,1,b)], dx, d) ; r(1,i,1,1,b) = -d
      enddo
   enddo
   endsubroutine residual
endprogram stepper_tour
