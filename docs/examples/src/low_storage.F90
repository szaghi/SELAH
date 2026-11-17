!run low_storage low_storage
program low_storage
!< Tutorial, chapter 5: the same advection in the Butcher form and in a low-storage form of each scheme.
use selah
implicit none
integer(I4P), parameter :: NI=16, NB=4, NSTEPS=16
real(R_P),    parameter :: PI=acos(-1._R_P), T_END=0.25_R_P
character(14)           :: forms(3) = [character(14) :: FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_KETCHESON_2R]
character(7)            :: schemes(3) = [character(7) :: SCHEME_SSP_33, SCHEME_LSRK_54, SCHEME_SSP_104]
real(R_P), allocatable  :: q_b(:,:), q_l(:,:)
integer(I4P)            :: s

print '(A)', 'scheme   low-storage form   registers (Butcher, low-storage)   max |q_butcher - q_low-storage|'
do s=1, 3
   q_b = advect(trim(schemes(s)), FORM_BUTCHER)
   q_l = advect(trim(schemes(s)), trim(forms(s)))
   print '(A8,1X,A14,I12,I8,ES28.2)', schemes(s), forms(s), registers_needed(trim(schemes(s)), FORM_BUTCHER), &
         registers_needed(trim(schemes(s)), trim(forms(s))), maxval(abs(q_b - q_l))
enddo

contains
   function advect(scheme, form) result(q_final)
   !< NSTEPS steps of the periodic advection to T_END in the given storage form; the final solution in the cells.
   character(*), intent(in)       :: scheme, form
   real(R_P), allocatable         :: q_final(:,:)
   real(R_P), allocatable, target :: q(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                      :: dt, dx
   integer(I4P)                   :: step, k, b, i

   dx = 1._R_P / (NI * NB) ; dt = T_END / NSTEPS
   allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
   do b=1, NB
      do i=1, NI
         q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
      enddo
   enddo
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
   call st%initialize(scheme, form, lb=lbound(q), ub=ubound(q))
   do step=1, NSTEPS
      !region step
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)
         if (st%stage_state_is_q(k)) then                ! the low-storage forms keep the stage state in q itself
            u => q
         else                                            ! the Butcher form keeps it in a register
            u => st%stage_state()
         endif
         r => st%rhs_register()
         call fill_ghosts(u)
         call residual(u, r, dx)
         call st%end_stage(k, bx, dt, q)                 ! the low-storage forms advance q here
      enddo
      call st%finish_step(bx, dt, q)                     ! the Butcher form advances q here
      !endregion step
   enddo
   q_final = q(1,1:NI,1,1,:)                             ! (i, b): the cells of every block
   call st%destroy
   endfunction advect

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
endprogram low_storage
