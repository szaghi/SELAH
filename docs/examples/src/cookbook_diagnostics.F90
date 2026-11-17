!run diagnostics-steady cookbook_diagnostics steady
!run diagnostics-contract cookbook_diagnostics contract
!run -s diagnostics-wrong-form cookbook_diagnostics wrong-form
program cookbook_diagnostics
!< Cookbook, diagnostics: the step residual as a convergence monitor; the contract checks a caller can make; what the
!< stepper does with a request it cannot honour.
use selah
implicit none
character(16) :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('steady')     ; call steady
case('contract')   ; call contract
case('wrong-form') ; call wrong_form
case default       ; error stop 'usage: cookbook_diagnostics steady|contract|wrong-form'
endselect

contains
   subroutine steady
   !< March u' = -(u - g) to its steady state g, stopping when the step residual max |dq| falls below a tolerance.
   integer(I4P), parameter        :: N=8
   real(R_P),    parameter        :: TOL=1.e-10_R_P
   real(R_P), allocatable, target :: u(:,:,:,:,:), dq(:,:,:,:,:), g(:,:,:,:,:)
   real(R_P), pointer             :: s_(:,:,:,:,:), r(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   integer(I4P)                   :: step, k, i

   allocate(u(1:1, 0:N+1, 0:2, 0:2, 1:1), dq(1:1, 0:N+1, 0:2, 0:2, 1:1), g(1:1, 0:N+1, 0:2, 0:2, 1:1))
   u = 0._R_P ; dq = 0._R_P
   g = 0._R_P ; g(1,1:N,1,1,1) = [(real(i, R_P) / N, i=1, N)]
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, N, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(SCHEME_SSP_33, FORM_BUTCHER, lb=lbound(u), ub=ubound(u))
   !region steady
   do step=1, 10000
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, 0.5_R_P, u)
         s_ => st%stage_state() ; r => st%rhs_register()
         r(1,1:N,1,1,1) = -(s_(1,1:N,1,1,1) - g(1,1:N,1,1,1))
         call st%end_stage(k, bx, 0.5_R_P, u)
      enddo
      call st%finish_step(bx, 0.5_R_P, u, dq=dq)
      if (maxval(abs(dq(1,1:N,1,1,1))) < TOL) exit             ! converged: the step changes nothing any more
   enddo
   !endregion steady
   print '(A,I0,A,ES9.2,A,ES9.2)', 'converged in ', step, ' steps: max |dq| = ', maxval(abs(dq(1,1:N,1,1,1))), &
         ', max |u - g| = ', maxval(abs(u(1,1:N,1,1,1) - g(1,1:N,1,1,1)))
   call st%destroy
   endsubroutine steady

   subroutine contract
   !< What a caller can check, cheaply, while it drives the stepper.
   real(R_P), allocatable, target :: q(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   integer(I4P)                   :: k

   allocate(q(1:1, 0:5, 0:2, 0:2, 1:1)) ; q = 0._R_P
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, 4_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q), stages_on_q=.true.)
   !region contract
   call st%begin_step(bx, q)
   do k=1, st%stages()
      call st%begin_stage(k, bx, 0.1_R_P, q)
      u => st%stage_state() ; r => st%rhs_register()
      ! the stage state is q itself or a register, never both; never the residual register
      if (st%stage_state_is_q(k) .eqv. associated(u)) error stop 'stage state ambiguous'
      if (associated(u, r)) error stop 'residual aliases the stage state'
      print '(A,I0,A,L1,A,L1)', 'stage ', k, ': state is q ', st%stage_state_is_q(k), ', stage_state() associated ', &
            associated(u)
      r = 0._R_P
      call st%end_stage(k, bx, 0.1_R_P, q)
   enddo
   call st%finish_step(bx, 0.1_R_P, q)
   !endregion contract
   call st%destroy
   endsubroutine contract

   subroutine wrong_form
   !< rk4 has a Butcher form only: asking for its two-register form stops the program with a message.
   type(stepper) :: st
   !region wrong-form
   call st%initialize(SCHEME_RK4, FORM_SHU_OSHER_2R, lb=[1, 0, 0, 0, 1], ub=[1, 5, 2, 2, 1])
   !endregion wrong-form
   endsubroutine wrong_form
endprogram cookbook_diagnostics
