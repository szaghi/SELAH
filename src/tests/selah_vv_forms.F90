!< SELAH verification C2: every storage form of a scheme computes the solution of its Butcher form, to round-off.
program selah_vv_forms
!< SELAH verification C2: every storage form of a scheme computes the solution of its Butcher form, to round-off.
!<
!< The two-register (euler, ssp-22, ssp-33), Williamson 2N (lsrk-33, lsrk-54) and Ketcheson (ssp-104) forms store other
!< coefficients and other registers, and add the stage terms in another order: they must give the same step. Every
!< scheme with a second form runs 40 steps of the periodic advection (64 cells, central residual, CFL 0.5, a profile
!< with a jump) in both forms through the stepper. Asserted: the largest difference of the two solutions over the run
!< stays below 1000 eps(R_P) times the largest |q|. Report mode (SELAH_VV_DIR): forms.dat (scheme, form, step,
!< difference).
use selah
use selah_test_tools, only : check, failures
use selah_vv_cases
implicit none
integer(I4P), parameter :: N=64, NSTEPS=40
character(7)            :: schemes(6) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_LSRK_33, &
                                                         SCHEME_LSRK_54, SCHEME_SSP_104]
character(14)           :: forms(6) = [character(14) :: FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, &
                                                        FORM_WILLIAMSON_2N, FORM_WILLIAMSON_2N, FORM_KETCHESON_2R]
real(R_P), allocatable  :: hb(:,:), hf(:,:)
real(R_P)               :: worst
integer                 :: s, step, u

if (vv_report()) call vv_open('forms.dat', u, 'scheme form step max|q_form - q_butcher|/max|q|')
do s=1, size(schemes)
   call run(trim(schemes(s)), FORM_BUTCHER, hb)
   call run(trim(schemes(s)), trim(forms(s)), hf)
   worst = 0._R_P
   do step=1, NSTEPS
      worst = max(worst, maxval(abs(hf(:,step) - hb(:,step))) / maxval(abs(hb(:,step))))
      if (vv_report()) write(u, '(A,1X,A,I5,ES14.5)') trim(schemes(s)), trim(forms(s)), step, &
                                 maxval(abs(hf(:,step) - hb(:,step))) / maxval(abs(hb(:,step)))
   enddo
   print '(A,ES10.2)', '        '//trim(schemes(s))//' '//trim(forms(s))//' vs butcher, worst relative difference:', worst
   call check(worst <= 1000._R_P * epsilon(1._R_P), trim(schemes(s))//' '//trim(forms(s))//': the Butcher solution to round-off')
enddo
if (vv_report()) close(u)
if (failures > 0) error stop 'selah_vv_forms: failed'
print '(A)', 'selah_vv_forms: ok'

contains
   subroutine run(name, form, hist)
   !< NSTEPS steps; hist(:,step) the solution in the cells after each step.
   character(*),           intent(in)  :: name, form
   real(R_P), allocatable, intent(out) :: hist(:,:)
   real(R_P), allocatable, target      :: q(:,:,:,:,:)
   real(R_P), pointer                  :: u(:,:,:,:,:), r(:,:,:,:,:)
   type(stepper)                       :: st
   type(box)                           :: bx
   real(R_P)                           :: dx, dt
   integer(I4P)                        :: step, k, i

   dx = 1._R_P / N ; dt = 0.5_R_P * dx
   allocate(q(1:1, 0:N+1, 0:2, 0:2, 1:1), hist(N, NSTEPS)) ; q = 0._R_P
   do i=1, N
      q(1,i,1,1,1) = sin(6.283185307179586_R_P * (i - 0.5_R_P) * dx) + merge(1._R_P, 0._R_P, i > N/2)
   enddo
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, N, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(name, form, lb=lbound(q), ub=ubound(q))
   do step=1, NSTEPS
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)
         if (st%stage_state_is_q(k)) then ; u => q ; else ; u => st%stage_state() ; endif
         r => st%rhs_register()
         u(1,0,1,1,1) = u(1,N,1,1,1) ; u(1,N+1,1,1,1) = u(1,1,1,1,1)
         do i=1, N
            r(1,i,1,1,1) = -(u(1,i+1,1,1,1) - u(1,i-1,1,1,1)) / (2._R_P * dx)
         enddo
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
      hist(:,step) = q(1,1:N,1,1,1)
   enddo
   call st%destroy
   endsubroutine run
endprogram selah_vv_forms
