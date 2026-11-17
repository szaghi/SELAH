!< SELAH verification P3: over long times a symplectic splitting keeps the energy bounded; a Runge–Kutta scheme drifts.
program selah_vv_splitting_energy
!< SELAH verification P3: over long times a symplectic splitting keeps the energy bounded; a Runge–Kutta scheme drifts.
!<
!< The harmonic oscillator q' = p, p' = -q, H = (q^2 + p^2)/2, from (q, p) = (1, 0) with dt = 0.1 for 1000 periods
!< (62832 steps), through the stepper: the splittings (strang, yoshida-4, blanes-moan-s6) with their partitions (q, p),
!< and two Runge–Kutta schemes of order 4 (rk4, ssp-54) in Butcher form. The energy error of a symplectic method
!< oscillates with a bounded amplitude; a non-symplectic one drifts (for a linear problem |R(i dt)| < 1: the energy
!< decays every step). Asserted, comparing the largest |H - H0| over the last tenth of the run with that over the first
!< tenth: at most 2 for every splitting, at least 5 for the two Runge–Kutta schemes.
!< Report mode (SELAH_VV_DIR): splitting-energy-<scheme>.dat (t, H - H0, every 200 steps) and the phase portrait of
!< strang and rk4 at a large step, dt = 0.6 (splitting-phase-<scheme>.dat).
use selah
use selah_test_tools, only : check, failures
use selah_vv_cases
implicit none
real(R_P),    parameter :: PI=acos(-1._R_P), DT=0.1_R_P
integer(I4P), parameter :: NSTEPS=nint(2000*PI/DT)
character(15)           :: names(5) = [character(15) :: SCHEME_STRANG, SCHEME_YOSHIDA_4, SCHEME_BLANES_MOAN_S6, &
                                       SCHEME_RK4, SCHEME_SSP_54]
real(R_P)               :: first, last
integer                 :: s

do s=1, size(names)
   call run(trim(names(s)), DT, NSTEPS, first, last, 'splitting-energy-'//trim(names(s))//'.dat')
   print '(A,2ES10.2,A,F9.2)', '        '//trim(names(s))//' max |H - H0| first and last tenth:', first, last, &
         ', ratio', last / first
   if (s <= 3) then
      call check(last <= 2._R_P * first, trim(names(s))//': energy error bounded (last tenth <= 2 x first tenth)')
   else
      call check(last >= 5._R_P * first, trim(names(s))//': energy drifts (last tenth >= 5 x first tenth)')
   endif
enddo
if (vv_report()) then
   call run(SCHEME_STRANG, 0.6_R_P, 400_I4P, first, last, 'splitting-phase-strang.dat')
   call run(SCHEME_RK4, 0.6_R_P, 400_I4P, first, last, 'splitting-phase-rk4.dat')
endif
if (failures > 0) error stop 'selah_vv_splitting_energy: failed'
print '(A)', 'selah_vv_splitting_energy: ok'

contains
   subroutine run(name, dt, nsteps, first, last, report)
   !< nsteps steps of the oscillator; the largest energy error over the first and the last tenth of them.
   character(*), intent(in)       :: name, report
   real(R_P),    intent(in)       :: dt
   integer(I4P), intent(in)       :: nsteps
   real(R_P),    intent(out)      :: first, last
   real(R_P), allocatable, target :: y(:,:,:,:,:)
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P)                      :: h0, e
   integer(I4P)                   :: step, k
   integer                        :: ur
   logical                        :: split

   split = any(name == SPLITTING_SCHEMES)
   allocate(y(1:2, 0:2, 0:2, 0:2, 1:1)) ; y = 0._R_P ; y(1,1,1,1,1) = 1._R_P
   bx = box_interior(LAYOUT_V_FIRST, 2_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   if (split) then
      call st%initialize(name, FORM_SPLITTING, lb=lbound(y), ub=ubound(y), partitions=reshape([1, 1, 2, 2], [2, 2]))
   else
      call st%initialize(name, FORM_BUTCHER, lb=lbound(y), ub=ubound(y), stages_on_q=.true.)
   endif
   h0 = energy(y) ; first = 0._R_P ; last = 0._R_P
   if (vv_report()) call vv_open(report, ur, 't H-H0 q p, scheme '//name)
   do step=1, nsteps
      call st%begin_step(bx, y)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, y)
         if (st%stage_state_is_q(k)) then ; u => y ; else ; u => st%stage_state() ; endif
         r => st%rhs_register()
         select case(st%stage_partition(k))
         case(0) ; r(1,1,1,1,1) = u(2,1,1,1,1) ; r(2,1,1,1,1) = -u(1,1,1,1,1)   ! Runge-Kutta: both components
         case(1) ; r(1,1,1,1,1) = u(2,1,1,1,1)                                  ! the flow of q
         case(2) ; r(2,1,1,1,1) = -u(1,1,1,1,1)                                 ! the flow of p
         endselect
         call st%end_stage(k, bx, dt, y)
      enddo
      call st%finish_step(bx, dt, y)
      e = abs(energy(y) - h0)
      if (step <= nsteps / 10) first = max(first, e)
      if (step > nsteps - nsteps / 10) last = max(last, e)
      if (vv_report() .and. (mod(step, max(nsteps / 400, 1)) == 0)) &
         write(ur, '(4ES16.8)') step * dt, energy(y) - h0, y(1,1,1,1,1), y(2,1,1,1,1)
   enddo
   if (vv_report()) close(ur)
   call st%destroy
   endsubroutine run

   pure function energy(y) result(h)
   real(R_P), intent(in) :: y(1:,0:,0:,0:,1:)
   real(R_P)             :: h
   h = 0.5_R_P * (y(1,1,1,1,1)**2 + y(2,1,1,1,1)**2)
   endfunction energy
endprogram selah_vv_splitting_energy
