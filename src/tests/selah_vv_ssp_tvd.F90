!< SELAH verification P1: strong stability. An SSP scheme keeps a forward-Euler TVD operator TVD up to C times its step.
program selah_vv_ssp_tvd
!< SELAH verification P1: strong stability. An SSP scheme keeps a forward-Euler TVD operator TVD up to C times its step.
!<
!< The spatial operator is first-order upwind for u_t + u_x = 0 on a periodic grid of N cells: a forward-Euler step is
!< TVD (it does not increase the total variation TV(u) = sum |u_{i+1} - u_i|) for CFL nu <= 1. By Gottlieb, Shu &
!< Tadmor (2001), a scheme with SSP coefficient C keeps it TVD for nu <= C. The initial data is a square wave, every
!< scheme runs through the stepper (Butcher form), NSTEPS steps per CFL, and the CFL is swept as a multiple of C.
!< Asserted: for every SSP scheme the total variation never grows (by more than 1000 eps of itself, round-off) for
!< nu <= C. Measured and reported, not asserted: the CFL above which it does grow (the scheme's threshold on this
!< linear problem, which may exceed C: C is a sufficient bound for every TVD operator), for the SSP schemes and the
!< others (rk4 and the 2N schemes, which have no such bound).
!< Report mode (SELAH_VV_DIR): ssp-tvd-sweep.dat (scheme, nu/C or nu, max relative TV growth), ssp-tvd-summary.dat,
!< and the frames of ssp-54 at nu = C and at 1.15 times its observed threshold (ssp-tvd-frames-*.dat).
use selah
use selah_test_tools, only : check, failures
use gideon, only : gideon_fv_face_ul1
use selah_vv_cases
implicit none
integer(I4P), parameter   :: N=100, NSTEPS=60, NSWEEP=40
real(R_P),    parameter   :: TOL=1000._R_P*epsilon(1._R_P) ! round-off of a 100-term sum
! a run whose TV passes BLOWUP x TV0 stops: 12 decades below overflow, above one step's growth (|R| <= 1.4e4 at nu = 2C,
! ssp-104) times the residual's 1/dx. R4P overflowed otherwise (ssp-22 at nu = 1.9: 5^59); R8P never reaches it
real(R_P),    parameter   :: BLOWUP=1.e-12_R_P*huge(1._R_P)
type(butcher_tableau)     :: bt
real(R_P)                 :: nu, growth, threshold, cref
logical                   :: tvd_ok
integer                   :: s, k, us, um
character(:), allocatable :: name

if (vv_report()) then
   call vv_open('ssp-tvd-sweep.dat', us, 'scheme nu/C (nu for non-SSP schemes) max_relative_TV_growth')
   call vv_open('ssp-tvd-summary.dat', um, 'scheme ssp_coefficient observed_threshold_nu')
endif
do s=1, size(RK_SCHEMES)
   name = trim(RK_SCHEMES(s))
   bt = butcher(name)
   cref = merge(real(bt%ssp_coefficient, R_P), 1._R_P, bt%ssp_coefficient > 0._R_P)  ! the sweep's unit
   tvd_ok = .true. ; threshold = -1._R_P
   do k=1, NSWEEP
      nu = cref * 0.05_R_P * k                                    ! 0.05 C .. 2 C
      growth = tv_growth(name, nu)
      if (bt%ssp_coefficient > 0._R_P .and. nu <= cref * (1._R_P + 1.e-12_R_P) .and. growth > TOL) tvd_ok = .false.
      if (threshold < 0._R_P .and. growth > TOL) threshold = nu
      if (vv_report()) write(us, '(A,2ES17.8E3)') name, nu / cref, growth
   enddo
   if (bt%ssp_coefficient > 0._R_P) then
      call check(tvd_ok, name//': TVD for every nu <= C (the SSP bound)')
      print '(A,F7.4,A,F7.3)', '        '//name//': C = ', cref, ', TV first grows at nu = ', threshold
   else
      print '(A,F7.3)', '        '//name//' (not SSP): TV first grows at nu = ', threshold
   endif
   if (vv_report()) write(um, '(A,2ES17.8E3)') name, real(bt%ssp_coefficient, R_P), threshold
   if (vv_report() .and. name == 'ssp-54') then
      growth = tv_growth(name, cref, 'ssp-tvd-frames-at-c.dat')
      growth = tv_growth(name, 1.15_R_P * threshold, 'ssp-tvd-frames-above.dat')
   endif
enddo
if (vv_report()) then
   close(us) ; close(um)
endif
if (failures > 0) error stop 'selah_vv_ssp_tvd: failed'
print '(A)', 'selah_vv_ssp_tvd: ok'

contains
   function tv_growth(name, nu, frames) result(growth)
   !< NSTEPS steps at CFL nu; the largest growth of TV over one step, relative to the initial TV.
   character(*), intent(in)           :: name
   real(R_P),    intent(in)           :: nu
   character(*), intent(in), optional :: frames
   real(R_P)                          :: growth
   real(R_P), allocatable, target     :: q(:,:,:,:,:)
   real(R_P), pointer                 :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                          :: fl, fr
   type(stepper)                      :: st
   type(box)                          :: bx
   real(R_P)                          :: tv0, tv_old, tv_new, dx, dt
   integer(I4P)                       :: step, k, i, uf

   dx = 1._R_P / N ; dt = nu * dx
   allocate(q(1:1, 0:N+1, 0:2, 0:2, 1:1)) ; q = 0._R_P
   do i=1, N
      q(1,i,1,1,1) = merge(1._R_P, 0._R_P, i > N/4 .and. i <= N/2)
   enddo
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, N, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(name, FORM_BUTCHER, lb=lbound(q), ub=ubound(q), stages_on_q=.true.)
   tv0 = tv(q) ; tv_old = tv0 ; growth = 0._R_P
   if (present(frames)) then
      call vv_open(frames, uf, 'line 1: x; then t, u(1..N), exact(1..N); nu = '//trim(adjustl(str(nu))))
      write(uf, '(*(ES14.6))') [((i - 0.5_R_P) * dx, i=1, N)]
      call write_frame(uf, 0._R_P, q, dx)
   endif
   do step=1, NSTEPS
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)
         if (st%stage_state_is_q(k)) then ; u => q ; else ; u => st%stage_state() ; endif
         r => st%rhs_register()
         u(1,0,1,1,1) = u(1,N,1,1,1)                               ! periodic
         do i=1, N                                                 ! first-order upwind, GIDEON's face values
            call gideon_fv_face_ul1([u(1,i-1,1,1,1)], fl) ; call gideon_fv_face_ul1([u(1,i,1,1,1)], fr)
            r(1,i,1,1,1) = -(fr - fl) / dx
         enddo
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
      tv_new = tv(q)
      growth = max(growth, (tv_new - tv_old) / tv0)
      tv_old = tv_new
      if (present(frames)) call write_frame(uf, step * dt, q, dx)
      if (tv_new > BLOWUP * tv0) exit
   enddo
   if (present(frames)) close(uf)
   call st%destroy
   endfunction tv_growth

   subroutine write_frame(uf, t, q, dx)
   !< One frame: t, the solution in the cells, the exact (translated) square wave.
   integer,   intent(in) :: uf
   real(R_P), intent(in) :: t, dx
   real(R_P), intent(in) :: q(1:,0:,0:,0:,1:)
   real(R_P)             :: x, ex(N)
   integer(I4P)          :: i
   do i=1, N
      x = modulo((i - 0.5_R_P) * dx - t, 1._R_P)
      ex(i) = merge(1._R_P, 0._R_P, x > 0.25_R_P .and. x <= 0.5_R_P)
   enddo
   write(uf, '(*(ES14.6))') t, q(1,1:N,1,1,1), ex
   endsubroutine write_frame

   pure function tv(q) result(v)
   real(R_P), intent(in) :: q(1:,0:,0:,0:,1:)
   real(R_P)             :: v
   integer(I4P)          :: i
   v = abs(q(1,1,1,1,1) - q(1,N,1,1,1))
   do i=1, N-1
      v = v + abs(q(1,i+1,1,1,1) - q(1,i,1,1,1))
   enddo
   endfunction tv

   function str(x) result(s)
   real(R_P), intent(in) :: x
   character(16)         :: s
   write(s, '(F8.4)') x
   endfunction str
endprogram selah_vv_ssp_tvd
