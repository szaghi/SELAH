!< SELAH test: the per-stage forward-Euler step and the SSP coefficient, and a positivity limiter that uses them.

#include "fundal.H"
program selah_test_stage_euler
!< SELAH test: the per-stage forward-Euler step and the SSP coefficient, and a positivity limiter that uses them.
!<
!< A limiter inside an SSP Runge–Kutta step (Zhang & Shu 2010) makes the state of each stage plus a forward-Euler step
!< of its residual admissible; the step it must use is the one the residual enters in the Shu–Osher form,
!< f_k dt = max_i beta_ik / alpha_ik dt (`stage_euler_fraction(k)`), not dt (ADAM's FLUME took dt for every stage).
!<
!<+ **Values.** f_k for euler, ssp-22, ssp-33 (all 1), ssp-104 (all 1/6, Ketcheson 2008), ssp-54 (the ratios of the
!<  Spiteri–Ruuth coefficients, row by row); the SSP coefficient 1, 1, 1, 6 and 1.508 (Spiteri & Ruuth 2002, three
!<  digits); max_k f_k = 1/C; the same values in every storage form of a scheme; 0 for rk4 and the 2N schemes (not
!<  SSP); a splitting's flows c_k (Strang: 1/2, 1, 1/2); planted Shu–Osher forms: the largest of a stage's steps, and
!<  0 for a form that is not a convex combination of forward-Euler steps.
!<+ **A positivity limiter through the stepper** (host or device). The decay u_i' = -lambda_i u_i, lambda_i dt = 0.5,
!<  2, 5, 20 (forward Euler is positive for lambda tau <= 1 only), the residual of stage k scaled by
!<  theta = min(1, 1/(lambda tau_k)) so that u + tau_k R(u) >= 0. With tau_k = f_k dt every stage state and every step
!<  stays non-negative, in every SSP scheme and form; with tau_k = f_k dt / 2 some stage state goes negative (the step is
!<  necessary, not only sufficient); with tau_k = dt (ADAM's choice) the states stay positive, limited more than needed
!<  for ssp-54 and ssp-104 (f_k < 1).
use penf, only : I4P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
character(7),  parameter :: SSP_SCHEMES(5) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, &
                                              SCHEME_SSP_54, SCHEME_SSP_104]
integer(I4P),  parameter :: N = 4                                       !< Cells.
real(R_P),     parameter :: LAMBDA_DT(N) = [0.5_R_P, 2._R_P, 5._R_P, 20._R_P] !< lambda_i dt.
real(R_P),     parameter :: TOL = 100._R_P * epsilon(1._R_P)            !< Round-off.
integer(I4P)             :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'forward-Euler fractions and SSP coefficients:'
call values
print '(A)', 'positivity limiter (tau_k = f_k dt; f_k dt / 2; dt):'
do s=1, size(SSP_SCHEMES)
   call limited(trim(SSP_SCHEMES(s)), FORM_BUTCHER)
enddo
call limited(SCHEME_EULER, FORM_SHU_OSHER_2R)
call limited(SCHEME_SSP_22, FORM_SHU_OSHER_2R)
call limited(SCHEME_SSP_33, FORM_SHU_OSHER_2R)
call limited(SCHEME_SSP_104, FORM_KETCHESON_2R)
if (failures > 0) error stop 'selah_test_stage_euler: failed'
print '(A)', 'selah_test_stage_euler: ok'

contains
   subroutine values
   !< The fractions and the SSP coefficients against the literature.
   real(R_P), allocatable :: f(:)
   real(R_P)              :: a(6,5), b(6,5)
   type(stepper)          :: st
   type(shu_osher_form)   :: so
   integer(I4P)           :: k

   call expect(SCHEME_EULER, FORM_BUTCHER, [1._R_P], 1._R_P)
   call expect(SCHEME_SSP_22, FORM_BUTCHER, [1._R_P, 1._R_P], 1._R_P)
   call expect(SCHEME_SSP_33, FORM_BUTCHER, [1._R_P, 1._R_P, 1._R_P], 1._R_P)
   call expect(SCHEME_SSP_104, FORM_BUTCHER, [(1._R_P/6._R_P, k=1, 10)], 6._R_P)
   ! Spiteri & Ruuth (2002) SSPRK(5,4), Shu–Osher coefficients (alpha, beta), row by row
   a = 0._R_P ; b = 0._R_P
   a(2,1) = 1._R_P                ; b(2,1) = 0.391752226571890_R_P
   a(3,2) = 0.555629506348765_R_P ; b(3,2) = 0.368410593050371_R_P
   a(4,3) = 0.379898148511597_R_P ; b(4,3) = 0.251891774271694_R_P
   a(5,4) = 0.821920045606868_R_P ; b(5,4) = 0.544974750228521_R_P
   a(6,4) = 0.096059710526147_R_P ; b(6,4) = 0.063692468666290_R_P
   a(6,5) = 0.386708617503268_R_P ; b(6,5) = 0.226007483236906_R_P
   f = [b(2,1)/a(2,1), b(3,2)/a(3,2), b(4,3)/a(4,3), max(b(5,4)/a(5,4), b(6,4)/a(6,4)), b(6,5)/a(6,5)]
   call expect(SCHEME_SSP_54, FORM_BUTCHER, f, 1.508_R_P, ctol=1.e-3_R_P)
   ! every storage form of a scheme: the same stages, the same fractions
   call expect(SCHEME_EULER, FORM_SHU_OSHER_2R, [1._R_P], 1._R_P)
   call expect(SCHEME_SSP_22, FORM_SHU_OSHER_2R, [1._R_P, 1._R_P], 1._R_P)
   call expect(SCHEME_SSP_33, FORM_SHU_OSHER_2R, [1._R_P, 1._R_P, 1._R_P], 1._R_P)
   call expect(SCHEME_SSP_104, FORM_KETCHESON_2R, [(1._R_P/6._R_P, k=1, 10)], 6._R_P)
   ! not SSP: no forward-Euler decomposition
   call expect(SCHEME_RK4, FORM_BUTCHER, [(0._R_P, k=1, 4)], 0._R_P)
   call expect(SCHEME_LSRK_33, FORM_WILLIAMSON_2N, [(0._R_P, k=1, 3)], 0._R_P)
   call expect(SCHEME_LSRK_54, FORM_BUTCHER, [(0._R_P, k=1, 5)], 0._R_P)
   ! planted forms: stage 1 enters two steps (ratios 1 and 1/2): the larger one; a negative coefficient, and a residual
   ! without its state (beta > 0, alpha = 0): no decomposition
   so%stages = 2
   allocate(so%alpha(3,2), so%beta(3,2))
   so%alpha = 0._R_P ; so%beta = 0._R_P
   so%alpha(2,1) = 1._R_P   ; so%beta(2,1) = 1._R_P
   so%alpha(3,1) = 0.5_R_P  ; so%beta(3,1) = 0.25_R_P
   so%alpha(3,2) = 0.5_R_P  ; so%beta(3,2) = 0.5_R_P
   call check(all(euler_step_fractions(so) == [1._R_P, 1._R_P]), 'planted form: the largest step of a stage''s residual')
   so%beta(3,1) = -0.25_R_P
   call check(all(euler_step_fractions(so) == 0._R_P), 'planted form, a negative coefficient: no decomposition')
   so%beta(3,1) = 0.25_R_P ; so%alpha(3,1) = 0._R_P ; so%alpha(3,2) = 1._R_P
   call check(all(euler_step_fractions(so) == 0._R_P), 'planted form, beta > 0 with alpha = 0: no decomposition')
   ! a splitting: flow k is an Euler step of fraction c_k on its partition
   call st%initialize(SCHEME_STRANG, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], components=.true., &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   call check(all([(st%stage_euler_fraction(k), k=1, 3)] == [0.5_R_P, 1._R_P, 0.5_R_P]) .and. &
              st%ssp_coefficient() == 0._R_P, 'strang: flows 1/2, 1, 1/2; no SSP coefficient')
   call st%destroy
   endsubroutine values

   subroutine expect(name, form, f, c, ctol)
   !< The stepper's fractions equal f (to round-off), its SSP coefficient c (to ctol), and max f = 1/C when C > 0.
   character(*), intent(in)           :: name, form
   real(R_P),    intent(in)           :: f(:)
   real(R_P),    intent(in)           :: c
   real(R_P),    intent(in), optional :: ctol
   type(stepper)                      :: st
   real(R_P), allocatable             :: g(:)
   real(R_P)                          :: tc
   integer(I4P)                       :: k

   tc = TOL ; if (present(ctol)) tc = ctol
   call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
   g = [(st%stage_euler_fraction(k), k=1, st%stages())]
   print '(A,*(F9.6))', '        '//name//' '//form//': f_k =', g
   print '(A,F9.6)', '        '//repeat(' ', len(name//form))//'  C =', st%ssp_coefficient()
   call check(size(g) == size(f), name//' '//form//': one fraction per stage')
   if (size(g) == size(f)) call check(all(abs(g - f) <= TOL), name//' '//form//': forward-Euler fractions')
   call check(abs(st%ssp_coefficient() - c) <= tc, name//' '//form//': SSP coefficient')
   if (st%ssp_coefficient() > 0._R_P) &
      call check(abs(maxval(g) * st%ssp_coefficient() - 1._R_P) <= TOL, name//' '//form//': max f_k = 1/C')
   call st%destroy
   endsubroutine expect

   subroutine limited(name, form)
   !< The limited decay with tau_k = f_k dt (positive), f_k dt / 2 (some stage negative), dt (positive).
   character(*), intent(in) :: name, form
   real(R_P)                :: lo_exact, lo_half, lo_dt

   lo_exact = limited_min(name, form, 1._R_P, .false.)
   lo_half = limited_min(name, form, 0.5_R_P, .false.)
   lo_dt = limited_min(name, form, 1._R_P, .true.)
   print '(A,3ES11.2)', '        '//name//' '//form//': min state:', lo_exact, lo_half, lo_dt
   call check(lo_exact >= -TOL, name//' '//form//': positive with tau_k = f_k dt')
   call check(lo_half < -0.1_R_P, name//' '//form//': negative with tau_k = f_k dt / 2')
   call check(lo_dt >= -TOL, name//' '//form//': positive with tau_k = dt')
   endsubroutine limited

   function limited_min(name, form, scale, whole_dt) result(lo)
   !< Four steps of the limited decay; the least value of every stage state and step, relative to u(0) = 1.
   character(*), intent(in) :: name, form
   real(R_P),    intent(in) :: scale    !< tau_k = scale f_k dt (or scale dt).
   logical,      intent(in) :: whole_dt !< tau_k = scale dt for every stage.
   real(R_P)                :: lo
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q(:)=>null(), u(:), r(:), lam(:)=>null()
   real(R_P), allocatable   :: h(:), hl(:)
   real(R_P)                :: dt, tau
   integer(I4P)             :: step, k, ierr

   dt = 1._R_P
   allocate(h(N), hl(N))
   h = 1._R_P
   hl = LAMBDA_DT / dt
   call dev_alloc(fptr_dev=q, ubounds=[N], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=lam, ubounds=[N], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=q, src=h)
   call dev_memcpy_to_device(dst=lam, src=hl)
   bx = box_interior(lb=[1_I4P], ub=[N], ngc=0_I4P, components=.false.)
   call st%initialize(name, form, lb=[1_I4P], ub=[N], components=.false.)
   lo = 1._R_P
   do step=1, 4
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)
         if (st%stage_state_is_q(k)) then
            u => q
         else
            call st%stage_state_view(u)
         endif
         call dev_memcpy_from_device(dst=h, src=u)
         lo = min(lo, minval(h))
         call st%rhs_register_view(r)
         tau = scale * dt * merge(1._R_P, st%stage_euler_fraction(k), whole_dt)
         call limited_decay(tau, lam, u, r)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
      call dev_memcpy_from_device(dst=h, src=q)
      lo = min(lo, minval(h))
   enddo
   call st%destroy
   call dev_free(q) ; call dev_free(lam)
   endfunction limited_min

   subroutine limited_decay(tau, lam, u, r)
   !< r = theta (-lambda u), theta = min(1, 1/(lambda tau)): u + tau r >= 0 for u >= 0.
   real(R_P), intent(in)    :: tau
   real(R_P), intent(in)    :: lam(1:), u(1:)
   real(R_P), intent(inout) :: r(1:)
   integer(I4P)             :: i

   !$acc parallel loop independent DEVICEVAR(lam, u, r)
   !$omp OMPLOOP DEVICEPTR(lam, u, r)
   do i=1, N
      r(i) = -min(1._R_P, 1._R_P / (lam(i) * tau)) * lam(i) * u(i)
   enddo
   endsubroutine limited_decay
endprogram selah_test_stage_euler
