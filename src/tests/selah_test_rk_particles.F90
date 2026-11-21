!< SELAH test: particles pushed by the field's Runge–Kutta scheme on a companion, and the axis stability bounds.

#include "fundal.H"
program selah_test_rk_particles
!< SELAH test: particles pushed by the field's Runge–Kutta scheme on a companion, and the axis stability bounds.
!<
!< PRISM pushes its PIC particles with the field's scheme, in step with it (now: a companion, selah_test_companion).
!< The properties such a push has, from the stability polynomial R of the scheme:
!<
!<+ **Axis bounds.** `imag_axis_bound()` against the closed forms (euler and ssp-22 0, ssp-33 and lsrk-33 sqrt(3), rk4
!<  2 sqrt(2)) and the computed values (ssp-54 3.2783556, PRISM's 3.278; ssp-104 4.9214531; lsrk-54 3.3407180);
!<  `real_axis_bound()` (2, 2.512745, 2.785294, 5.331473, 13.917047, 4.656757); rk4 with a rounded weight keeps its
!<  interval.
!<+ **Gyration.** A charge in a uniform B (w = v_x + i v_y, w' = -i omega w), pushed on a companion of a field
!<  stepper: after n steps w_n = R(-i omega dt)^n w_0 to round-off (R from the tableau by a complex triangular solve,
!<  not the library's polynomial): the energy decays as |R|^(2n), the phase advances by arg R per step.
!<+ **Stability at the bound.** At omega dt = beta (1 - 1e-6) the speed never grows over 10^4 steps; at
!<  beta (1 + 1e-3) it grows; euler and ssp-22 grow at any omega dt (0.1 here). R8P builds only (R4P cannot resolve
!<  the bound to 1e-6).
!<+ **ExB drift.** In uniform E perpendicular to B the drift velocity is an equilibrium of every explicit scheme: a
!<  particle started at v_E keeps it to round-off (8 eps; bit for bit where the form adds dt times a zero residual, not
!<  in the Ketcheson form, whose (u^n + 9 q1)/25 rounds) and its position advances by dt v_E per step, every scheme and
!<  form.
use penf, only : I4P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
character(7),  parameter :: OTHER(3) = [character(7) :: SCHEME_SSP_33, SCHEME_LSRK_54, SCHEME_SSP_104]
character(14), parameter :: OTHER_FORMS(3) = [character(14) :: FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_KETCHESON_2R]
integer(I4P)             :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'axis bounds:'
call bounds
print '(A)', 'gyration on a companion (w_n = R(-i omega dt)^n w_0):'
do s=1, size(RK_SCHEMES)
   call gyration(trim(RK_SCHEMES(s)), FORM_BUTCHER)
enddo
do s=1, size(OTHER)
   call gyration(trim(OTHER(s)), trim(OTHER_FORMS(s)))
enddo
print '(A)', 'ExB drift:'
do s=1, size(RK_SCHEMES)
   call drift(trim(RK_SCHEMES(s)), FORM_BUTCHER)
enddo
do s=1, size(OTHER)
   call drift(trim(OTHER(s)), trim(OTHER_FORMS(s)))
enddo
if (failures > 0) error stop 'selah_test_rk_particles: failed'
print '(A)', 'selah_test_rk_particles: ok'

contains
   subroutine bounds
   !< The stepper's axis bounds against the closed forms and the computed values.
   character(7), parameter :: NAMES(8) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_LSRK_33, &
                                          SCHEME_RK4, SCHEME_SSP_54, SCHEME_SSP_104, SCHEME_LSRK_54]
   real(R_P) :: bi(8), br(8), tol, tol_r
   type(stepper) :: st
   type(butcher_tableau) :: bt
   integer(I4P)  :: k

   bi = [0._R_P, 0._R_P, sqrt(3._R_P), sqrt(3._R_P), 2._R_P*sqrt(2._R_P), 3.2783556_R_P, 4.9214531_R_P, 3.3407180_R_P]
   br = [2._R_P, 2._R_P, 2.512745_R_P, 2.512745_R_P, 2.785294_R_P, 5.331473_R_P, 13.917047_R_P, 4.656757_R_P]
   tol = max(1.e-6_R_P, 1.e3_R_P * epsilon(1._R_P))
   do k=1, size(NAMES)
      call st%initialize(trim(NAMES(k)), FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
      print '(A,2F12.7)', '        '//trim(NAMES(k))//': imaginary and real axis bounds', st%imag_axis_bound(), &
                          st%real_axis_bound()
      call check(abs(st%imag_axis_bound() - bi(k)) <= tol, trim(NAMES(k))//': imag_axis_bound')
      ! the real bound is a root of R(-x) = +-1: rounding resolves it only to eps sum|g_i| beta^i / |R'(-beta)|, 1.5e-3
      ! for ssp-104 in R4P (degree 10, terms of 4e4 cancelling to 1), above 1e3 eps
      tol_r = max(tol, epsilon(1._R_P) * root_condition(stability_polynomial(butcher(trim(NAMES(k)))), br(k)))
      call check(abs(st%real_axis_bound() - br(k)) <= tol_r, trim(NAMES(k))//': real_axis_bound')
      call st%destroy
   enddo
   ! a tableau whose digits round: rk4 with b_1 off by 4 eps. |R(iy)|^2 - 1 then has a tiny positive y^2 term, which
   ! must not close rk4's interval (the order's exact zeros are imposed)
   bt = butcher(SCHEME_RK4)
   bt%b(1) = bt%b(1) + 4._R_P * epsilon(1._R_P)
   call check(abs(imag_axis_bound(bt) - 2._R_P*sqrt(2._R_P)) <= tol, 'rk4, b_1 rounded up: the interval kept')
   call st%initialize(SCHEME_STRANG, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], components=.true., &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   call check(st%imag_axis_bound() == 0._R_P .and. st%real_axis_bound() == 0._R_P, 'strang: no axis bounds (splitting)')
   call st%destroy
   endsubroutine bounds

   pure function root_condition(g, beta) result(c)
   !< The condition of the root beta of |R(-x)| = 1: sum_i |g_i| beta^i / |R'(-beta)|.
   real(R_P), intent(in) :: g(0:), beta
   real(R_P)             :: c, s, d
   integer(I4P)          :: i

   s = abs(g(0)) ; d = 0._R_P
   do i=1, ubound(g, 1)
      s = s + abs(g(i)) * beta**i
      d = d + i * g(i) * (-beta)**(i - 1)
   enddo
   c = s / abs(d)
   endfunction root_condition

   function r_of(name, z) result(r)
   !< R(z) = 1 + z b^T (I - z A)^{-1} 1, by forward substitution (A strictly lower triangular).
   character(*), intent(in) :: name
   complex(R_P), intent(in) :: z
   complex(R_P)             :: r
   type(butcher_tableau)    :: bt
   complex(R_P), allocatable :: k(:)
   integer(I4P)             :: m

   bt = butcher(name)
   allocate(k(bt%stages))
   do m=1, bt%stages ! k_m = 1 + z sum_{j<m} a_mj k_j (m: a stage)
      k(m) = 1._R_P + z * sum(bt%a(m,1:m-1) * k(1:m-1))
   enddo
   r = 1._R_P + z * sum(bt%b * k)
   endfunction r_of

   subroutine push(st, sp, dt, y, p, w2, e0, n)
   !< n steps of a trivial field (residual 0) and its companion, particles p = (x, y, vx, vy): x' = v,
   !< v' = (E + v x B) with B = (0, 0, w2) (q/m = 1) and E = (0, e0).
   type(stepper),      intent(inout) :: st, sp
   real(R_P),          intent(in)    :: dt, w2, e0
   real(R_P), pointer, intent(inout) :: y(:), p(:,:)
   integer(I4P),       intent(in)    :: n
   type(box)                         :: by, bp
   real(R_P), pointer                :: q(:,:), ry(:), rp(:,:)
   integer(I4P)                      :: step, k

   by = box_interior(lb=[1_I4P], ub=[1_I4P], ngc=0_I4P, components=.false.)
   bp = box_interior(lb=[1_I4P, 1_I4P], ub=[4_I4P, 1_I4P], ngc=0_I4P, components=.true.)
   do step=1, n
      call st%begin_step(by, y) ; call sp%begin_step(bp, p)
      do k=1, st%stages()
         call st%begin_stage(k, by, dt, y) ; call sp%begin_stage(k, bp, dt, p)
         if (sp%stage_state_is_q(k)) then
            q => p
         else
            call sp%stage_state_view(q)
         endif
         call st%rhs_register_view(ry) ; call sp%rhs_register_view(rp)
         call lorentz(w2, e0, q, rp, ry)
         call st%end_stage(k, by, dt, y) ; call sp%end_stage(k, bp, dt, p)
      enddo
      call st%finish_step(by, dt, y) ; call sp%finish_step(bp, dt, p)
   enddo
   endsubroutine push

   subroutine lorentz(w2, e0, q, rp, ry)
   !< The particle's residual (x' = v, v' = E + v x B) and the trivial field's (0).
   real(R_P), intent(in)    :: w2, e0
   real(R_P), intent(in)    :: q(1:,1:)
   real(R_P), intent(inout) :: rp(1:,1:), ry(1:)

   !$acc serial DEVICEVAR(q, rp, ry)
   rp(1,1) = q(3,1)
   rp(2,1) = q(4,1)
   rp(3,1) = q(4,1) * w2
   rp(4,1) = e0 - q(3,1) * w2
   ry(1) = 0._R_P
   !$acc end serial
   endsubroutine lorentz

   subroutine gyration(name, form)
   !< w_n = R(-i y)^n w_0, y = omega dt = 0.5 (when stable) or 0.1; the stability at the bound.
   character(*), intent(in) :: name, form
   type(stepper)            :: st, sp
   real(R_P), pointer       :: y(:)=>null(), p(:,:)=>null()
   real(R_P)                :: h(4,1), yy, beta, s0, s1, smax
   complex(R_P)             :: w0, wn, we
   integer(I4P)             :: ierr, n, i

   call dev_alloc(fptr_dev=y, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=p, ubounds=[4_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
   call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[4_I4P, 1_I4P], components=.true.)
   beta = st%imag_axis_bound()
   yy = merge(0.5_R_P, 0.1_R_P, beta > 0.5_R_P)
   n = 50
   h(:,1) = [0._R_P, 0._R_P, 1._R_P, 0.5_R_P]
   call dev_memcpy_to_device(dst=p, src=h) ; call dev_memcpy_to_device(dst=y, src=[0._R_P])
   call push(st, sp, 1._R_P, y, p, yy, 0._R_P, n)
   call dev_memcpy_from_device(dst=h, src=p)
   w0 = cmplx(1._R_P, 0.5_R_P, R_P)
   wn = cmplx(h(3,1), h(4,1), R_P)
   we = r_of(name, cmplx(0._R_P, -yy, R_P))**n * w0
   print '(A,F4.1,A,2ES11.3)', '        '//name//' '//form//': omega dt =', yy, ': |w_n|/|w_0| and |w_n - R^n w_0|', &
                               abs(wn)/abs(w0), abs(wn - we)
   call check(abs(wn - we) <= 1.e3_R_P * epsilon(1._R_P) * abs(w0), name//' '//form//': w_n = R(-i omega dt)^n w_0')
#if defined _R_P_IS_R4P
   ! stability at the bound: below R4P's resolution of the bound (1e-6) and of the growth per step
   if (.false.) then
#endif
   ! stability at the bound: 10^4 steps
   do i=1, 2
      if (beta > 0._R_P) then
         yy = beta * merge(1._R_P - 1.e-6_R_P, 1._R_P + 1.e-3_R_P, i == 1)
      else
         yy = 0.1_R_P
      endif
      h(:,1) = [0._R_P, 0._R_P, 1._R_P, 0._R_P]
      call dev_memcpy_to_device(dst=p, src=h)
      s0 = 1._R_P ; smax = 0._R_P
      do n=1, 10
         call push(st, sp, 1._R_P, y, p, yy, 0._R_P, 1000_I4P)
         call dev_memcpy_from_device(dst=h, src=p)
         s1 = sqrt(h(3,1)**2 + h(4,1)**2)
         smax = max(smax, s1)
      enddo
      if (beta > 0._R_P .and. i == 1) then
         call check(smax <= s0 * (1._R_P + 1.e3_R_P * epsilon(1._R_P)), name//' '//form// &
                    ': no growth at omega dt = beta (1 - 1e-6) over 10^4 steps')
      elseif (beta > 0._R_P) then
         call check(s1 > s0, name//' '//form//': growth at omega dt = beta (1 + 1e-3)')
      elseif (i == 1) then
         call check(s1 > s0, name//' '//form//': no imaginary-axis interval, growth at omega dt = 0.1')
      endif
   enddo
#if defined _R_P_IS_R4P
   endif
#endif
   call st%destroy ; call sp%destroy
   call dev_free(y) ; call dev_free(p)
   endsubroutine gyration

   subroutine drift(name, form)
   !< Started at v_E = (e0/w2, 0), the particle keeps v_E bit for bit and advances by dt v_E per step.
   character(*), intent(in) :: name, form
   real(R_P),    parameter  :: W2 = 2._R_P, E0 = 0.5_R_P, DT = 0.25_R_P
   type(stepper)            :: st, sp
   real(R_P), pointer       :: y(:)=>null(), p(:,:)=>null()
   real(R_P)                :: h(4,1)
   integer(I4P)             :: ierr

   call dev_alloc(fptr_dev=y, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=p, ubounds=[4_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
   call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[4_I4P, 1_I4P], components=.true.)
   h(:,1) = [0._R_P, 0._R_P, E0/W2, 0._R_P]
   call dev_memcpy_to_device(dst=p, src=h) ; call dev_memcpy_to_device(dst=y, src=[0._R_P])
   call push(st, sp, DT, y, p, W2, E0, 64_I4P)
   call dev_memcpy_from_device(dst=h, src=p)
   print '(A,ES10.2)', '        '//name//' '//form//': |v - v_E|/|v_E| after 64 steps', &
                       sqrt((h(3,1) - E0/W2)**2 + h(4,1)**2) / (E0/W2)
   call check(sqrt((h(3,1) - E0/W2)**2 + h(4,1)**2) <= 64*epsilon(1._R_P)*E0/W2, &
              name//' '//form//': v_E kept to round-off (64 eps over 64 steps)')
   call check(abs(h(1,1) - 64*DT*E0/W2) <= 64*epsilon(1._R_P)*64*DT*E0/W2 .and. abs(h(2,1)) <= 64*epsilon(1._R_P), &
              name//' '//form//': the position advances by dt v_E per step')
   call st%destroy ; call sp%destroy
   call dev_free(y) ; call dev_free(p)
   endsubroutine drift
endprogram selah_test_rk_particles
