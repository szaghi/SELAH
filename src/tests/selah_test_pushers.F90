!< SELAH test: the two-level particle pushers (Boris, relativistic Boris, Vay, Higuera–Cary), on the host or the device.

#include "fundal.H"
program selah_test_pushers
!< SELAH test: the two-level particle pushers (Boris, relativistic Boris, Vay, Higuera–Cary), on the host or the device.
!<
!< One particle (or a few) per call, q/m = 1, c = 1 (unless stated); x at components 1:3, v or u at 4:6, q and m at 7, 8.
!<
!<+ **Boris.** In a pure B the speed is kept to round-off over 1000 steps and the velocity turns by 2 atan(omega dt/2)
!<  per step; started at the ExB velocity, the particle keeps it; the map v^(n-1/2) -> v^(n+1/2) has determinant 1
!<  (centred differences); the two layouts give the same result bit for bit.
!<+ **Defining equations.** Each pusher's result satisfies its implicit defining equation (selah_pushers) to
!<  round-off, with generic fields and momentum.
!<+ **Force-free (Vay 2008).** gamma_0 = 10 moving along x in B = z with E = -v x B, q B dt/m = 0.1: Vay and
!<  Higuera–Cary keep the momentum to round-off; Boris does not (1.24e-7, the research's reference value).
!<+ **ExB drift, relativistic.** E/cB = 0.5, q B dt/m = 1, started at u_E: Higuera–Cary and Vay keep it; Boris moves it.
!<+ **Volume.** The momentum map's determinant: 1 for relativistic Boris and Higuera–Cary, 1 - 1.22e-4 for Vay at
!<  q dt/m = 1 (the research's reference case).
!<+ **Non-relativistic limit.** With c = 1e6 every relativistic pusher gives Boris's result to 1e-10.
use penf, only : I4P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
character(18), parameter :: NAMES(4) = [character(18) :: 'boris', 'boris-relativistic', 'vay', 'higuera-cary']

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
call boris
call defining_equations
call force_free
call drift
call volume
call nonrelativistic_limit
if (failures > 0) error stop 'selah_test_pushers: failed'
print '(A)', 'selah_test_pushers: ok'

contains
   function push1(pusher, x, u, e, b, dt, c, layout, steps) result(y)
   !< One particle pushed `steps` times (default 1) through device arrays: y = (x, u) after.
   integer(I4P), intent(in)           :: pusher
   real(R_P),    intent(in)           :: x(3), u(3), e(3), b(3), dt
   real(R_P),    intent(in), optional :: c
   integer(I4P), intent(in), optional :: layout, steps
   real(R_P)                          :: y(6)
   real(R_P), pointer                 :: p(:,:)=>null(), pe(:,:)=>null(), pb(:,:)=>null()
   real(R_P)                          :: h(8,1)
   integer(I4P)                       :: ly, ns, k, ierr

   ly = LAYOUT_V_FIRST ; if (present(layout)) ly = layout
   ns = 1 ; if (present(steps)) ns = steps
   h(:,1) = [x, u, 1._R_P, 1._R_P]
   if (ly == LAYOUT_V_FIRST) then
      call dev_alloc(fptr_dev=p, ubounds=[8_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_alloc(fptr_dev=pe, ubounds=[3_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_alloc(fptr_dev=pb, ubounds=[3_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_memcpy_to_device(dst=p, src=h)
      call dev_memcpy_to_device(dst=pe, src=reshape(e, [3, 1]))
      call dev_memcpy_to_device(dst=pb, src=reshape(b, [3, 1]))
   else
      call dev_alloc(fptr_dev=p, ubounds=[1_I4P, 8_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_alloc(fptr_dev=pe, ubounds=[1_I4P, 3_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_alloc(fptr_dev=pb, ubounds=[1_I4P, 3_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_memcpy_to_device(dst=p, src=transpose(h))
      call dev_memcpy_to_device(dst=pe, src=reshape(e, [1, 3]))
      call dev_memcpy_to_device(dst=pb, src=reshape(b, [1, 3]))
   endif
   do k=1, ns
      call selah_push(pusher, ly, 1_I4P, dt, 1_I4P, 4_I4P, 7_I4P, 8_I4P, p, pe, pb, c)
   enddo
   if (ly == LAYOUT_V_FIRST) then
      call dev_memcpy_from_device(dst=h, src=p)
   else
      block
         real(R_P) :: ht(1,8)
         call dev_memcpy_from_device(dst=ht, src=p)
         h = transpose(ht)
      endblock
   endif
   y = h(1:6,1)
   call dev_free(p) ; call dev_free(pe) ; call dev_free(pb)
   endfunction push1

   pure function gam(u) result(g)
   real(R_P), intent(in) :: u(3)
   real(R_P)             :: g
   g = sqrt(1._R_P + sum(u**2))
   endfunction gam

   pure function cross(a, b) result(c)
   real(R_P), intent(in) :: a(3), b(3)
   real(R_P)             :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   endfunction cross

   subroutine boris
   !< Speed and turning angle in pure B; the ExB velocity; the layouts.
   real(R_P)    :: y(6), z(6), w, ang
   integer(I4P) :: k

   print '(A)', 'Boris (non-relativistic):'
   do k=1, 3
      w = merge(0.1_R_P, merge(1._R_P, 3._R_P, k == 2), k == 1)
      y = push1(PUSHER_BORIS, [0._R_P, 0._R_P, 0._R_P], [1._R_P, 0._R_P, 0._R_P], [0._R_P, 0._R_P, 0._R_P], &
                [0._R_P, 0._R_P, w], 1._R_P)
      ang = atan2(-y(5), y(4))
      print '(A,F4.1,A,2ES24.16)', '        omega dt =', w, ': angle, 2 atan(omega dt/2)', ang, 2._R_P*atan(w/2._R_P)
      call check(abs(ang - 2._R_P*atan(w/2._R_P)) <= 8*epsilon(1._R_P), 'boris: turns by 2 atan(omega dt/2)')
   enddo
   y = push1(PUSHER_BORIS, [0._R_P, 0._R_P, 0._R_P], [1._R_P, 0.3_R_P, -0.2_R_P], [0._R_P, 0._R_P, 0._R_P], &
             [0.3_R_P, -0.4_R_P, 1.2_R_P], 0.7_R_P, steps=1000_I4P)
   call check(abs(sqrt(sum(y(4:6)**2)) - sqrt(1.13_R_P)) <= 1.e3_R_P*epsilon(1._R_P), &
              'boris: the speed kept in a pure B over 1000 steps')
   y = push1(PUSHER_BORIS, [0._R_P, 0._R_P, 0._R_P], [0.25_R_P, 0._R_P, 0._R_P], [0._R_P, 0.5_R_P, 0._R_P], &
             [0._R_P, 0._R_P, 2._R_P], 0.5_R_P, steps=16_I4P)
   call check(all(abs(y(4:6) - [0.25_R_P, 0._R_P, 0._R_P]) <= 8*epsilon(1._R_P)), 'boris: the ExB velocity kept')
   y = push1(PUSHER_BORIS, [0.1_R_P, 0.2_R_P, 0.3_R_P], [1._R_P, 0.3_R_P, -0.2_R_P], [0.3_R_P, -0.2_R_P, 0.15_R_P], &
             [0.1_R_P, 0.4_R_P, 1.3_R_P], 0.7_R_P, steps=7_I4P)
   z = push1(PUSHER_BORIS, [0.1_R_P, 0.2_R_P, 0.3_R_P], [1._R_P, 0.3_R_P, -0.2_R_P], [0.3_R_P, -0.2_R_P, 0.15_R_P], &
             [0.1_R_P, 0.4_R_P, 1.3_R_P], 0.7_R_P, layout=LAYOUT_V_LAST, steps=7_I4P)
   call check(all(y == z), 'boris: the (np, ncomp) layout bitwise the (ncomp, np) one')
   endsubroutine boris

   subroutine defining_equations
   !< Each pusher's u^+ against its implicit equation.
   real(R_P) :: e(3), b(3), u(3), y(6), up(3), ub(3), lhs(3), rhs(3), dt, qh
   integer(I4P) :: k

   print '(A)', 'defining equations (relative residual):'
   e = [0.3_R_P, -0.2_R_P, 0.15_R_P] ; b = [0.1_R_P, 0.4_R_P, 1.3_R_P] ; u = [0.7_R_P, 1.9_R_P, -0.4_R_P]
   dt = 0.8_R_P ; qh = 0.5_R_P * dt
   do k=1, 4
      y = push1(k, [0._R_P, 0._R_P, 0._R_P], u, e, b, dt)
      up = y(4:6)
      ub = 0.5_R_P * (up + u)
      lhs = up - u
      select case(k)
      case(PUSHER_BORIS)
         rhs = 2._R_P*qh * (e + cross(ub, b))
      case(PUSHER_BORIS_RELATIVISTIC)
         rhs = 2._R_P*qh * (e + cross(ub, b) / gam(u + qh*e))
      case(PUSHER_VAY)
         rhs = 2._R_P*qh * (e + cross(0.5_R_P*(up/gam(up) + u/gam(u)), b))
      case(PUSHER_HIGUERA_CARY)
         rhs = 2._R_P*qh * (e + cross(ub, b) / gam(ub))
      endselect
      print '(A,ES10.2)', '        '//trim(NAMES(k))//':', maxval(abs(lhs - rhs)) / sqrt(sum(u**2))
      call check(maxval(abs(lhs - rhs)) <= 64*epsilon(1._R_P)*sqrt(sum(u**2)), trim(NAMES(k))//': its defining equation')
   enddo
   endsubroutine defining_equations

   subroutine force_free
   !< gamma_0 = 10 along x, B = z, E = -v x B, q B dt/m = 0.1: one step.
   real(R_P) :: u(3), v(3), e(3), b(3), y(6), d
   integer(I4P) :: k

   print '(A)', 'force-free, gamma_0 = 10, one step (|u^1 - u^0|/|u^0|):'
   u = [sqrt(99._R_P), 0._R_P, 0._R_P] ; v = u / gam(u) ; b = [0._R_P, 0._R_P, 1._R_P] ; e = -cross(v, b)
   do k=2, 4
      y = push1(k, [0._R_P, 0._R_P, 0._R_P], u, e, b, 0.1_R_P)
      d = sqrt(sum((y(4:6) - u)**2)) / sqrt(sum(u**2))
      print '(A,ES10.2)', '        '//trim(NAMES(k))//':', d
      if (k == PUSHER_BORIS_RELATIVISTIC) then
         call check(abs(d / 1.237e-7_R_P - 1._R_P) < 0.05_R_P, trim(NAMES(k))//': not force-free (1.24e-7, reference)')
      else
         call check(d <= 64*epsilon(1._R_P), trim(NAMES(k))//': force-free, the momentum kept')
      endif
   enddo
   endsubroutine force_free

   subroutine drift
   !< E = (0, 0.5, 0), B = z, q B dt/m = 1, from u_E = gamma_E v_E: one step.
   real(R_P) :: u(3), e(3), b(3), y(6), d
   integer(I4P) :: k

   print '(A)', 'relativistic ExB drift, E/cB = 0.5, q B dt/m = 1, one step from u_E (|u^1 - u_E|/|u_E|):'
   e = [0._R_P, 0.5_R_P, 0._R_P] ; b = [0._R_P, 0._R_P, 1._R_P]
   u = [0.5_R_P, 0._R_P, 0._R_P] / sqrt(1._R_P - 0.25_R_P)
   do k=2, 4
      y = push1(k, [0._R_P, 0._R_P, 0._R_P], u, e, b, 1._R_P)
      d = sqrt(sum((y(4:6) - u)**2)) / sqrt(sum(u**2))
      print '(A,ES10.2)', '        '//trim(NAMES(k))//':', d
      if (k == PUSHER_BORIS_RELATIVISTIC) then
         call check(d > 1.e-3_R_P, trim(NAMES(k))//': u_E not kept (the wrong drift velocity)')
      else
         call check(d <= 64*epsilon(1._R_P), trim(NAMES(k))//': u_E kept')
      endif
   enddo
   endsubroutine drift

   subroutine volume
   !< det of the momentum map by centred differences (h = 1e-5): 1 for Boris (both) and HC, 1 - 1.22e-4 for Vay.
   real(R_P) :: u(3), e(3), b(3), jac(3,3), y(6), z(6), d, h
   integer(I4P) :: k, j

   print '(A)', 'volume: det of the momentum map - 1 (q dt/m = 1):'
   e = [0.3_R_P, -0.2_R_P, 0.15_R_P] ; b = [0.1_R_P, 0.4_R_P, 1.3_R_P] ; u = [0.7_R_P, 1.9_R_P, -0.4_R_P]
   h = 1.e-5_R_P
#if defined _R_P_IS_R4P
   h = 1.e-2_R_P
#endif
   do k=1, 4
      do j=1, 3
         y = push1(k, [0._R_P, 0._R_P, 0._R_P], u + h*unit(j), e, b, 1._R_P)
         z = push1(k, [0._R_P, 0._R_P, 0._R_P], u - h*unit(j), e, b, 1._R_P)
         jac(:,j) = (y(4:6) - z(4:6)) / (2._R_P*h)
      enddo
      d = jac(1,1)*(jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2)) - jac(1,2)*(jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1)) + &
          jac(1,3)*(jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1)) - 1._R_P
      print '(A,ES11.3)', '        '//trim(NAMES(k))//':', d
#if !defined _R_P_IS_R4P
      if (k == PUSHER_VAY) then
         call check(abs(d + 1.22e-4_R_P) < 0.02e-4_R_P, trim(NAMES(k))//': not volume-preserving (-1.22e-4, reference)')
      else
         call check(abs(d) < 1.e-8_R_P, trim(NAMES(k))//': volume-preserving')
      endif
#endif
   enddo
   endsubroutine volume

   pure function unit(j) result(v)
   integer(I4P), intent(in) :: j
   real(R_P)                :: v(3)
   v = 0._R_P ; v(j) = 1._R_P
   endfunction unit

   subroutine nonrelativistic_limit
   !< c = 1e6: every relativistic pusher gives Boris's result.
   real(R_P) :: y0(6), y(6)
   integer(I4P) :: k

   y0 = push1(PUSHER_BORIS, [0.1_R_P, 0.2_R_P, 0.3_R_P], [1._R_P, 0.3_R_P, -0.2_R_P], [0.3_R_P, -0.2_R_P, 0.15_R_P], &
              [0.1_R_P, 0.4_R_P, 1.3_R_P], 0.7_R_P, c=1.e6_R_P)
   do k=2, 4
      y = push1(k, [0.1_R_P, 0.2_R_P, 0.3_R_P], [1._R_P, 0.3_R_P, -0.2_R_P], [0.3_R_P, -0.2_R_P, 0.15_R_P], &
                [0.1_R_P, 0.4_R_P, 1.3_R_P], 0.7_R_P, c=1.e6_R_P)
      call check(maxval(abs(y - y0)) <= max(1.e-10_R_P, 1.e3_R_P*epsilon(1._R_P)), &
                 trim(NAMES(k))//': the non-relativistic limit is Boris')
   enddo
   endsubroutine nonrelativistic_limit
endprogram selah_test_pushers
