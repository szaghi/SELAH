!< SELAH test: a companion with an exactly integrated linear decay (a PML's auxiliary fields), on the host or the device.

#include "fundal.H"
program selah_test_decay
!< SELAH test: a companion with an exactly integrated linear decay (a PML's auxiliary fields), on the host or the device.
!<
!< psi' = g - lambda psi with lambda >= 0 per element: the caller's residual is g, the stepper integrates the decay
!< (`decay=lambda` at begin_stage and finish_step; Butcher form): the state of stage i is
!< e^(c_i z) psi^n + phi_1(c_i z) dt sum_j a_ij g_j, z = -lambda dt (notes/companions-pushers-pml.md §5, variant F).
!<
!<+ **Frozen source.** g constant: psi^n = g/lambda + (psi^0 - g/lambda) e^(-lambda t^n) to 1e-13 relative, for lambda dt
!<  = 0.1, 1, 10, 1e3, 1e6, every Runge–Kutta scheme (the explicit stages blow up for lambda dt above the real-axis
!<  bound).
!<+ **lambda = 0.** Bitwise the run without `decay` (every scheme); so is lambda the smallest subnormal, whose
!<  lambda dt/2 rounds to 0 (the kernels' guard: never 0/0).
!<+ **Coupled, two-way.** A field u and its auxiliary psi, the PML's structure: u' = -psi + J, psi' = sigma u - lambda psi
!<  (u: the field stepper; psi: a companion with decay). (i) The static equilibrium u* = lambda J / sigma, psi* = J is
!<  kept to round-off at lambda dt = 50; (ii) the observed order with a smooth start is 2 (1 for euler; the variant's
!<  order, notes §5) against a fine reference; (iii) at lambda dt = 50 the solution stays bounded where the explicit
!<  stages of the same scheme blow up.
use penf, only : I4P, R8P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures, observed_orders
implicit none
integer(I4P) :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'frozen source, lambda dt = 0.1, 1, 10, 1e3, 1e6 (largest relative error over 8 steps):'
do s=1, size(RK_SCHEMES)
   call frozen(trim(RK_SCHEMES(s)))
enddo
print '(A)', 'coupled field and auxiliary field:'
do s=1, size(RK_SCHEMES)
   call coupled(trim(RK_SCHEMES(s)))
enddo
if (failures > 0) error stop 'selah_test_decay: failed'
print '(A)', 'selah_test_decay: ok'

contains
   subroutine frozen(name)
   !< 5 cells, lambda dt = 0.1 ... 1e6, g = 1, psi^0 = 3: 8 steps against the exact solution; lambda = 0 bitwise.
   character(*), intent(in) :: name
   real(R_P),    parameter  :: LDT(5) = [0.1_R_P, 1._R_P, 10._R_P, 1.e3_R_P, 1.e6_R_P], G = 1._R_P, P0 = 3._R_P
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q(:)=>null(), lam(:)=>null(), r(:), q0(:)=>null()
   real(R_P)                :: h(5), hl(5), ex(5), worst, dt
   integer(I4P)             :: step, k, ierr

   dt = 1._R_P
   call dev_alloc(fptr_dev=q, ubounds=[5_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=q0, ubounds=[5_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=lam, ubounds=[5_I4P], lbounds=[1_I4P], ierr=ierr)
   bx = box_interior(lb=[1_I4P], ub=[5_I4P], ngc=0_I4P, components=.false.)
   call st%initialize(name, FORM_BUTCHER, lb=[1_I4P], ub=[5_I4P], components=.false.)
   hl = LDT / dt ; h = P0
   call dev_memcpy_to_device(dst=lam, src=hl) ; call dev_memcpy_to_device(dst=q, src=h)
   worst = 0._R_P
   do step=1, 8
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q, decay=lam)
         call st%rhs_register_view(r)
         call constant(G, r)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q, decay=lam)
      call dev_memcpy_from_device(dst=h, src=q)
      ex = G/hl + (P0 - G/hl) * exp(-hl * step * dt)
      worst = max(worst, maxval(abs(h - ex) / abs(ex)))
   enddo
   print '(A,ES10.2)', '        '//name//':', worst
   call check(worst <= 1.e-13_R_P .or. (worst <= 1.e3_R_P*epsilon(1._R_P)), name//': the frozen source exact for any lambda dt')
   ! lambda = 0: bitwise the run without decay
   hl = 0._R_P ; h = P0
   call dev_memcpy_to_device(dst=lam, src=hl) ; call dev_memcpy_to_device(dst=q, src=h)
   call dev_memcpy_to_device(dst=q0, src=h)
   do step=1, 3
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q, decay=lam)
         call st%rhs_register_view(r)
         call constant(G, r)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q, decay=lam)
      call st%begin_step(bx, q0)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q0)
         call st%rhs_register_view(r)
         call constant(G, r)
         call st%end_stage(k, bx, dt, q0)
      enddo
      call st%finish_step(bx, dt, q0)
   enddo
   call dev_memcpy_from_device(dst=h, src=q) ; call dev_memcpy_from_device(dst=ex, src=q0)
   call check(all(h == ex), name//': lambda = 0 bitwise the explicit step')
   ! lambda the smallest subnormal: lambda dt / 2 rounds to 0, the kernels take the explicit branch (never 0/0)
   hl = nearest(0._R_P, 1._R_P) ; h = P0
   call dev_memcpy_to_device(dst=lam, src=hl) ; call dev_memcpy_to_device(dst=q, src=h)
   do step=1, 3
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q, decay=lam)
         call st%rhs_register_view(r)
         call constant(G, r)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q, decay=lam)
   enddo
   call dev_memcpy_from_device(dst=h, src=q)
   call check(all(h == ex), name//': lambda subnormal, bitwise the explicit step (no NaN)')
   call st%destroy
   call dev_free(q) ; call dev_free(q0) ; call dev_free(lam)
   endsubroutine frozen

   subroutine constant(g, r)
   !< r = g everywhere.
   real(R_P), intent(in)    :: g
   real(R_P), intent(inout) :: r(1:)
   integer(I4P)             :: i

   !$acc parallel loop independent DEVICEVAR(r)
   !$omp OMPLOOP DEVICEPTR(r)
   do i=1, size(r)
      r(i) = g
   enddo
   endsubroutine constant

   subroutine coupled(name)
   !< u' = -psi + J, psi' = sigma u - lambda psi: the equilibrium, the order, the stability.
   character(*), intent(in) :: name
   real(R_P),    parameter  :: SIG = 3._R_P, J = 0.5_R_P
   real(R_P)                :: err(4), y(2), lam, o(3), ref(2)
   integer(I4P)             :: l, n
   type(butcher_tableau)    :: bt

   bt = butcher(name)
   ! the static equilibrium at lambda dt = 50
   lam = 50._R_P
   y = run(name, lam*J/SIG, J, SIG, lam, J, 1._R_P, 4_I4P, .true.)
   call check(all(abs(y - [lam*J/SIG, J]) <= 64*epsilon(1._R_P)*[lam*J/SIG, J]), &
              name//': the coupled static equilibrium kept at lambda dt = 50')
   ! stability at lambda dt = 50, from a perturbed state: bounded; the explicit stages blow up
   y = run(name, 1._R_P, 0._R_P, SIG, lam, 0._R_P, 1._R_P, 20_I4P, .true.)
   call check(all(abs(y) < 10._R_P), name//': bounded at lambda dt = 50 with the decay integrated')
   y = run(name, 1._R_P, 0._R_P, SIG, lam, 0._R_P, 1._R_P, 20_I4P, .false.)
   call check(.not. all(abs(y) < 10._R_P), name//': the explicit stages blow up at lambda dt = 50')
#if !defined _R_P_IS_R4P
   ! the order, lambda = 2, smooth start, t = 1
   lam = 2._R_P
   ref = reference(1._R_P, 0._R_P, SIG, lam, J)
   do l=1, 4
      n = 16 * 2**(l-1)
      y = run(name, 1._R_P, 0._R_P, SIG, lam, J, 1._R_P/n, n, .true.)
      err(l) = maxval(abs(y - ref))
   enddo
   o = observed_orders(real(err, R8P))
   print '(A,4ES10.2,A,3F6.2)', '        '//name//': err', err, '  orders', o
   call check(all(abs(o(2:) - min(bt%order, 2_I4P)) <= 0.15_R_P), name//': order min(p, 2), coupled')
#endif
   endsubroutine coupled

   function run(name, u0, p0, sig, lam, j, dt, n, exponential) result(y)
   !< n steps of the coupled pair from (u0, p0): u the field, psi a companion (decay = lam when exponential, else
   !< -lam psi in its residual).
   character(*), intent(in) :: name
   real(R_P),    intent(in) :: u0, p0, sig, lam, j, dt
   integer(I4P), intent(in) :: n
   logical,      intent(in) :: exponential
   real(R_P)                :: y(2)
   type(stepper)            :: st, sp
   type(box)                :: bx
   real(R_P), pointer       :: u(:)=>null(), p(:)=>null(), dl(:)=>null(), uu(:), pp(:), ru(:), rp(:)
   real(R_P)                :: h(1)
   integer(I4P)             :: step, k, ierr

   call dev_alloc(fptr_dev=u, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=p, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=dl, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=u, src=[u0]) ; call dev_memcpy_to_device(dst=p, src=[p0])
   call dev_memcpy_to_device(dst=dl, src=[lam])
   bx = box_interior(lb=[1_I4P], ub=[1_I4P], ngc=0_I4P, components=.false.)
   call st%initialize(name, FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
   call sp%initialize_companion(st, lb=[1_I4P], ub=[1_I4P], components=.false.)
   do step=1, n
      call st%begin_step(bx, u) ; call sp%begin_step(bx, p)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         if (exponential) then
            call sp%begin_stage(k, bx, dt, p, decay=dl)
         else
            call sp%begin_stage(k, bx, dt, p)
         endif
         if (st%stage_state_is_q(k)) then ; uu => u ; else ; call st%stage_state_view(uu) ; endif
         if (sp%stage_state_is_q(k)) then ; pp => p ; else ; call sp%stage_state_view(pp) ; endif
         call st%rhs_register_view(ru) ; call sp%rhs_register_view(rp)
         call pml_residual(sig, merge(0._R_P, lam, exponential), j, uu, pp, ru, rp)
         call st%end_stage(k, bx, dt, u) ; call sp%end_stage(k, bx, dt, p)
      enddo
      call st%finish_step(bx, dt, u)
      if (exponential) then
         call sp%finish_step(bx, dt, p, decay=dl)
      else
         call sp%finish_step(bx, dt, p)
         ! the explicit stages' blow-up stops 12 decades below overflow, above one step's growth (|R(-50)| <= 1.1e8,
         ! ssp-104): R4P overflowed otherwise (ssp-22: 1200^13); R8P never reaches it
         call dev_memcpy_from_device(dst=h, src=u) ; y(1) = h(1)
         call dev_memcpy_from_device(dst=h, src=p) ; y(2) = h(1)
         if (maxval(abs(y)) > 1.e-12_R_P * huge(1._R_P)) exit
      endif
   enddo
   call dev_memcpy_from_device(dst=h, src=u) ; y(1) = h(1)
   call dev_memcpy_from_device(dst=h, src=p) ; y(2) = h(1)
   call st%destroy ; call sp%destroy
   call dev_free(u) ; call dev_free(p) ; call dev_free(dl)
   endfunction run

   subroutine pml_residual(sig, lam, j, u, p, ru, rp)
   !< ru = -psi + J; rp = sigma u - lam psi (lam = 0 when the stepper integrates the decay).
   real(R_P), intent(in)    :: sig, lam, j
   real(R_P), intent(in)    :: u(1:), p(1:)
   real(R_P), intent(inout) :: ru(1:), rp(1:)

   !$acc serial DEVICEVAR(u, p, ru, rp)
   ru(1) = -p(1) + j
   rp(1) = sig * u(1) - lam * p(1)
   !$acc end serial
   endsubroutine pml_residual

   function reference(u0, p0, sig, lam, j) result(y)
   !< (u, psi) at t = 1 by rk4, 100000 steps.
   real(R_P), intent(in) :: u0, p0, sig, lam, j
   real(R_P)             :: y(2), k1(2), k2(2), k3(2), k4(2), dt
   integer(I4P)          :: step

   y = [u0, p0] ; dt = 1._R_P / 100000._R_P
   do step=1, 100000
      k1 = f(y, sig, lam, j) ; k2 = f(y + 0.5_R_P*dt*k1, sig, lam, j) ; k3 = f(y + 0.5_R_P*dt*k2, sig, lam, j)
      k4 = f(y + dt*k3, sig, lam, j)
      y = y + dt * (k1 + 2._R_P*k2 + 2._R_P*k3 + k4) / 6._R_P
   enddo
   endfunction reference

   pure function f(z, sig, lam, j) result(d)
   !< The coupled pair's right-hand side.
   real(R_P), intent(in) :: z(2), sig, lam, j
   real(R_P)             :: d(2)

   d = [-z(2) + j, sig * z(1) - lam * z(2)]
   endfunction f
endprogram selah_test_decay
