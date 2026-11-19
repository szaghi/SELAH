!< SELAH test: time-dependent residuals through the stepper, every scheme and form, on the host or the device.

#include "fundal.H"
program selah_test_nonautonomous
!< SELAH test: time-dependent residuals through the stepper, every scheme and form, on the host or the device.
!<
!< The residual of stage k is evaluated at t + stage_time(k) dt, as a caller with time-dependent boundary conditions or
!< sources evaluates it (FOODIE's convergence and polynomial-exactness problems, ported through the stepper):
!<
!<+ **Quadrature (T2).** u_q' = q t^(q-1), u_q(0) = 0, one cell per q = 1, ..., p + 1, three steps to t = 1: a scheme of
!<  order p integrates t^(q-1) exactly for q <= p, so u_q(1) = 1 to round-off; a wrong stage time or weight breaks it
!<  for q >= 2. Every Runge–Kutta scheme in the Butcher form and in its other forms (two-register, 2N, Ketcheson), and
!<  every splitting twice: with the time-dependent residual in partition 2 (time carried by partition 1) and in
!<  partition 1 (time carried by partition 2). The error at q = p + 1 is printed, not asserted (a symmetric
!<  composition gains a degree).
!<+ **Observed order (T1), R8P builds.** The Riccati equation u' = t u² (nonlinear and non-autonomous), four cells
!<  u(0) = 1/4, 1/2, 3/4, 1, exact u = 1/(1/u(0) - t²/2): every Runge–Kutta scheme and form reaches its declared order
!<  (the error the largest over the steps and cells, at four halvings of dt). The splittings, on the forced oscillator
!<  x' = y, y' = -x - 3 sin 2t (exact x = sin 2t) and on its mirror x' = y + 3 sin 2t, y' = -x (exact y = sin 2t):
!<  every sub-flow exact, the declared order with either partition time-dependent.
use penf, only : I4P, R8P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures, observed_orders
implicit none
character(7),  parameter :: OTHER_SCHEMES(6) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, &
                                                SCHEME_LSRK_33, SCHEME_LSRK_54, SCHEME_SSP_104]
character(14), parameter :: OTHER_FORMS(6) = [character(14) :: FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, &
                                              FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_WILLIAMSON_2N, &
                                              FORM_KETCHESON_2R]
real(R_P),     parameter :: TOL_EXACT = 1000._R_P * epsilon(1._R_P) !< Quadrature exactness, relative to u_q(1) = 1.
integer(I4P) :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
#endif
print '(A)', 'quadrature, Runge–Kutta (u_q'' = q t^(q-1), exact for q <= p):'
do s=1, size(RK_SCHEMES)
   call quadrature_rk(trim(RK_SCHEMES(s)), FORM_BUTCHER)
enddo
do s=1, size(OTHER_SCHEMES)
   call quadrature_rk(trim(OTHER_SCHEMES(s)), trim(OTHER_FORMS(s)))
enddo
print '(A)', 'quadrature, splittings (the time-dependent residual in partition 2, then in partition 1):'
do s=1, size(SPLITTING_SCHEMES)
   call quadrature_splitting(trim(SPLITTING_SCHEMES(s)), 2_I4P)
   call quadrature_splitting(trim(SPLITTING_SCHEMES(s)), 1_I4P)
enddo
#if defined _R_P_IS_R4P
print '(A)', '  skip  observed orders (R4P coefficients cannot show the order of the exact ones)'
#else
print '(A)', 'observed order, Runge–Kutta (Riccati u'' = t u^2):'
do s=1, size(RK_SCHEMES)
   call order_rk(trim(RK_SCHEMES(s)), FORM_BUTCHER)
enddo
do s=1, size(OTHER_SCHEMES)
   call order_rk(trim(OTHER_SCHEMES(s)), trim(OTHER_FORMS(s)))
enddo
print '(A)', 'observed order, splittings (the forced oscillator, either partition time-dependent):'
do s=1, size(SPLITTING_SCHEMES)
   call order_splitting(trim(SPLITTING_SCHEMES(s)), 2_I4P)
   call order_splitting(trim(SPLITTING_SCHEMES(s)), 1_I4P)
enddo
#endif

if (failures > 0) error stop 'selah_test_nonautonomous: failed'
print '(A)', 'selah_test_nonautonomous: ok'

contains
   subroutine quadrature_rk(name, form)
   !< u_q' = q t^(q-1), q = 1, ..., p + 1 (one cell each), three steps to t = 1 through the stepper.
   character(*), intent(in) :: name !< Scheme.
   character(*), intent(in) :: form !< Storage form.
   type(butcher_tableau)    :: bt
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: u(:)=>null(), r(:)
   real(R_P), allocatable   :: h(:)
   real(R_P)                :: dt, t, tau
   integer(I4P)             :: p, m, step, k, ierr

   bt = butcher(name)
   p = bt%order
   m = p + 1_I4P
   allocate(h(m)) ; h = 0._R_P
   call dev_alloc(fptr_dev=u, ubounds=[m], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=u, src=h)
   bx = box_interior(lb=[1_I4P], ub=[m], ngc=0_I4P, components=.false.)
   call st%initialize(name, form, lb=[1_I4P], ub=[m], components=.false.)
   dt = 1._R_P / 3._R_P ; t = 0._R_P
   do step=1, 3
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         tau = t + st%stage_time(k) * dt
         call power_residual(m, tau, r)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
      t = step * dt
   enddo
   call dev_memcpy_from_device(dst=h, src=u)
   call st%destroy
   call dev_free(u)
   print '(A,ES10.2,A,ES10.2)', '        '//name//' '//form//': max |u_q(1) - 1|, q <= p:', maxval(abs(h(1:p) - 1._R_P)), &
                                '; q = p + 1:', abs(h(m) - 1._R_P)
   call check(maxval(abs(h(1:p) - 1._R_P)) <= TOL_EXACT, name//' '//form//': t^(q-1) integrated exactly for q <= p')
   endsubroutine quadrature_rk

   subroutine quadrature_splitting(name, tp)
   !< u_q' = q t^(q-1) in partition tp (component tp), 0 in the other, q = 1, ..., p + 1 (one cell each).
   character(*), intent(in) :: name !< Scheme.
   integer(I4P), intent(in) :: tp   !< Partition whose residual depends on time.
   type(splitting_composition) :: sc
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: u(:,:)=>null(), r(:,:)
   real(R_P), allocatable   :: h(:,:)
   real(R_P)                :: dt, t, tau
   integer(I4P)             :: p, m, step, k, ierr
   character(1)             :: tag

   sc = splitting(name)
   p = sc%order
   m = p + 1_I4P
   write(tag, '(I1)') tp
   allocate(h(2,m)) ; h = 0._R_P
   call dev_alloc(fptr_dev=u, ubounds=[2_I4P, m], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=u, src=h)
   bx = box_interior(lb=[1_I4P, 1_I4P], ub=[2_I4P, m], ngc=0_I4P, components=.true.)
   call st%initialize(name, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[2_I4P, m], components=.true., &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   dt = 1._R_P / 3._R_P ; t = 0._R_P
   do step=1, 3
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         tau = t + st%stage_time(k) * dt
         call split_power_residual(m, tau, merge(1._R_P, 0._R_P, st%stage_partition(k) == tp), st%stage_partition(k), r)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
      t = step * dt
   enddo
   call dev_memcpy_from_device(dst=h, src=u)
   call st%destroy
   call dev_free(u)
   print '(A,ES10.2,A,ES10.2)', '        '//name//', partition '//tag//' time-dependent: max |u_q(1) - 1|, q <= p:', &
                                maxval(abs(h(tp,1:p) - 1._R_P)), '; q = p + 1:', abs(h(tp,m) - 1._R_P)
   call check(maxval(abs(h(tp,1:p) - 1._R_P)) <= TOL_EXACT .and. all(h(3-tp,:) == 0._R_P), &
              name//', partition '//tag//' time-dependent: t^(q-1) integrated exactly for q <= p')
   endsubroutine quadrature_splitting

   subroutine power_residual(m, tau, r)
   !< r_q = q tau^(q-1), q = 1, ..., m.
   integer(I4P), intent(in)    :: m   !< Cells.
   real(R_P),    intent(in)    :: tau !< Time of the stage.
   real(R_P),    intent(inout) :: r(1:)
   integer(I4P)                :: q

   !$acc parallel loop independent DEVICEVAR(r)
   !$omp OMPLOOP DEVICEPTR(r)
   do q=1, m
      r(q) = q * tau**(q - 1)
   enddo
   endsubroutine power_residual

   subroutine split_power_residual(m, tau, on, part, r)
   !< The residual of partition part (component part): on q tau^(q-1), q = 1, ..., m (on = 1 for the time-dependent
   !< partition, 0 for the other).
   integer(I4P), intent(in)    :: m    !< Cells.
   real(R_P),    intent(in)    :: tau  !< Time of the flow.
   real(R_P),    intent(in)    :: on   !< 1 or 0.
   integer(I4P), intent(in)    :: part !< Partition of the flow.
   real(R_P),    intent(inout) :: r(1:,1:)
   integer(I4P)                :: q

   !$acc parallel loop independent DEVICEVAR(r)
   !$omp OMPLOOP DEVICEPTR(r)
   do q=1, m
      r(part,q) = on * q * tau**(q - 1)
   enddo
   endsubroutine split_power_residual

   subroutine order_rk(name, form)
   !< The observed order on the Riccati equation, four halvings of dt.
   character(*), intent(in) :: name !< Scheme.
   character(*), intent(in) :: form !< Storage form.
   type(butcher_tableau)    :: bt
   real(R8P)                :: errors(4), orders(3)
   integer(I4P)             :: p, n0, i

   bt = butcher(name)
   p = bt%order
   n0 = merge(4_I4P, 16_I4P, p >= 3)
   do i=1, 4
      errors(i) = riccati_error(name, form, n0 * 2_I4P**(i-1))
   enddo
   orders = observed_orders(errors)
   print '(A,3F7.3,A,ES9.2)', '        '//name//' '//form//' observed orders:', orders, '  finest error', errors(4)
   call check(orders(3) >= real(p, R8P) - 0.2_R8P, name//' '//form//': the declared order on the Riccati equation')
   endsubroutine order_rk

   function riccati_error(name, form, n) result(err)
   !< max over the steps and cells of |u - exact| for n steps of u' = t u² to t = 1.
   character(*), intent(in) :: name, form
   integer(I4P), intent(in) :: n   !< Steps.
   real(R8P)                :: err !< Error.
   integer(I4P), parameter  :: M = 4
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: u(:)=>null(), us(:), r(:)
   real(R_P)                :: h(M), u0(M), dt, t, tau
   integer(I4P)             :: step, k, j, ierr

   u0 = [(0.25_R_P * j, j=1, M)]
   h = u0
   call dev_alloc(fptr_dev=u, ubounds=[M], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=u, src=h)
   bx = box_interior(lb=[1_I4P], ub=[M], ngc=0_I4P, components=.false.)
   call st%initialize(name, form, lb=[1_I4P], ub=[M], components=.false.)
   dt = 1._R_P / real(n, R_P) ; t = 0._R_P ; err = 0._R8P
   do step=1, n
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         if (st%stage_state_is_q(k)) then
            us => u
         else
            call st%stage_state_view(us)
         endif
         call st%rhs_register_view(r)
         tau = t + st%stage_time(k) * dt
         call riccati_residual(M, tau, us, r)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
      t = step * dt
      call dev_memcpy_from_device(dst=h, src=u)
      err = max(err, real(maxval(abs(h - 1._R_P / (1._R_P / u0 - t**2 / 2._R_P))), R8P))
   enddo
   call st%destroy
   call dev_free(u)
   endfunction riccati_error

   subroutine riccati_residual(m, tau, us, r)
   !< r = tau us².
   integer(I4P), intent(in)    :: m   !< Cells.
   real(R_P),    intent(in)    :: tau !< Time of the stage.
   real(R_P),    intent(in)    :: us(1:)
   real(R_P),    intent(inout) :: r(1:)
   integer(I4P)                :: j

   !$acc parallel loop independent DEVICEVAR(us, r)
   !$omp OMPLOOP DEVICEPTR(us, r)
   do j=1, m
      r(j) = tau * us(j)**2
   enddo
   endsubroutine riccati_residual

   subroutine order_splitting(name, tp)
   !< The observed order on the forced oscillator with partition tp time-dependent, four halvings of dt.
   character(*), intent(in) :: name !< Scheme.
   integer(I4P), intent(in) :: tp   !< Partition whose residual depends on time.
   type(splitting_composition) :: sc
   real(R8P)                :: errors(4), orders(3)
   integer(I4P)             :: p, n0, i
   character(1)             :: tag

   sc = splitting(name)
   p = sc%order
   write(tag, '(I1)') tp
   n0 = merge(4_I4P, 16_I4P, p >= 4)
   if (p >= 6) n0 = 2_I4P   ! order 6 reaches the R8P floor sooner: its finest step keeps the error near 1e-11
   do i=1, 4
      errors(i) = forced_error(name, tp, n0 * 2_I4P**(i-1))
   enddo
   orders = observed_orders(errors)
   print '(A,3F7.3,A,ES9.2)', '        '//name//', partition '//tag//' time-dependent, observed orders:', orders, &
                              '  finest error', errors(4)
   call check(orders(3) >= real(p, R8P) - 0.2_R8P, name//', partition '//tag//' time-dependent: the declared order')
   endsubroutine order_splitting

   function forced_error(name, tp, n) result(err)
   !< max over the steps of |(x, y) - exact| for n steps to t = 1 of the forced oscillator, partition tp forced:
   !< tp = 2: x' = y, y' = -x - 3 sin 2t, (x, y)(0) = (0, 2), exact (sin 2t, 2 cos 2t);
   !< tp = 1: x' = y + 3 sin 2t, y' = -x, (x, y)(0) = (-2, 0), exact (-2 cos 2t, sin 2t).
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: tp  !< Forced partition.
   integer(I4P), intent(in) :: n   !< Steps.
   real(R8P)                :: err !< Error.
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: u(:,:)=>null(), r(:,:)
   real(R_P)                :: h(2,1), ex(2), dt, t, tau
   integer(I4P)             :: step, k, ierr

   h(:,1) = merge([0._R_P, 2._R_P], [-2._R_P, 0._R_P], tp == 2)
   call dev_alloc(fptr_dev=u, ubounds=[2_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=u, src=h)
   bx = box_interior(lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], ngc=0_I4P, components=.true.)
   call st%initialize(name, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], components=.true., &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   dt = 1._R_P / real(n, R_P) ; t = 0._R_P ; err = 0._R8P
   do step=1, n
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         tau = t + st%stage_time(k) * dt
         call forced_residual(tp, st%stage_partition(k), tau, u, r)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
      t = step * dt
      call dev_memcpy_from_device(dst=h, src=u)
      ex = merge([sin(2*t), 2*cos(2*t)], [-2*cos(2*t), sin(2*t)], tp == 2)
      err = max(err, real(maxval(abs(h(:,1) - ex)), R8P))
   enddo
   call st%destroy
   call dev_free(u)
   endfunction forced_error

   subroutine forced_residual(tp, part, tau, u, r)
   !< The residual of partition part: x' = y (+ 3 sin 2tau if tp = 1), y' = -x (- 3 sin 2tau if tp = 2).
   integer(I4P), intent(in)    :: tp, part
   real(R_P),    intent(in)    :: tau
   real(R_P),    intent(in)    :: u(1:,1:)
   real(R_P),    intent(inout) :: r(1:,1:)
   real(R_P)                   :: f
   integer(I4P)                :: j

   f = 3._R_P * sin(2._R_P * tau)
   !$acc parallel loop independent DEVICEVAR(u, r)
   !$omp OMPLOOP DEVICEPTR(u, r)
   do j=1, 1
      if (part == 1) then
         r(1,j) = u(2,j) + merge(f, 0._R_P, tp == 1)
      else
         r(2,j) = -u(1,j) - merge(f, 0._R_P, tp == 2)
      endif
   enddo
   endsubroutine forced_residual
endprogram selah_test_nonautonomous
