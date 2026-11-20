!< SELAH test: additive splitting of a field step and a caller's flow, on the host or the device.

#include "fundal.H"
program selah_test_additive_splitting
!< SELAH test: additive splitting of a field step and a caller's flow, on the host or the device.
!<
!< u' = F(u) + G(u) with F the rotation x' = y, y' = -x and G the decay x' = -x, y' = -2 y (they do not commute), exact
!< solution exp((A + B) t) u_0 in closed form.
!<
!<+ **Composition (exact sub-flows).** Both flows exact: the observed order of every scheme (lie 1, strang 2, with
!<  either operator first; yoshida-4 4, blanes-moan-s6 4, blanes-moan-s10 6, declared reversible), the queries
!<  (substeps, fractions summing to 1 per operator, `backward()` false for lie and strang only, substep times). Without
!<  `reversible` those three stop the program: selah_test_xfail_additive_backward.
!<+ **The field step through the stepper.** F integrated by the stepper (rk4, Butcher form; ssp-54, Butcher form), G
!<  the exact decay: lie and strang reach orders 1 and 2.
!<+ **Time-dependent terms, both.** F + (sin 2t, 0) and G with the rate (1 + t)(1, 2): each operator integrates its
!<  own time over its substep (the stepper through stage_time, the decay exactly from the substep's start time), from
!<  substep_time: lie 1, strang 2, against a fine reference (rk4 on the whole system).
!<+ **A stiff flow (the cooling model).** u' = a u - k u^2: the field step integrates a u (rk4), the flow is the exact
!<  solution of u' = -k u^2 (Townsend's exact cooling has this form), the exact solution the logistic curve. With
!<  k u_0 from 0.01 to 1000 (k u_0 dt from 1e-3 to 125): lie and strang keep their orders 1 and 2 in the whole range,
!<  the error constant bounded as k grows (no order reduction with this exact, dissipative flow; measured, not
!<  assumed).
use penf, only : I4P, R8P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures, observed_orders
implicit none
real(R_P), parameter :: U0(2) = [1._R_P, 0.5_R_P]   !< Initial state.
real(R_P), parameter :: T_END = 1._R_P              !< Final time.
integer(I4P)         :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'queries:'
call queries
#if defined _R_P_IS_R4P
print '(A)', '  skip  observed orders (R4P)'
#else
print '(A)', 'composition, exact sub-flows (observed order):'
do s=1, size(ADDITIVE_SCHEMES)
   call composition(trim(ADDITIVE_SCHEMES(s)), OPERATOR_FLOW)
   call composition(trim(ADDITIVE_SCHEMES(s)), OPERATOR_FIELD)
enddo
print '(A)', 'the field step through the stepper:'
call with_stepper(SCHEME_LIE, SCHEME_RK4, .false.)
call with_stepper(SCHEME_STRANG, SCHEME_RK4, .false.)
call with_stepper(SCHEME_LIE, SCHEME_SSP_54, .false.)
call with_stepper(SCHEME_STRANG, SCHEME_SSP_54, .false.)
print '(A)', 'both terms time-dependent:'
call with_stepper(SCHEME_LIE, SCHEME_RK4, .true.)
call with_stepper(SCHEME_STRANG, SCHEME_RK4, .true.)
print '(A)', 'a stiff flow (u'' = a u - k u^2):'
call stiff(SCHEME_LIE, OPERATOR_FIELD)
call stiff(SCHEME_STRANG, OPERATOR_FLOW)
call stiff(SCHEME_STRANG, OPERATOR_FIELD)
#endif
if (failures > 0) error stop 'selah_test_additive_splitting: failed'
print '(A)', 'selah_test_additive_splitting: ok'

contains
   subroutine queries
   !< Substeps, fractions, times, backward substeps.
   type(additive_splitting) :: sp
   integer(I4P)             :: m

   call sp%initialize(SCHEME_LIE)
   call check(sp%substeps() == 2 .and. sp%order() == 1 .and. sp%substep_operator(1) == OPERATOR_FIELD .and. &
              sp%substep_operator(2) == OPERATOR_FLOW .and. .not. sp%backward(), 'lie: field then flow, order 1')
   call check(all([(sp%substep_time(m), m=1, 2)] == 0._R_P), 'lie: both substeps start at t')
   call sp%initialize(SCHEME_STRANG)
   call check(sp%substeps() == 3 .and. sp%order() == 2 .and. sp%substep_operator(1) == OPERATOR_FLOW .and. &
              all([(sp%substep_fraction(m), m=1, 3)] == [0.5_R_P, 1._R_P, 0.5_R_P]), &
              'strang: flow h/2, field h, flow h/2')
   call check(all([(sp%substep_time(m), m=1, 3)] == [0._R_P, 0._R_P, 0.5_R_P]), &
              'strang: substep times 0, 0, 1/2 (the same operator''s fractions before)')
   do s=1, size(ADDITIVE_SCHEMES)
      call sp%initialize(trim(ADDITIVE_SCHEMES(s)), reversible=.true.)
      call check(abs(sum(sp%sc%coefficient, mask=sp%sc%partition == OPERATOR_FIELD) - 1._R_P) <= 100*epsilon(1._R_P) &
                 .and. abs(sum(sp%sc%coefficient, mask=sp%sc%partition == OPERATOR_FLOW) - 1._R_P) <= &
                 100*epsilon(1._R_P), trim(ADDITIVE_SCHEMES(s))//': each operator covers the step')
      call check(sp%backward() .eqv. (s > 2), trim(ADDITIVE_SCHEMES(s))//': backward substeps above order 2 only')
   enddo
   endsubroutine queries

   pure function exact(t) result(u)
   !< exp((A + B) t) U0, A + B = [[-1, 1], [-1, -2]]: eigenvalues -3/2 +- i sqrt(3)/2.
   real(R_P), intent(in) :: t
   real(R_P)             :: u(2)
   real(R_P)             :: w, m(2,2)

   w = sqrt(3._R_P) / 2._R_P
   m = reshape([-1._R_P + 1.5_R_P, -1._R_P, 1._R_P, -2._R_P + 1.5_R_P], [2, 2])
   u = exp(-1.5_R_P*t) * (cos(w*t) * U0 + sin(w*t) / w * matmul(m, U0))
   endfunction exact

   pure function rotate(u, h) result(v)
   !< The exact flow of F: x' = y, y' = -x.
   real(R_P), intent(in) :: u(2), h
   real(R_P)             :: v(2)

   v = [cos(h)*u(1) + sin(h)*u(2), -sin(h)*u(1) + cos(h)*u(2)]
   endfunction rotate

   pure function decay(u, h, t0, timed) result(v)
   !< The exact flow of G: x' = -r x, y' = -2 r y, r = 1 (or 1 + t, from t0, when timed).
   real(R_P), intent(in) :: u(2), h, t0
   logical,   intent(in) :: timed
   real(R_P)             :: v(2), g

   g = h
   if (timed) g = h + 0.5_R_P * ((t0 + h)**2 - t0**2)
   v = [exp(-g)*u(1), exp(-2._R_P*g)*u(2)]
   endfunction decay

   subroutine composition(name, first)
   !< Both flows exact: the observed order at four halvings of dt.
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: first
   type(additive_splitting) :: sp
   real(R_P)                :: u(2), err(5), h
   integer(I4P)             :: l, n, step, m

   call sp%initialize(name, first=first, reversible=.true.)
   do l=1, 5
      n = 4 * 2**(l-1)
      u = U0
      do step=1, n
         do m=1, sp%substeps()
            h = sp%substep_fraction(m) * T_END / n
            if (sp%substep_operator(m) == OPERATOR_FIELD) then
               u = rotate(u, h)
            else
               u = decay(u, h, 0._R_P, .false.)
            endif
         enddo
      enddo
      err(l) = maxval(abs(u - exact(T_END)))
   enddo
   call report(name//merge(' (flow first) ', ' (field first)', first == OPERATOR_FLOW), err, sp%order())
   endsubroutine composition

   subroutine with_stepper(name, rk, timed)
   !< F through the stepper (rk, Butcher form, on q_d), G the exact decay; timed: both terms time-dependent.
   character(*), intent(in) :: name, rk
   logical,      intent(in) :: timed
   type(additive_splitting) :: sp
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q(:)=>null(), u(:), r(:)
   real(R_P)                :: h(2), err(5), dt, t, tm, sub
   integer(I4P)             :: l, n, step, m, k, ierr

   call sp%initialize(name)
   call dev_alloc(fptr_dev=q, ubounds=[2_I4P], lbounds=[1_I4P], ierr=ierr)
   bx = box_interior(lb=[1_I4P], ub=[2_I4P], ngc=0_I4P, components=.false.)
   call st%initialize(rk, FORM_BUTCHER, lb=[1_I4P], ub=[2_I4P], components=.false.)
   do l=1, 5
      n = 4 * 2**(l-1)
      dt = T_END / n
      h = U0
      call dev_memcpy_to_device(dst=q, src=h)
      do step=1, n
         t = (step - 1) * dt
         do m=1, sp%substeps()
            sub = sp%substep_fraction(m) * dt
            tm = t + sp%substep_time(m) * dt
            if (sp%substep_operator(m) == OPERATOR_FIELD) then
               call st%begin_step(bx, q)
               do k=1, st%stages()
                  call st%begin_stage(k, bx, sub, q)
                  if (st%stage_state_is_q(k)) then
                     u => q
                  else
                     call st%stage_state_view(u)
                  endif
                  call st%rhs_register_view(r)
                  call rotation_residual(u, r, merge(sin(2._R_P*(tm + st%stage_time(k)*sub)), 0._R_P, timed))
                  call st%end_stage(k, bx, sub, q)
               enddo
               call st%finish_step(bx, sub, q)
            else
               call dev_memcpy_from_device(dst=h, src=q)
               h = decay(h, sub, tm, timed)
               call dev_memcpy_to_device(dst=q, src=h)
            endif
         enddo
      enddo
      call dev_memcpy_from_device(dst=h, src=q)
      if (timed) then
         err(l) = maxval(abs(h - reference_timed()))
      else
         err(l) = maxval(abs(h - exact(T_END)))
      endif
   enddo
   call st%destroy
   call dev_free(q)
   call report(name//' + '//rk//merge(' (timed)', '        ', timed), err, sp%order())
   endsubroutine with_stepper

   subroutine rotation_residual(u, r, f)
   !< r = (y + f, -x).
   real(R_P), intent(in)    :: u(1:)
   real(R_P), intent(inout) :: r(1:)
   real(R_P), intent(in)    :: f

   !$acc serial DEVICEVAR(u, r)
   r(1) = u(2) + f
   r(2) = -u(1)
   !$acc end serial
   endsubroutine rotation_residual

   function reference_timed() result(u)
   !< The time-dependent problem to T_END by rk4 on the whole system, 20000 steps (error ~1e-17).
   real(R_P)    :: u(2), k1(2), k2(2), k3(2), k4(2), t, dt
   integer(I4P) :: step

   u = U0 ; dt = T_END / 20000._R_P
   do step=1, 20000
      t = (step - 1) * dt
      k1 = full(u, t)
      k2 = full(u + 0.5_R_P*dt*k1, t + 0.5_R_P*dt)
      k3 = full(u + 0.5_R_P*dt*k2, t + 0.5_R_P*dt)
      k4 = full(u + dt*k3, t + dt)
      u = u + dt * (k1 + 2._R_P*k2 + 2._R_P*k3 + k4) / 6._R_P
   enddo
   endfunction reference_timed

   pure function full(u, t) result(f)
   !< The whole time-dependent system: F + G.
   real(R_P), intent(in) :: u(2), t
   real(R_P)             :: f(2)

   f = [u(2) + sin(2._R_P*t) - (1._R_P + t)*u(1), -u(1) - 2._R_P*(1._R_P + t)*u(2)]
   endfunction full

   subroutine stiff(name, first)
   !< u' = a u - k u^2, a = 1: the field step a u by the stepper (rk4), the flow u' = -k u^2 exact.
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: first
   real(R_P),    parameter  :: A = 1._R_P, UI = 1._R_P, KS(6) = [0.01_R_P, 0.1_R_P, 1._R_P, 10._R_P, 100._R_P, &
                                                                   1000._R_P]
   type(additive_splitting) :: sp
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q(:)=>null(), u(:), r(:)
   real(R_P)                :: h(1), err(5), dt, sub, kk, ex
   integer(I4P)             :: j, l, n, step, m, k, ierr

   call sp%initialize(name, first=first)
   call dev_alloc(fptr_dev=q, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   bx = box_interior(lb=[1_I4P], ub=[1_I4P], ngc=0_I4P, components=.false.)
   call st%initialize(SCHEME_RK4, FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
   do j=1, size(KS)
      kk = KS(j)
      do l=1, 5
         n = 8 * 2**(l-1)
         dt = T_END / n
         h = UI
         call dev_memcpy_to_device(dst=q, src=h)
         do step=1, n
            do m=1, sp%substeps()
               sub = sp%substep_fraction(m) * dt
               if (sp%substep_operator(m) == OPERATOR_FIELD) then
                  call st%begin_step(bx, q)
                  do k=1, st%stages()
                     call st%begin_stage(k, bx, sub, q)
                     if (st%stage_state_is_q(k)) then
                        u => q
                     else
                        call st%stage_state_view(u)
                     endif
                     call st%rhs_register_view(r)
                     call growth_residual(A, u, r)
                     call st%end_stage(k, bx, sub, q)
                  enddo
                  call st%finish_step(bx, sub, q)
               else
                  call dev_memcpy_from_device(dst=h, src=q)
                  h = h / (1._R_P + kk * sub * h)
                  call dev_memcpy_to_device(dst=q, src=h)
               endif
            enddo
         enddo
         call dev_memcpy_from_device(dst=h, src=q)
         ex = A * UI * exp(A*T_END) / (A + kk * UI * (exp(A*T_END) - 1._R_P))
         err(l) = abs(h(1) - ex) / ex
      enddo
      call report(name//merge(' (flow first) ', ' (field first)', first == OPERATOR_FLOW)//', k u0 = '// &
                  trim(real_str(kk)), err, sp%order())
   enddo
   call st%destroy
   call dev_free(q)
   endsubroutine stiff

   subroutine growth_residual(a, u, r)
   !< r = a u.
   real(R_P), intent(in)    :: a
   real(R_P), intent(in)    :: u(1:)
   real(R_P), intent(inout) :: r(1:)

   !$acc serial DEVICEVAR(u, r)
   r(1) = a * u(1)
   !$acc end serial
   endsubroutine growth_residual

   function real_str(x) result(s)
   real(R_P), intent(in) :: x
   character(16)         :: s

   write(s, '(ES8.1)') x
   endfunction real_str

   subroutine report(label, err, p)
   !< Print the errors and orders; check the order: the last two halvings whose finer error is above round-off
   !< (1e-12), within 0.15 of p.
   character(*), intent(in) :: label
   real(R_P),    intent(in) :: err(:)
   integer(I4P), intent(in) :: p
   real(R_P)                :: o(size(err)-1)
   integer(I4P)             :: last

   o = observed_orders(real(err, R8P))
   print '(A,5ES10.2,A,4F6.2)', '        '//label//': err', err, '  orders', o
   last = size(o)
   do while (last > 2 .and. err(last+1) < 1.e-12_R_P)
      last = last - 1
   enddo
   call check(all(abs(o(last-1:last) - p) <= 0.15_R_P), label//': order')
   endsubroutine report
endprogram selah_test_additive_splitting
