!< SELAH test: companion states stepped with the field stepper, on the host or the device.

#include "fundal.H"
program selah_test_companion
!< SELAH test: companion states stepped with the field stepper, on the host or the device.
!<
!< A field y (one cell) and NP particles z = (z1, z2, c) coupled both ways inside each stage, as PRISM's fields and PIC
!< particles are: y' = -y + mean_i(c_i z1_i^2), z1_i' = z2_i, z2_i' = -c_i (1 + y) z1_i (c_i constant: the charge, a
!< component that does not evolve). The field has its stepper, the particles a companion (`initialize_companion`),
!< driven through the same calls in the same stage window.
!<
!<+ **Orders.** Every Runge–Kutta scheme in the Butcher form, and the two-register, 2N and Ketcheson forms: the coupled
!<  system reaches the scheme's order (reference: rk4 on the whole system, 20 000 steps).
!<+ **Bitwise one state.** The same steps with field and particles as one state of one stepper: bitwise equal, every
!<  form (the same kernels, element by element).
!<+ **Active extent and capacity.** The companion allocated for 8 particles, the box covering n = 8, 5, 6 of them in
!<  turn (particles leaving and arriving between steps): the particles past n and the constant component c keep their
!<  values bit for bit, and the active ones are bitwise those of a companion allocated for exactly n.
!<+ **Both layouts.** The particles as (3, np) (components first, the host's layout) and as (np, 3) (a scalar rank-2
!<  field whose second index is the component, the device's layout): bitwise equal.
!<+ **Rebuild.** The companion re-initialised between two steps with a new capacity (a regrid): bitwise the run
!<  without the rebuild.
use penf, only : I4P, R8P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures, observed_orders
implicit none
integer(I4P), parameter :: NP = 4                                   !< Particles.
real(R_P),    parameter :: T_END = 1._R_P                           !< Final time.
character(7), parameter :: OTHER(4) = [character(7) :: SCHEME_SSP_33, SCHEME_LSRK_54, SCHEME_SSP_104, SCHEME_SSP_54]
character(14),parameter :: OTHER_FORMS(4) = [character(14) :: FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, &
                                             FORM_KETCHESON_2R, FORM_BUTCHER]
integer(I4P)            :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'field and companion against one state (bitwise), both layouts:'
do s=1, size(RK_SCHEMES)
   call one_state(trim(RK_SCHEMES(s)), FORM_BUTCHER, .false.)
enddo
call one_state(SCHEME_SSP_54, FORM_BUTCHER, .true.)
do s=1, 3
   call one_state(trim(OTHER(s)), trim(OTHER_FORMS(s)), .false.)
enddo
print '(A)', 'active extent, capacity and rebuild:'
call extent(SCHEME_SSP_54, FORM_BUTCHER)
call extent(SCHEME_SSP_33, FORM_SHU_OSHER_2R)
call extent(SCHEME_LSRK_54, FORM_WILLIAMSON_2N)
#if defined _R_P_IS_R4P
print '(A)', '  skip  observed orders (R4P)'
#else
print '(A)', 'observed order of the coupled system:'
do s=1, size(RK_SCHEMES)
   call order(trim(RK_SCHEMES(s)), FORM_BUTCHER)
enddo
do s=1, 3
   call order(trim(OTHER(s)), trim(OTHER_FORMS(s)))
enddo
#endif
if (failures > 0) error stop 'selah_test_companion: failed'
print '(A)', 'selah_test_companion: ok'

contains
   pure function initial_particles(n) result(h)
   !< z1, z2, c of n particles (host layout).
   integer(I4P), intent(in) :: n
   real(R_P)                :: h(3,n)
   integer(I4P)             :: i

   do i=1, n
      h(:,i) = [0.5_R_P + 0.1_R_P*i, 0.2_R_P - 0.05_R_P*i, 1._R_P + 0.25_R_P*i]
   enddo
   endfunction initial_particles

   subroutine field_residual(n, y, w, z, r)
   !< r = -y + mean_i(c_i z1_i^2), particles in the host layout (3, n): z1 of the stage state w, c of the solution z (the
   !< stage state holds the evolving components only).
   integer(I4P), intent(in)    :: n
   real(R_P),    intent(in)    :: y(1:), w(1:,1:), z(1:,1:)
   real(R_P),    intent(inout) :: r(1:)
   real(R_P)                   :: acc
   integer(I4P)                :: i

   acc = 0._R_P
   !$acc parallel loop DEVICEVAR(w, z) reduction(+:acc)
   !$omp OMPLOOP DEVICEPTR(w, z) reduction(+:acc)
   do i=1, n
      acc = acc + z(3,i) * w(1,i)**2
   enddo
   !$acc serial DEVICEVAR(y, r)
   r(1) = -y(1) + acc / n
   !$acc end serial
   endsubroutine field_residual

   subroutine particle_residual(n, y, w, z, r)
   !< r = (z2, -c (1 + y) z1) of each particle, host layout (z1, z2 of the stage state w, c of the solution z).
   integer(I4P), intent(in)    :: n
   real(R_P),    intent(in)    :: y(1:), w(1:,1:), z(1:,1:)
   real(R_P),    intent(inout) :: r(1:,1:)
   integer(I4P)                :: i

   !$acc parallel loop independent DEVICEVAR(y, w, z, r)
   !$omp OMPLOOP DEVICEPTR(y, w, z, r)
   do i=1, n
      r(1,i) = w(2,i)
      r(2,i) = -z(3,i) * (1._R_P + y(1)) * w(1,i)
   enddo
   endsubroutine particle_residual

   subroutine particle_residual_t(n, y, w, z, r)
   !< particle_residual in the device layout (n, 3).
   integer(I4P), intent(in)    :: n
   real(R_P),    intent(in)    :: y(1:), w(1:,1:), z(1:,1:)
   real(R_P),    intent(inout) :: r(1:,1:)
   integer(I4P)                :: i

   !$acc parallel loop independent DEVICEVAR(y, w, z, r)
   !$omp OMPLOOP DEVICEPTR(y, w, z, r)
   do i=1, n
      r(i,1) = w(i,2)
      r(i,2) = -z(i,3) * (1._R_P + y(1)) * w(i,1)
   enddo
   endsubroutine particle_residual_t

   subroutine field_residual_t(n, y, w, z, r)
   !< field_residual with the particles in the device layout (n, 3).
   integer(I4P), intent(in)    :: n
   real(R_P),    intent(in)    :: y(1:), w(1:,1:), z(1:,1:)
   real(R_P),    intent(inout) :: r(1:)
   real(R_P)                   :: acc
   integer(I4P)                :: i

   acc = 0._R_P
   !$acc parallel loop DEVICEVAR(w, z) reduction(+:acc)
   !$omp OMPLOOP DEVICEPTR(w, z) reduction(+:acc)
   do i=1, n
      acc = acc + z(i,3) * w(i,1)**2
   enddo
   !$acc serial DEVICEVAR(y, r)
   r(1) = -y(1) + acc / n
   !$acc end serial
   endsubroutine field_residual_t

   subroutine coupled_step(st, sp, dt, y, z, n, layout_t)
   !< One step of field and companion: the same calls in the same stage window.
   type(stepper),      intent(inout) :: st, sp
   real(R_P),          intent(in)    :: dt
   real(R_P), pointer, intent(inout) :: y(:), z(:,:)
   integer(I4P),       intent(in)    :: n        !< Active particles.
   logical,            intent(in)    :: layout_t !< Particles as (np, 3).
   type(box)                         :: by, bz
   real(R_P), pointer                :: u(:), w(:,:), ry(:), rz(:,:)
   integer(I4P)                      :: k

   by = box_interior(lb=[1_I4P], ub=[1_I4P], ngc=0_I4P, components=.false.)
   if (layout_t) then
      bz = box_range(lb=[1_I4P, 1_I4P], lo=[1_I4P, 1_I4P], hi=[n, 2_I4P], components=.false.)
   else
      bz = box_range(lb=[1_I4P, 1_I4P], lo=[1_I4P, 1_I4P], hi=[2_I4P, n], components=.true.)
   endif
   call st%begin_step(by, y) ; call sp%begin_step(bz, z)
   do k=1, st%stages()
      call st%begin_stage(k, by, dt, y) ; call sp%begin_stage(k, bz, dt, z)
      if (st%stage_state_is_q(k)) then
         u => y
      else
         call st%stage_state_view(u)
      endif
      if (sp%stage_state_is_q(k)) then
         w => z
      else
         call sp%stage_state_view(w)
      endif
      call st%rhs_register_view(ry) ; call sp%rhs_register_view(rz)
      if (layout_t) then
         call field_residual_t(n, u, w, z, ry) ; call particle_residual_t(n, u, w, z, rz)
      else
         call field_residual(n, u, w, z, ry) ; call particle_residual(n, u, w, z, rz)
      endif
      call st%end_stage(k, by, dt, y) ; call sp%end_stage(k, bz, dt, z)
   enddo
   call st%finish_step(by, dt, y) ; call sp%finish_step(bz, dt, z)
   endsubroutine coupled_step

   subroutine one_state_step(sw, dt, w)
   !< One step of the same system as one state w = (y, z1_1..z1_NP, z2_1..z2_NP, c_1..c_NP), c not stepped.
   type(stepper),      intent(inout) :: sw
   real(R_P),          intent(in)    :: dt
   real(R_P), pointer, intent(inout) :: w(:)
   type(box)                         :: bw
   real(R_P), pointer                :: u(:), r(:)
   integer(I4P)                      :: k

   bw = box_range(lb=[1_I4P], lo=[1_I4P], hi=[1_I4P + 2_I4P*NP], components=.false.)
   call sw%begin_step(bw, w)
   do k=1, sw%stages()
      call sw%begin_stage(k, bw, dt, w)
      if (sw%stage_state_is_q(k)) then
         u => w
      else
         call sw%stage_state_view(u)
      endif
      call sw%rhs_register_view(r)
      call one_state_residual(u, w, r)
      call sw%end_stage(k, bw, dt, w)
   enddo
   call sw%finish_step(bw, dt, w)
   endsubroutine one_state_step

   subroutine one_state_residual(u, w, r)
   !< The residual of the one-state system (the same expressions as field_residual and particle_residual); c from w.
   real(R_P), intent(in)    :: u(1:), w(1:)
   real(R_P), intent(inout) :: r(1:)
   real(R_P)                :: acc
   integer(I4P)             :: i, n

   ! the particle count at run time, as field_residual's n: with the constant NP nvfortran (R4P, -O2 and above) sums
   ! the 4 terms of acc in another order than the runtime loop, and the one state is no longer bitwise the field
   n = (size(w) - 1) / 3

   acc = 0._R_P
   !$acc parallel loop DEVICEVAR(u, w) reduction(+:acc)
   !$omp OMPLOOP DEVICEPTR(u, w) reduction(+:acc)
   do i=1, n
      acc = acc + w(1+2*n+i) * u(1+i)**2
   enddo
   !$acc serial DEVICEVAR(u, r)
   r(1) = -u(1) + acc / n
   !$acc end serial
   !$acc parallel loop independent DEVICEVAR(u, w, r)
   !$omp OMPLOOP DEVICEPTR(u, w, r)
   do i=1, n
      r(1+i) = u(1+n+i)
      r(1+n+i) = -w(1+2*n+i) * (1._R_P + u(1)) * u(1+i)
   enddo
   endsubroutine one_state_residual

   subroutine new_state(y, z, zt, w)
   !< The initial field, particles (both layouts) and one state, on the device.
   real(R_P), pointer, intent(inout) :: y(:), z(:,:), zt(:,:), w(:)
   real(R_P)                         :: h(3,NP)
   integer(I4P)                      :: ierr

   h = initial_particles(NP)
   call dev_alloc(fptr_dev=y, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=z, ubounds=[3_I4P, NP], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=zt, ubounds=[NP, 3_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=w, ubounds=[1_I4P + 3_I4P*NP], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=y, src=[0.3_R_P])
   call dev_memcpy_to_device(dst=z, src=h)
   call dev_memcpy_to_device(dst=zt, src=transpose(h))
   call dev_memcpy_to_device(dst=w, src=[0.3_R_P, h(1,:), h(2,:), h(3,:)])
   endsubroutine new_state

   subroutine one_state(name, form, on_q)
   !< Field + companion (both layouts) against one state, 8 steps: bitwise.
   character(*), intent(in) :: name, form
   logical,      intent(in) :: on_q
   type(stepper)            :: st, sp, st2, spt, sw
   real(R_P), pointer       :: y(:)=>null(), y2(:)=>null(), z(:,:)=>null(), zt(:,:)=>null(), w(:)=>null()
   real(R_P)                :: hy(1), hz(3,NP), hzt(NP,3), hw(1+3*NP), dt
   integer(I4P)             :: step, ierr

   call new_state(y, z, zt, w)
   call dev_alloc(fptr_dev=y2, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=y2, src=[0.3_R_P])
   call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false., stages_on_q=on_q)
   call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[3_I4P, NP], components=.true.)
   call st2%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false., stages_on_q=on_q)
   call spt%initialize_companion(st2, lb=[1_I4P, 1_I4P], ub=[NP, 3_I4P], components=.false.)
   call sw%initialize(name, form, lb=[1_I4P], ub=[1_I4P + 3_I4P*NP], components=.false., stages_on_q=on_q)
   dt = 0.125_R_P
   do step=1, 8
      call coupled_step(st, sp, dt, y, z, NP, .false.)
      call coupled_step(st2, spt, dt, y2, zt, NP, .true.)
      call one_state_step(sw, dt, w)
   enddo
   call dev_memcpy_from_device(dst=hy, src=y) ; call dev_memcpy_from_device(dst=hz, src=z)
   call dev_memcpy_from_device(dst=hzt, src=zt) ; call dev_memcpy_from_device(dst=hw, src=w)
   call check(hy(1) == hw(1) .and. all(hz(1,:) == hw(2:1+NP)) .and. all(hz(2,:) == hw(2+NP:1+2*NP)) .and. &
              all(hz(3,:) == hw(2+2*NP:)), name//' '//form//merge(' (stages on q)', '              ', on_q)// &
              ': field and companion bitwise one state')
   call check(all(hzt == transpose(hz)), name//' '//form//merge(' (stages on q)', '              ', on_q)// &
              ': the (np, 3) layout bitwise the (3, np) one')
   call check(sp%stages() == st%stages() .and. (sp%stage_state_is_q(1) .eqv. st%stage_state_is_q(1)), &
              name//' '//form//': the companion has the field''s stages')
   call st%destroy ; call sp%destroy ; call st2%destroy ; call spt%destroy ; call sw%destroy
   call dev_free(y) ; call dev_free(y2) ; call dev_free(z) ; call dev_free(zt) ; call dev_free(w)
   endsubroutine one_state

   subroutine extent(name, form)
   !< Capacity 8, active 8, 5, 6 in turn; the rebuild at a new capacity; against a companion of exactly n.
   character(*), intent(in) :: name, form
   integer(I4P), parameter  :: CAP = 8, NS(3) = [8, 5, 6]
   type(stepper)            :: st, sp, st1, sp1
   real(R_P), pointer       :: y(:)=>null(), z(:,:)=>null(), y1(:)=>null(), z1(:,:)=>null()
   real(R_P)                :: hz(3,CAP), h1(3,CAP), hy(1), hy1(1), dt
   integer(I4P)             :: step, n, ierr
   logical                  :: same, idle

   h1 = initial_particles(CAP)
   call dev_alloc(fptr_dev=y, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=z, ubounds=[3_I4P, CAP], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=y, src=[0.3_R_P]) ; call dev_memcpy_to_device(dst=z, src=h1)
   call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
   call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[3_I4P, CAP], components=.true.)
   call st1%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
   dt = 0.125_R_P
   same = .true. ; idle = .true.
   do step=1, 3
      n = NS(step)
      ! the reference: field and a companion of exactly n particles, from the same state
      call dev_memcpy_from_device(dst=hz, src=z) ; call dev_memcpy_from_device(dst=hy, src=y)
      call dev_alloc(fptr_dev=y1, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
      call dev_alloc(fptr_dev=z1, ubounds=[3_I4P, n], lbounds=[1_I4P, 1_I4P], ierr=ierr)
      call dev_memcpy_to_device(dst=y1, src=hy) ; call dev_memcpy_to_device(dst=z1, src=hz(:,1:n))
      call sp1%initialize_companion(st1, lb=[1_I4P, 1_I4P], ub=[3_I4P, n], components=.true.)
      call coupled_step(st1, sp1, dt, y1, z1, n, .false.)
      if (step == 2) then ! a regrid: the companion rebuilt at a new capacity between two steps
         call sp%destroy
         call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[3_I4P, CAP], components=.true.)
      endif
      call coupled_step(st, sp, dt, y, z, n, .false.)
      call dev_memcpy_from_device(dst=h1, src=z) ; call dev_memcpy_from_device(dst=hy1, src=y)
      block
         real(R_P) :: r1(3,n), ry(1)
         call dev_memcpy_from_device(dst=r1, src=z1) ; call dev_memcpy_from_device(dst=ry, src=y1)
         same = same .and. all(h1(:,1:n) == r1) .and. hy1(1) == ry(1)
      endblock
      idle = idle .and. all(h1(:,n+1:) == hz(:,n+1:)) .and. all(h1(3,:) == hz(3,:))
      call sp1%destroy ; call dev_free(y1) ; call dev_free(z1)
   enddo
   call check(same, name//' '//form//': n of 8 particles active, bitwise a companion of exactly n (and the rebuild)')
   call check(idle, name//' '//form//': the particles past n and the constant component untouched')
   call st%destroy ; call sp%destroy ; call st1%destroy
   call dev_free(y) ; call dev_free(z)
   endsubroutine extent

   subroutine order(name, form)
   !< The observed order of the coupled system at four halvings of dt.
   character(*), intent(in) :: name, form
   type(stepper)            :: st, sp
   type(butcher_tableau)    :: bt
   real(R_P), pointer       :: y(:)=>null(), z(:,:)=>null(), zt(:,:)=>null(), w(:)=>null()
   real(R_P)                :: hy(1), hz(3,NP), err(4), ref(1+2*NP), o(3)
   integer(I4P)             :: l, n, step

   bt = butcher(name)
   ref = reference()
   do l=1, 4
      n = merge(8, 16, bt%order >= 4) * 2**(l-1)
      call new_state(y, z, zt, w)
      call st%initialize(name, form, lb=[1_I4P], ub=[1_I4P], components=.false.)
      call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[3_I4P, NP], components=.true.)
      do step=1, n
         call coupled_step(st, sp, T_END/n, y, z, NP, .false.)
      enddo
      call dev_memcpy_from_device(dst=hy, src=y) ; call dev_memcpy_from_device(dst=hz, src=z)
      err(l) = maxval(abs([hy(1), hz(1,:), hz(2,:)] - ref))
      call st%destroy ; call sp%destroy
      call dev_free(y) ; call dev_free(z) ; call dev_free(zt) ; call dev_free(w)
   enddo
   o = observed_orders(real(err, R8P))
   print '(A,4ES10.2,A,3F6.2)', '        '//name//' '//form//': err', err, '  orders', o
   call check(all(abs(o(2:) - bt%order) <= 0.2_R_P), name//' '//form//': the scheme''s order, coupled')
   endsubroutine order

   function reference() result(v)
   !< (y, z1, z2) at T_END by rk4 on the whole system, 20000 steps.
   real(R_P)    :: v(1+2*NP), k1(1+2*NP), k2(1+2*NP), k3(1+2*NP), k4(1+2*NP), c(NP), h(3,NP), dt
   integer(I4P) :: step

   h = initial_particles(NP)
   v = [0.3_R_P, h(1,:), h(2,:)] ; c = h(3,:)
   dt = T_END / 20000._R_P
   do step=1, 20000
      k1 = f(v, c) ; k2 = f(v + 0.5_R_P*dt*k1, c) ; k3 = f(v + 0.5_R_P*dt*k2, c) ; k4 = f(v + dt*k3, c)
      v = v + dt * (k1 + 2._R_P*k2 + 2._R_P*k3 + k4) / 6._R_P
   enddo
   endfunction reference

   pure function f(v, c) result(d)
   !< The whole system's right-hand side.
   real(R_P), intent(in) :: v(1+2*NP), c(NP)
   real(R_P)             :: d(1+2*NP)

   d(1) = -v(1) + sum(c * v(2:1+NP)**2) / NP
   d(2:1+NP) = v(2+NP:)
   d(2+NP:) = -c * (1._R_P + v(1)) * v(2:1+NP)
   endfunction f
endprogram selah_test_companion
