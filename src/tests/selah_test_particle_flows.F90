!< SELAH test: particle kick and drift sub-flows composed by SELAH's splittings, on the host or the device.

#include "fundal.H"
program selah_test_particle_flows
!< SELAH test: particle kick and drift sub-flows composed by SELAH's splittings, on the host or the device.
!<
!< A particle array (6, np), positions 1:3 and velocities 4:6, stepped by the stepper's splitting form with the
!< velocities as partition 1 and the positions as partition 2: the kick v += c dt a(x) and the drift x += c dt v are
!< exact sub-flows (each partition's residual does not depend on its own variables), and `strang` is kick-drift-kick
!< (KDK) leapfrog, with SELAH's fractions instead of a hard-coded 1/2.
!<
!<+ **Harmonic oscillator, the discrete solution.** x'' = -w^2 x, x(0) = 1, v(0) = 0: KDK's map has trace
!<  2 - (w h)^2 and determinant 1, so x_n = cos(n theta), cos(theta) = 1 - (w h)^2/2, to round-off over 100 periods;
!<  the phase error per period, 2 pi (theta/(w h) - 1) ~ pi (w h)^2/12, against that formula.
!<+ **Kepler orbit** (GM = 1, eccentricity 0.5, one particle in 3-D): the observed order of the position after one
!<  orbit (strang 2, yoshida-4 4, blanes-moan-s6 4, blanes-moan-s10 6); the energy error bounded over 50 orbits (no
!<  drift: its largest value over every step of the last 10 orbits within 1.2 times that over the first 10; sampled
!<  once per orbit it would follow the precessing perihelion instead), every scheme symplectic.
!<+ **Gas and a sink, coupled.** A gas mass m accreted by a sink in a harmonic potential: m' = -kappa x^2 m,
!<  x'' = -(1 + m) x. The field step (the gas, rk4 through the stepper) and the sink's flow (one KDK step, the
!<  splitting stepper, over the substep) composed by selah_additive_splitting: lie order 1, strang order 2, against a
!<  fine rk4 reference of the whole system.
use penf, only : I4P, R8P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures, observed_orders
implicit none
real(R_P),    parameter :: PI = acos(-1._R_P)
integer(I4P), parameter :: KICK_DRIFT(2,2) = reshape([4_I4P, 6_I4P, 1_I4P, 3_I4P], [2, 2]) !< v: partition 1; x: 2.
integer(I4P)            :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'harmonic oscillator, KDK = strang (kick first):'
call oscillator(0.1_R_P)
call oscillator(0.5_R_P)
#if defined _R_P_IS_R4P
print '(A)', '  skip  observed orders (R4P)'
#else
print '(A)', 'Kepler orbit, e = 0.5:'
do s=1, size(SPLITTING_SCHEMES)
   call kepler(trim(SPLITTING_SCHEMES(s)))
enddo
print '(A)', 'gas and a sink, coupled through the additive splitting:'
call gas_and_sink(SCHEME_LIE)
call gas_and_sink(SCHEME_STRANG)
#endif
if (failures > 0) error stop 'selah_test_particle_flows: failed'
print '(A)', 'selah_test_particle_flows: ok'

contains
   subroutine particle_step(st, bx, dt, p, w2)
   !< One step of the particle array p by the splitting stepper: the kick (partition 1) with a = -w2 x, the drift
   !< (partition 2) with v; w2 = 0 for the Kepler force (a = -x/|x|^3).
   type(stepper),      intent(inout) :: st
   type(box),          intent(in)    :: bx
   real(R_P),          intent(in)    :: dt, w2
   real(R_P), pointer, intent(inout) :: p(:,:)
   real(R_P), pointer                :: r(:,:)
   integer(I4P)                      :: k

   call st%begin_step(bx, p)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, p)
      call st%rhs_register_view(r)
      call particle_residual(st%stage_partition(k), w2, p, r)
      call st%end_stage(k, bx, dt, p)
   enddo
   call st%finish_step(bx, dt, p)
   endsubroutine particle_step

   subroutine particle_residual(part, w2, p, r)
   !< Partition 1 (velocities): the acceleration; partition 2 (positions): the velocity.
   integer(I4P), intent(in)    :: part
   real(R_P),    intent(in)    :: w2
   real(R_P),    intent(in)    :: p(1:,1:)
   real(R_P),    intent(inout) :: r(1:,1:)
   integer(I4P)                :: i
   real(R_P)                   :: d3

   !$acc parallel loop independent DEVICEVAR(p, r) private(d3)
   !$omp OMPLOOP DEVICEPTR(p, r) private(d3)
   do i=1, size(p, 2)
      if (part == 1_I4P) then
         if (w2 > 0._R_P) then
            r(4:6,i) = -w2 * p(1:3,i)
         else
            d3 = sqrt(p(1,i)**2 + p(2,i)**2 + p(3,i)**2)**3
            r(4:6,i) = -p(1:3,i) / d3
         endif
      else
         r(1:3,i) = p(4:6,i)
      endif
   enddo
   endsubroutine particle_residual

   subroutine new_particles(p, h, st, bx, name)
   !< Device particles from the host array h, and their splitting stepper.
   real(R_P), pointer, intent(inout) :: p(:,:)
   real(R_P),          intent(in)    :: h(:,:)
   type(stepper),      intent(inout) :: st
   type(box),          intent(out)   :: bx
   character(*),       intent(in)    :: name
   integer(I4P)                      :: ierr, np

   np = size(h, 2, kind=I4P)
   call dev_alloc(fptr_dev=p, ubounds=[6_I4P, np], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_memcpy_to_device(dst=p, src=h)
   bx = box_interior(lb=[1_I4P, 1_I4P], ub=[6_I4P, np], ngc=0_I4P, components=.true.)
   call st%initialize(name, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[6_I4P, np], components=.true., &
                      partitions=KICK_DRIFT)
   endsubroutine new_particles

   subroutine oscillator(wh)
   !< KDK's discrete solution x_n = cos(n theta) over 100 periods, and its phase error per period.
   real(R_P), intent(in) :: wh   !< w dt (w = 1).
   type(stepper)         :: st
   type(box)             :: bx
   real(R_P), pointer    :: p(:,:)=>null()
   real(R_P)             :: h(6,1), theta, worst, phase
   integer(I4P)          :: n, nsteps
   character(8)          :: tag

   write(tag, '(F4.2)') wh
   h = 0._R_P ; h(1,1) = 1._R_P
   call new_particles(p, h, st, bx, SCHEME_STRANG)
   call check(st%stage_partition(1) == 1 .and. st%stage_euler_fraction(1) == 0.5_R_P, &
              'w dt = '//trim(tag)//': strang with the velocities first is a half kick')
   theta = acos(1._R_P - 0.5_R_P * wh**2)
   nsteps = nint(100._R_P * 2._R_P * PI / wh)
   worst = 0._R_P
   do n=1, nsteps
      call particle_step(st, bx, wh, p, 1._R_P)
      call dev_memcpy_from_device(dst=h, src=p)
      worst = max(worst, abs(h(1,1) - cos(n * theta)))
   enddo
   phase = 2._R_P * PI * (theta / wh - 1._R_P)
   print '(A,ES10.2,A,ES11.3,A,ES11.3)', '        w dt = '//trim(tag)//': max |x_n - cos(n theta)|', worst, &
                                         '; phase error per period', phase, ', pi (w dt)^2/12', PI * wh**2 / 12._R_P
   call check(worst <= 1.e4_R_P * epsilon(1._R_P), 'w dt = '//trim(tag)//': x_n = cos(n theta) (KDK''s discrete solution)')
   call check(abs(phase / (PI * wh**2 / 12._R_P) - 1._R_P) < 0.05_R_P, &
              'w dt = '//trim(tag)//': phase error per period pi (w dt)^2/12 to 5 %')
   call st%destroy
   call dev_free(p)
   endsubroutine oscillator

   subroutine kepler(name)
   !< One orbit at four halvings of dt (observed order); 50 orbits for the energy.
   character(*), intent(in) :: name
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: p(:,:)=>null()
   real(R_P)                :: h(6,1), h0(6,1), err(4), e0, de(50), o(3)
   integer(I4P)             :: l, n, step, orb, p_order
   type(splitting_composition) :: sc

   sc = splitting(name)
   p_order = sc%order
   ! perihelion of an orbit of semi-major axis 1, e = 0.5: r = 0.5, v = sqrt((1 + e)/(1 - e)) = sqrt(3); period 2 pi
   h0 = 0._R_P ; h0(1,1) = 0.5_R_P ; h0(5,1) = sqrt(3._R_P)
   do l=1, 4
      n = merge(64, 128, p_order >= 6) * 2**(l-1)
      call new_particles(p, h0, st, bx, name)
      do step=1, n
         call particle_step(st, bx, 2._R_P*PI/n, p, 0._R_P)
      enddo
      call dev_memcpy_from_device(dst=h, src=p)
      err(l) = maxval(abs(h(1:3,1) - h0(1:3,1)))
      call st%destroy
      call dev_free(p)
   enddo
   o = observed_orders(real(err, R8P))
   print '(A,4ES10.2,A,3F6.2)', '        '//name//': position error after one orbit', err, '  orders', o
   call check(all(abs(o(2:) - p_order) <= 0.2_R_P), name//': order on the Kepler orbit')
   ! energy over 50 orbits, 256 steps per orbit
   call new_particles(p, h0, st, bx, name)
   e0 = energy(h0(:,1))
   do orb=1, 50
      de(orb) = 0._R_P
      do step=1, 256
         call particle_step(st, bx, 2._R_P*PI/256._R_P, p, 0._R_P)
         call dev_memcpy_from_device(dst=h, src=p)
         de(orb) = max(de(orb), abs(energy(h(:,1)) - e0))
      enddo
   enddo
   print '(A,ES10.2,A,ES10.2)', '        '//name//': energy error, largest over orbits 1-10', maxval(de(1:10)), &
                                ', over orbits 41-50', maxval(de(41:50))
   call check(maxval(de(41:50)) <= 1.2_R_P * maxval(de(1:10)) + 1.e-13_R_P, name//': energy error bounded, no drift')
   call st%destroy
   call dev_free(p)
   endsubroutine kepler

   pure function energy(y) result(e)
   !< Kepler energy, GM = 1.
   real(R_P), intent(in) :: y(6)
   real(R_P)             :: e

   e = 0.5_R_P * sum(y(4:6)**2) - 1._R_P / sqrt(sum(y(1:3)**2))
   endfunction energy

   subroutine gas_and_sink(name)
   !< The gas (rk4, the field step) and the sink (one KDK step over the substep, the flow): observed order.
   character(*), intent(in)  :: name
   real(R_P),    parameter   :: KAPPA = 2._R_P, TE = 2._R_P
   type(additive_splitting)  :: sp
   type(stepper)             :: gas, sink
   type(box)                 :: bg, bs
   real(R_P), pointer        :: m(:)=>null(), p(:,:)=>null(), u(:), r(:)
   real(R_P)                 :: hm(1), hp(6,1), err(4), dt, sub, ref(2)
   integer(I4P)              :: l, n, step, j, k, ierr

   call sp%initialize(name)
   ref = reference()
   call dev_alloc(fptr_dev=m, ubounds=[1_I4P], lbounds=[1_I4P], ierr=ierr)
   bg = box_interior(lb=[1_I4P], ub=[1_I4P], ngc=0_I4P, components=.false.)
   call gas%initialize(SCHEME_RK4, FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
   do l=1, 4
      n = 16 * 2**(l-1)
      dt = TE / n
      hm = 1._R_P
      call dev_memcpy_to_device(dst=m, src=hm)
      hp = 0._R_P ; hp(1,1) = 1._R_P
      call new_particles(p, hp, sink, bs, SCHEME_STRANG)
      do step=1, n
         do j=1, sp%substeps()
            sub = sp%substep_fraction(j) * dt
            if (sp%substep_operator(j) == OPERATOR_FIELD) then ! the gas, the sink frozen
               call dev_memcpy_from_device(dst=hp, src=p)
               call gas%begin_step(bg, m)
               do k=1, gas%stages()
                  call gas%begin_stage(k, bg, sub, m)
                  if (gas%stage_state_is_q(k)) then
                     u => m
                  else
                     call gas%stage_state_view(u)
                  endif
                  call gas%rhs_register_view(r)
                  call accretion(KAPPA * hp(1,1)**2, u, r)
                  call gas%end_stage(k, bg, sub, m)
               enddo
               call gas%finish_step(bg, sub, m)
            else ! the sink, the gas frozen: w^2 = 1 + m
               call dev_memcpy_from_device(dst=hm, src=m)
               call particle_step(sink, bs, sub, p, 1._R_P + hm(1))
            endif
         enddo
      enddo
      call dev_memcpy_from_device(dst=hm, src=m)
      call dev_memcpy_from_device(dst=hp, src=p)
      err(l) = max(abs(hm(1) - ref(1)), abs(hp(1,1) - ref(2)))
      call sink%destroy
      call dev_free(p)
   enddo
   call gas%destroy
   call dev_free(m)
   print '(A,4ES10.2,A,3F6.2)', '        '//name//': error', err, '  orders', observed_orders(real(err, R8P))
   call check(all(abs(observed_orders(real(err, R8P)) - sp%order()) <= 0.15_R_P), name//': gas and sink, order')
   endsubroutine gas_and_sink

   subroutine accretion(rate, u, r)
   !< r = -rate u.
   real(R_P), intent(in)    :: rate
   real(R_P), intent(in)    :: u(1:)
   real(R_P), intent(inout) :: r(1:)

   !$acc serial DEVICEVAR(u, r)
   r(1) = -rate * u(1)
   !$acc end serial
   endsubroutine accretion

   function reference() result(y)
   !< (m, x) at t = 2 by rk4 on the whole system, 40000 steps.
   real(R_P)    :: y(2), z(3), k1(3), k2(3), k3(3), k4(3), dt
   integer(I4P) :: step

   z = [1._R_P, 1._R_P, 0._R_P] ; dt = 2._R_P / 40000._R_P
   do step=1, 40000
      k1 = f(z) ; k2 = f(z + 0.5_R_P*dt*k1) ; k3 = f(z + 0.5_R_P*dt*k2) ; k4 = f(z + dt*k3)
      z = z + dt * (k1 + 2._R_P*k2 + 2._R_P*k3 + k4) / 6._R_P
   enddo
   y = z(1:2)
   endfunction reference

   pure function f(z) result(d)
   !< (m, x, v)' = (-kappa x^2 m, v, -(1 + m) x), kappa = 2.
   real(R_P), intent(in) :: z(3)
   real(R_P)             :: d(3)

   d = [-2._R_P * z(2)**2 * z(1), z(3), -(1._R_P + z(1)) * z(2)]
   endfunction f
endprogram selah_test_particle_flows
