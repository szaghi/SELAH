!< SELAH test: the splitting schemes.
program selah_test_splittings
!< SELAH test: the splitting schemes.
!<
!< For every scheme: each partition's coefficients sum to 1, the composition is a palindrome and alternates the
!< partitions. In R8P builds, the observed order on two separable problems q' = p, p' = -V'(q), with exact sub-flows
!< (drift q += c dt p, kick p -= c dt V'(q)), reaches the declared order, with either partition first: the harmonic
!< oscillator (V = q²/2, exact solution) and the pendulum (V = -cos q, reference by the classic RK4 at 40000 steps in
!< QP).
use penf, only : I4P, R8P
use selah_kind_parameters, only : R_P
use selah_splittings
use selah_test_tools, only : QP, TOL_COEF, check, failures, observed_orders
implicit none
type(splitting_composition) :: sc
character(:), allocatable   :: name ! not associate(name=>trim(...)): gfortran 15.1-16.2 frees that temporary
                                    ! twice in a loop (a reduced reproducer is in the night log of 2026-11-17)
integer(I4P)                :: s, i

do s=1, size(SPLITTING_SCHEMES)
   sc = splitting(SPLITTING_SCHEMES(s))
   name = trim(SPLITTING_SCHEMES(s))
   print '(A)', name//':'
   call check(sc%flows > 0 .and. sc%order >= 2, name//': known, order declared')
   call check(abs(sum(sc%coefficient, mask=sc%partition==1) - 1._R_P) <= TOL_COEF, name//': partition 1 sums to 1')
   call check(abs(sum(sc%coefficient, mask=sc%partition==2) - 1._R_P) <= TOL_COEF, name//': partition 2 sums to 1')
   call check(all(sc%coefficient == sc%coefficient(sc%flows:1:-1)) .and. &
              all(sc%partition == sc%partition(sc%flows:1:-1)), name//': palindrome')
   call check(all([(sc%partition(i) /= sc%partition(i+1), i=1, sc%flows-1)]), name//': partitions alternate')
#if defined _R_P_IS_R4P
   print '(A)', '  skip  '//name//': observed order (R4P coefficients cannot show the order of the exact ones)'
#else
   call check_orders(name, sc%order)
#endif
enddo
sc = splitting('no-such-scheme')
call check(sc%flows == 0, 'unknown name: 0 flows')

if (failures > 0) error stop 'selah_test_splittings: failed'
print '(A)', 'selah_test_splittings: ok'

contains
   subroutine check_orders(name, p)
   !< Observed order with partition 1 = drift and with partition 1 = kick, on both problems.
   character(*), intent(in) :: name !< Scheme name.
   integer(I4P), intent(in) :: p    !< Declared order.
   type(splitting_composition) :: sc
   real(R8P)    :: errors(4), orders(3)
   integer(I4P) :: problem, first, k, n0
   character(1) :: tag(2) = ['1', '2']

   n0 = merge(2_I4P, 8_I4P, p >= 6)   ! order 6 reaches the R8P floor sooner: its finest step (16 per unit time) keeps
                                      ! the error near 1e-13, above the rounding of an R8P-only build (nvfortran)
   do first=1, 2
      sc = splitting(name, first=first)
      do problem=1, 2
         do k=1, 4
            errors(k) = error_at_one(sc, problem, n0 * 2_I4P**(k-1))
         enddo
         orders = observed_orders(errors)
         print '(A,I1,A,I1,A,3F7.3,A,ES9.2)', '        '//name//' first=', first, ' problem ', problem, &
                                             ' observed orders:', orders, '  finest error', errors(4)
         call check(orders(3) >= real(p, R8P) - 0.2_R8P, &
                    name//': observed order, first partition '//tag(first)//', problem '//tag(problem))
      enddo
   enddo
   endsubroutine check_orders

   function error_at_one(sc, problem, n) result(err)
   !< |(q, p)(1) - reference| after n steps of the composition; partition 1 = drift, 2 = kick.
   type(splitting_composition), intent(in) :: sc      !< Composition.
   integer(I4P),                intent(in) :: problem !< 1: harmonic oscillator, 2: pendulum.
   integer(I4P),                intent(in) :: n       !< Steps.
   real(R8P)                               :: err     !< Error at t = 1.
   real(QP)                                :: q, pm, dt, qr, pr
   integer(I4P)                            :: step, m

   dt = 1._QP / real(n, QP)
   q = 1._QP ; pm = 0._QP
   do step=1, n
      do m=1, sc%flows
         if (sc%partition(m) == 1) then
            q = q + real(sc%coefficient(m), QP) * dt * pm
         else
            pm = pm - real(sc%coefficient(m), QP) * dt * dvdq(problem, q)
         endif
      enddo
   enddo
   call reference(problem, qr, pr)
   err = real(sqrt((q - qr)**2 + (pm - pr)**2), R8P)
   endfunction error_at_one

   pure function dvdq(problem, q) result(d)
   !< V'(q): q (harmonic) or sin q (pendulum).
   integer(I4P), intent(in) :: problem !< Problem.
   real(QP),     intent(in) :: q       !< Position.
   real(QP)                 :: d       !< V'(q).

   d = merge(q, sin(q), problem == 1)
   endfunction dvdq

   subroutine reference(problem, qr, pr)
   !< (q, p)(1) from (1, 0): exact for the oscillator, classic RK4 at 40000 steps in QP for the pendulum.
   integer(I4P), intent(in)  :: problem !< Problem.
   real(QP),     intent(out) :: qr, pr  !< Reference state.
   real(QP)                  :: h, k1q, k1p, k2q, k2p, k3q, k3p, k4q, k4p
   integer(I4P)              :: step

   if (problem == 1) then
      qr = cos(1._QP) ; pr = -sin(1._QP)
      return
   endif
   h = 1._QP / 40000._QP
   qr = 1._QP ; pr = 0._QP
   do step=1, 40000
      k1q = pr                  ; k1p = -sin(qr)
      k2q = pr + 0.5_QP*h*k1p   ; k2p = -sin(qr + 0.5_QP*h*k1q)
      k3q = pr + 0.5_QP*h*k2p   ; k3p = -sin(qr + 0.5_QP*h*k2q)
      k4q = pr + h*k3p          ; k4p = -sin(qr + h*k3q)
      qr = qr + h/6._QP * (k1q + 2._QP*k2q + 2._QP*k3q + k4q)
      pr = pr + h/6._QP * (k1p + 2._QP*k2p + 2._QP*k3p + k4p)
   enddo
   endsubroutine reference
endprogram selah_test_splittings
