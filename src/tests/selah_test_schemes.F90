!< SELAH test: the Runge–Kutta catalogue, scheme by scheme.
program selah_test_schemes
!< SELAH test: the Runge–Kutta catalogue, scheme by scheme.
!<
!< For every scheme: the order conditions up to the declared order vanish (<= TOL_COEF, evaluated in QP) and some
!< condition of the next order does not (the order is exact); c = A·1; the tableau is explicit; the SSP coefficient,
!< computed from the Shu–Osher form, matches the literature; every other form converts to the same tableau. Then the
!< published facts a typo would break: the abscissas of Carpenter & Kennedy, the weights 1/10 of SSPRK(10,4). In R8P
!< builds, the observed order on y' = -y and y' = -y² (4 halvings from 8 steps) reaches the declared order.
use penf, only : I4P, R8P
use selah_kind_parameters, only : R_P
use selah_tableaux
use selah_schemes
use selah_test_tools
implicit none
type(butcher_tableau) :: bt
integer(I4P)          :: s

do s=1, size(RK_SCHEMES)
   bt = butcher(RK_SCHEMES(s))
   print '(A)', trim(RK_SCHEMES(s))//':'
   call check(bt%stages > 0 .and. bt%order >= 1, trim(RK_SCHEMES(s))//': known, order declared')
   call check_scheme(trim(RK_SCHEMES(s)), bt)
enddo

print '(A)', 'published values:'
bt = butcher(SCHEME_SSP_54)
call check(abs(bt%ssp_coefficient - 1.508_R_P) < 5.e-4_R_P, 'ssp-54: SSP coefficient 1.508 (Spiteri & Ruuth)')
print '(A,F12.9)', '        ssp-54 SSP coefficient = ', bt%ssp_coefficient
bt = butcher(SCHEME_SSP_104)
call check(abs(bt%ssp_coefficient - 6._R_P) <= TOL_COEF * 6._R_P, 'ssp-104: SSP coefficient 6 (Ketcheson)')
call check(maxval(abs(bt%b - 0.1_R_P)) <= TOL_COEF, 'ssp-104: every weight 1/10 (Ketcheson)')
call check(maxval(abs(bt%c - [0._R_P, 1._R_P/6._R_P, 1._R_P/3._R_P, 0.5_R_P, 2._R_P/3._R_P, 1._R_P/3._R_P, 0.5_R_P, &
                              2._R_P/3._R_P, 5._R_P/6._R_P, 1._R_P])) <= TOL_COEF, 'ssp-104: abscissas (Ketcheson)')
bt = butcher(SCHEME_LSRK_54)
call check(maxval(abs(bt%c - [0._R_P, 1432997174477._R_P/9575080441755._R_P, 2526269341429._R_P/6820363962896._R_P, &
                              2006345519317._R_P/3224310063776._R_P, 2802321613138._R_P/2924317926251._R_P])) &
           <= TOL_COEF, 'lsrk-54: abscissas (Carpenter & Kennedy)')
bt = butcher('no-such-scheme')
call check(bt%stages == 0, 'unknown name: 0 stages')

if (failures > 0) error stop 'selah_test_schemes: failed'
print '(A)', 'selah_test_schemes: ok'

contains
   subroutine check_scheme(name, bt)
   !< The checks every scheme passes.
   character(*),          intent(in) :: name !< Scheme name.
   type(butcher_tableau), intent(in) :: bt   !< Its tableau.
   type(shu_osher_form)              :: so                           !< Shu-Osher form.
   real(R8P)                         :: res(8), errors(5), orders(4) !< Residuals, errors, observed orders.
   integer(I4P)                      :: p, problem, k                !< Order, problem, counter.

   p = bt%order
   res = order_residuals(bt)
   call check(maxval(res(1:conditions_number(p))) <= TOL_COEF, name//': order conditions up to the declared order')
   print '(A,ES9.2)', '        '//name//' largest residual up to the declared order:', maxval(res(1:conditions_number(p)))
   if (p < 4) call check(maxval(res(conditions_number(p)+1:conditions_number(p+1))) > 1.e-6_R8P, &
                         name//': some condition of the next order fails (the order is exact)')
   call check(maxval(abs(bt%c - sum(bt%a, dim=2))) <= TOL_COEF, name//': c = A.1')
   call check(is_explicit(bt), name//': explicit')
   select case(name)
   case(SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33)
      call check(abs(bt%ssp_coefficient - 1._R_P) <= TOL_COEF, name//': SSP coefficient 1')
      call check(same(bt, butcher_from_shu_osher(shu_osher(name))), name//': Shu-Osher form -> same tableau')
      call check(same(bt, butcher_from_shu_osher_2r(shu_osher_2r(name))), name//': two-register form -> same tableau')
   case(SCHEME_RK4, SCHEME_LSRK_33, SCHEME_LSRK_54)
      call check(bt%ssp_coefficient == 0._R_P, name//': not SSP (coefficient 0)')
      so = shu_osher(name)
      call check(so%stages == 0, name//': no Shu-Osher form in the catalogue')
   endselect
#if defined _R_P_IS_R4P
   print '(A)', '  skip  '//name//': observed order (R4P coefficients cannot show the order of the exact ones)'
#else
   do problem=1, 2
      do k=1, 5
         errors(k) = integrate_scalar(bt, problem, 8_I4P * 2_I4P**(k-1))
      enddo
      orders = observed_orders(errors)
      print '(A,I1,A,4F7.3)', '        '//name//' problem ', problem, ' observed orders:', orders
      call check(orders(4) >= real(p, R8P) - 0.2_R8P, name//': observed order on problem '//merge('1','2',problem==1))
   enddo
#endif
   endsubroutine check_scheme

   logical function same(x, y)
   !< Same stages and coefficients within TOL_COEF.
   type(butcher_tableau), intent(in) :: x, y !< Tableaux.

   same = x%stages == y%stages
   if (.not. same) return
   same = maxval(abs(x%a - y%a)) <= TOL_COEF .and. maxval(abs(x%b - y%b)) <= TOL_COEF .and. &
          maxval(abs(x%c - y%c)) <= TOL_COEF
   endfunction same
endprogram selah_test_schemes
