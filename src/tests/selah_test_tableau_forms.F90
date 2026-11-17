!< SELAH test: the conversions of the coefficient forms to the Butcher tableau, on tableaux known exactly.
program selah_test_tableau_forms
!< SELAH test: the conversions of the coefficient forms to the Butcher tableau, on tableaux known exactly.
!<
!< SSP-22 (Heun) and SSP-33 in two-register and general Shu–Osher form, and Williamson's 3-stage 2N scheme, convert
!< to their textbook Butcher tableaux; the SSP coefficient of the Shu–Osher forms is 1; the order conditions of the
!< exact classic RK4 tableau vanish (a check of the checker).
use penf, only : I4P, R8P
use selah_kind_parameters, only : R_P
use selah_tableaux
use selah_test_tools
implicit none
type(butcher_tableau)    :: bt, ref
type(shu_osher_2r_form)  :: r
type(shu_osher_form)     :: so
type(williamson_2n_form) :: w

! SSP-22 (Heun), two-register: Y2 = Y1 + dt F1; Y3 = 1/2 Y1 + 1/2 Y2 + 1/2 dt F2
r = shu_osher_2r_form(stages=2, a=[1._R_P, 0.5_R_P], b=[0._R_P, 0.5_R_P], c=[1._R_P, 0.5_R_P])
bt = butcher_from_shu_osher_2r(r)
ref = exact_butcher(2, [0._R_P, 0._R_P, &
                        1._R_P, 0._R_P], [0.5_R_P, 0.5_R_P])
call check(same(bt, ref), 'SSP-22 two-register -> Butcher (Heun)')
call check(abs(ssp_coefficient(shu_osher_from_2r(r)) - 1._R_P) <= TOL_COEF, 'SSP-22 SSP coefficient = 1')

! SSP-33, two-register: Y2 = Y1 + dt F1; Y3 = 3/4 Y1 + 1/4 Y2 + 1/4 dt F2; Y4 = 1/3 Y1 + 2/3 Y3 + 2/3 dt F3
r = shu_osher_2r_form(stages=3, a=[1._R_P, 0.75_R_P, 1._R_P/3._R_P], b=[0._R_P, 0.25_R_P, 2._R_P/3._R_P], &
                      c=[1._R_P, 0.25_R_P, 2._R_P/3._R_P])
bt = butcher_from_shu_osher_2r(r)
ref = exact_butcher(3, [0._R_P,  0._R_P,  0._R_P, &
                        1._R_P,  0._R_P,  0._R_P, &
                        0.25_R_P, 0.25_R_P, 0._R_P], [1._R_P/6._R_P, 1._R_P/6._R_P, 2._R_P/3._R_P])
call check(same(bt, ref), 'SSP-33 two-register -> Butcher')
call check(abs(ssp_coefficient(shu_osher_from_2r(r)) - 1._R_P) <= TOL_COEF, 'SSP-33 SSP coefficient = 1')

! SSP-33, general Shu–Osher form written directly
so%stages = 3
allocate(so%alpha(4,3), so%beta(4,3))
so%alpha = 0._R_P ; so%beta = 0._R_P
so%alpha(2,1) = 1._R_P          ; so%beta(2,1) = 1._R_P
so%alpha(3,1) = 0.75_R_P        ; so%alpha(3,2) = 0.25_R_P      ; so%beta(3,2) = 0.25_R_P
so%alpha(4,1) = 1._R_P/3._R_P   ; so%alpha(4,3) = 2._R_P/3._R_P ; so%beta(4,3) = 2._R_P/3._R_P
bt = butcher_from_shu_osher(so, order=3)
call check(same(bt, ref), 'SSP-33 Shu-Osher -> Butcher')
call check(bt%order == 3 .and. abs(bt%ssp_coefficient - 1._R_P) <= TOL_COEF, 'SSP-33 order and SSP coefficient set')
call check(is_explicit(bt), 'SSP-33 explicit')

! Williamson (1980) 3-stage 2N: A = (0, -5/9, -153/128), B = (1/3, 15/16, 8/15)
w = williamson_2n_form(stages=3, a=[0._R_P, -5._R_P/9._R_P, -153._R_P/128._R_P], &
                       b=[1._R_P/3._R_P, 15._R_P/16._R_P, 8._R_P/15._R_P])
bt = butcher_from_williamson_2n(w, order=3)
ref = exact_butcher(3, [0._R_P,          0._R_P,          0._R_P, &
                        1._R_P/3._R_P,   0._R_P,          0._R_P, &
                        -3._R_P/16._R_P, 15._R_P/16._R_P, 0._R_P], [1._R_P/6._R_P, 3._R_P/10._R_P, 8._R_P/15._R_P])
call check(same(bt, ref), 'Williamson 2N -> Butcher')
call check(max_order_residual(bt, 3_I4P) <= TOL_COEF, 'Williamson 2N: order conditions up to 3')
call check(max_order_residual(bt, 4_I4P) > 1.e-3_R8P, 'Williamson 2N: not order 4 (the checker detects it)')

! the checker on the exact classic RK4
ref = exact_butcher(4, [0._R_P,   0._R_P,   0._R_P, 0._R_P, &
                        0.5_R_P,  0._R_P,   0._R_P, 0._R_P, &
                        0._R_P,   0.5_R_P,  0._R_P, 0._R_P, &
                        0._R_P,   0._R_P,   1._R_P, 0._R_P], &
                    [1._R_P/6._R_P, 1._R_P/3._R_P, 1._R_P/3._R_P, 1._R_P/6._R_P])
call check(max_order_residual(ref, 4_I4P) <= TOL_COEF, 'checker: classic RK4 satisfies the 8 conditions')
call check(.not. is_explicit(implicit_midpoint()), 'is_explicit: the implicit midpoint rule is not explicit')

if (failures > 0) error stop 'selah_test_tableau_forms: failed'
print '(A)', 'selah_test_tableau_forms: ok'

contains
   function exact_butcher(s, a_rows, b) result(bt)
   !< A Butcher tableau from A by rows and b; c = A·1.
   integer,   intent(in) :: s         !< Stages.
   real(R_P), intent(in) :: a_rows(:) !< A, row by row [s*s].
   real(R_P), intent(in) :: b(:)      !< Weights [s].
   type(butcher_tableau) :: bt        !< Tableau.

   bt%stages = s
   bt%a = transpose(reshape(a_rows, [s, s]))
   bt%b = b
   bt%c = sum(bt%a, dim=2)
   endfunction exact_butcher

   function implicit_midpoint() result(bt)
   !< The implicit midpoint rule: A = (1/2), b = (1), c = (1/2).
   type(butcher_tableau) :: bt !< Tableau.

   bt%stages = 1
   bt%a = reshape([0.5_R_P], [1, 1])
   bt%b = [1._R_P]
   bt%c = [0.5_R_P]
   endfunction implicit_midpoint

   logical function same(x, y)
   !< Same stages and coefficients within TOL_COEF (relative to 1).
   type(butcher_tableau), intent(in) :: x, y !< Tableaux.

   same = x%stages == y%stages
   if (.not. same) return
   same = maxval(abs(x%a - y%a)) <= TOL_COEF .and. maxval(abs(x%b - y%b)) <= TOL_COEF .and. &
          maxval(abs(x%c - y%c)) <= TOL_COEF
   endfunction same
endprogram selah_test_tableau_forms
