!< SELAH tests: shared tools (a quad kind, the order conditions, a scalar integrator for observed orders).
module selah_test_tools
!< SELAH tests: shared tools (a quad kind, the order conditions, a scalar integrator for observed orders).
!<
!< The order conditions are evaluated in QP, quadruple precision where the compiler has it (gfortran, ifx) and R8P
!< where it does not (nvfortran 26.9): the residual is then the tableau's, not the arithmetic's. The tolerances are
!< tied to the epsilon of R_P, the kind the coefficients are stored in.
use penf, only : I4P, R8P
use selah_kind_parameters, only : R_P
use selah_tableaux, only : butcher_tableau
implicit none
private
public :: QP
public :: TOL_COEF
public :: order_residuals
public :: max_order_residual
public :: conditions_number
public :: integrate_scalar
public :: observed_orders
public :: check
public :: failures

integer, parameter :: QP_CANDIDATE = selected_real_kind(33, 4931)                    !< Quad kind, negative if none.
integer, parameter :: QP = merge(QP_CANDIDATE, R8P, QP_CANDIDATE > 0)                 !< Quad kind, else R8P.
real(R8P), parameter :: TOL_COEF = 64._R8P * real(epsilon(1._R_P), R8P)              !< Coefficient tolerance.
integer(I4P), parameter :: CONDITIONS_UP_TO(4) = [1_I4P, 2_I4P, 4_I4P, 8_I4P]          !< Conditions up to order p.
integer(I4P), save :: failures = 0_I4P                                               !< Failed checks.

contains
   pure function conditions_number(p) result(n)
   !< Number of the order conditions up to order p (1 <= p <= 4): 1, 2, 4, 8.
   integer(I4P), intent(in) :: p !< Order.
   integer(I4P)             :: n !< Conditions number.

   n = CONDITIONS_UP_TO(p)
   endfunction conditions_number

   pure function order_residuals(bt) result(res)
   !< The eight order conditions up to order 4 (rooted trees of 1..4 nodes), evaluated in QP on the stored coefficients.
   !<
   !< 1: b·1 = 1; 2: b·c = 1/2; 3: b·c² = 1/3; 4: b·Ac = 1/6; 5: b·c³ = 1/4; 6: b·(c∘Ac) = 1/8; 7: b·Ac² = 1/12;
   !< 8: b·A²c = 1/24. The stored c is used: c = A·1 is checked apart.
   type(butcher_tableau), intent(in) :: bt     !< Butcher tableau.
   real(R8P)                         :: res(8) !< |residual| of each condition.
   real(QP), allocatable             :: a(:,:), b(:), c(:), ac(:), ac2(:), aac(:) !< QP copies and products.

   a = real(bt%a, QP)
   b = real(bt%b, QP)
   c = real(bt%c, QP)
   ac  = matmul(a, c)
   ac2 = matmul(a, c*c)
   aac = matmul(a, ac)
   res(1) = real(abs(sum(b)           - 1._QP        ), R8P)
   res(2) = real(abs(sum(b*c)         - 1._QP/2._QP  ), R8P)
   res(3) = real(abs(sum(b*c*c)       - 1._QP/3._QP  ), R8P)
   res(4) = real(abs(sum(b*ac)        - 1._QP/6._QP  ), R8P)
   res(5) = real(abs(sum(b*c*c*c)     - 1._QP/4._QP  ), R8P)
   res(6) = real(abs(sum(b*c*ac)      - 1._QP/8._QP  ), R8P)
   res(7) = real(abs(sum(b*ac2)       - 1._QP/12._QP ), R8P)
   res(8) = real(abs(sum(b*aac)       - 1._QP/24._QP ), R8P)
   endfunction order_residuals

   pure function max_order_residual(bt, p) result(r)
   !< Largest residual among the order conditions up to order p (1 <= p <= 4).
   type(butcher_tableau), intent(in) :: bt !< Butcher tableau.
   integer(I4P),          intent(in) :: p  !< Order.
   real(R8P)                         :: r  !< Largest residual.
   real(R8P)                         :: res(8) !< All residuals.

   res = order_residuals(bt)
   r = maxval(res(1:CONDITIONS_UP_TO(p)))
   endfunction max_order_residual

   pure function rhs(problem, y) result(f)
   !< Right-hand side of the scalar test problems: 1, y' = -y; 2, y' = -y².
   integer(I4P), intent(in) :: problem !< Problem.
   real(QP),     intent(in) :: y       !< State.
   real(QP)                 :: f       !< y'.

   select case(problem)
   case(1)
      f = -y
   case default
      f = -y * y
   endselect
   endfunction rhs

   pure function exact(problem, t) result(y)
   !< Exact solutions from y(0) = 1: e^{-t} and 1/(1 + t).
   integer(I4P), intent(in) :: problem !< Problem.
   real(QP),     intent(in) :: t       !< Time.
   real(QP)                 :: y       !< y(t).

   select case(problem)
   case(1)
      y = exp(-t)
   case default
      y = 1._QP / (1._QP + t)
   endselect
   endfunction exact

   pure function integrate_scalar(bt, problem, n) result(err)
   !< |error| at t = 1 of n steps of the explicit tableau on a scalar problem, in QP (coefficients converted from R_P).
   type(butcher_tableau), intent(in) :: bt      !< Explicit Butcher tableau.
   integer(I4P),          intent(in) :: problem !< Problem (see rhs).
   integer(I4P),          intent(in) :: n       !< Steps.
   real(R8P)                         :: err     !< |y_n - y(1)|.
   real(QP), allocatable             :: a(:,:), b(:), k(:) !< QP coefficients and stage residuals.
   real(QP)                          :: y, ys, dt          !< State, stage state, step.
   integer(I4P)                      :: step, i, j         !< Counters.

   a = real(bt%a, QP)
   b = real(bt%b, QP)
   allocate(k(bt%stages))
   dt = 1._QP / real(n, QP)
   y = 1._QP
   do step=1, n
      do i=1, bt%stages
         ys = y
         do j=1, i-1
            ys = ys + dt * a(i,j) * k(j)
         enddo
         k(i) = rhs(problem, ys)
      enddo
      do i=1, bt%stages
         y = y + dt * b(i) * k(i)
      enddo
   enddo
   err = real(abs(y - exact(problem, 1._QP)), R8P)
   endfunction integrate_scalar

   pure function observed_orders(errors) result(orders)
   !< Observed orders log2(e_k / e_{k+1}) of a sequence of errors at halved steps.
   real(R8P), intent(in) :: errors(:)                 !< Errors at n, 2n, 4n, ...
   real(R8P)             :: orders(size(errors) - 1)  !< Observed orders.

   orders = log(errors(1:size(errors)-1) / errors(2:)) / log(2._R8P)
   endfunction observed_orders

   subroutine check(condition, label)
   !< Print PASS/FAIL for a check and count the failures.
   logical,      intent(in) :: condition !< Check result.
   character(*), intent(in) :: label     !< What is checked.

   if (condition) then
      print '(A)', '  pass  '//label
   else
      print '(A)', '  FAIL  '//label
      failures = failures + 1_I4P
   endif
   endsubroutine check
endmodule selah_test_tools
