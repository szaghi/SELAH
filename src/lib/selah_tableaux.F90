!< SELAH, L0: the coefficient forms of a Runge–Kutta scheme and their conversion to the Butcher tableau.
module selah_tableaux
!< SELAH, L0: the coefficient forms of a Runge–Kutta scheme and their conversion to the Butcher tableau.
!<
!< Plain data, no behaviour beyond pure queries: the kernels (L1) and the stepper (L2) read these coefficients, the
!< tests verify them. Every form converts to the Butcher tableau, the form the order conditions and the reflux weights
!< (`b`) are stated in. For an autonomous system \(u' = F(u)\) the forms are:
!<
!<+ **Butcher**: \(Y_i = u^n + \Delta t \sum_j a_{ij} F(Y_j)\), \(u^{n+1} = u^n + \Delta t \sum_j b_j F(Y_j)\),
!<  \(c = A \mathbf{1}\);
!<+ **Shu–Osher**: \(Y_1 = u^n\), \(Y_i = \sum_{k<i} \left(\alpha_{ik} Y_k + \Delta t \beta_{ik} F(Y_k)\right)\) for
!<  \(i = 2, \dots, s+1\), \(u^{n+1} = Y_{s+1}\), each row of \(\alpha\) summing to 1;
!<+ **two-register Shu–Osher**: \(Y_{i+1} = a_i Y_1 + b_i Y_i + c_i \Delta t F(Y_i)\), the special case that needs
!<  only \(u^n\) and the current stage (ADAM's low-storage `ark`, `brk`, `crk`);
!<+ **Williamson 2N**: \(\Delta Q_i = A_i \Delta Q_{i-1} + \Delta t F(Y_i)\), \(Y_{i+1} = Y_i + B_i \Delta Q_i\),
!<  \(A_1 = 0\), \(Y_1 = u^n\), \(u^{n+1} = Y_{s+1}\) (Williamson 1980).
!<
!< The tableau `a` is not assumed strictly lower triangular: `is_explicit` tells (notes/PLAN.md, "Beyond v0.x").
use penf, only : I4P
use selah_kind_parameters, only : R_P
implicit none
private
public :: butcher_tableau
public :: shu_osher_form
public :: shu_osher_2r_form
public :: williamson_2n_form
public :: butcher_from_shu_osher
public :: butcher_from_shu_osher_2r
public :: butcher_from_williamson_2n
public :: shu_osher_from_2r
public :: is_explicit
public :: ssp_coefficient
public :: euler_step_fractions
public :: stability_polynomial
public :: imag_axis_bound
public :: real_axis_bound

type :: butcher_tableau
   !< Butcher tableau of an s-stage Runge–Kutta scheme.
   integer(I4P)           :: stages = 0_I4P          !< Stages number s.
   integer(I4P)           :: order  = 0_I4P          !< Declared order of accuracy.
   real(R_P), allocatable :: a(:,:)                  !< Matrix A [s, s].
   real(R_P), allocatable :: b(:)                    !< Weights b [s].
   real(R_P), allocatable :: c(:)                    !< Abscissas c [s].
   real(R_P), allocatable :: b_hat(:)                !< Embedded weights [s], unallocated when the scheme has none.
   real(R_P)              :: ssp_coefficient = 0._R_P !< SSP coefficient, 0 when the scheme is not SSP.
endtype butcher_tableau

type :: shu_osher_form
   !< Shu–Osher form: `alpha(i,k)` and `beta(i,k)` for i = 2..s+1, k = 1..i-1; row 1 and the upper triangle are zero.
   integer(I4P)           :: stages = 0_I4P !< Stages number s.
   real(R_P), allocatable :: alpha(:,:)     !< State coefficients [s+1, s].
   real(R_P), allocatable :: beta(:,:)      !< Residual coefficients [s+1, s].
endtype shu_osher_form

type :: shu_osher_2r_form
   !< Two-register Shu–Osher form: stage i+1 from u^n and stage i only.
   integer(I4P)           :: stages = 0_I4P !< Stages number s.
   real(R_P), allocatable :: a(:)           !< Coefficient of u^n [s].
   real(R_P), allocatable :: b(:)           !< Coefficient of the current stage [s].
   real(R_P), allocatable :: c(:)           !< Coefficient of dt F(current stage) [s].
endtype shu_osher_2r_form

type :: williamson_2n_form
   !< Williamson 2N low-storage form.
   integer(I4P)           :: stages = 0_I4P !< Stages number s.
   real(R_P), allocatable :: a(:)           !< A_i [s], A_1 = 0.
   real(R_P), allocatable :: b(:)           !< B_i [s].
endtype williamson_2n_form

contains
   pure function butcher_from_shu_osher(so, order) result(bt)
   !< Butcher tableau of a Shu–Osher form.
   !<
   !< Each stage is written as \(Y_i = u^n + \Delta t \sum_j B_{ij} F(Y_j)\): row i of B is the alpha-weighted sum of
   !< the rows k < i plus `beta(i,k)` on column k. Rows 1..s give A, row s+1 gives b, and c = A·1. The summation order
   !< is that of ADAM's `ssp_tableau` (k ascending, the alpha term before the beta term), so a Shu–Osher form with the
   !< same coefficients gives the same tableau bitwise.
   type(shu_osher_form),   intent(in)           :: so    !< Shu–Osher form.
   integer(I4P),           intent(in), optional :: order !< Declared order.
   type(butcher_tableau)                        :: bt    !< Butcher tableau.
   real(R_P), allocatable                       :: bb(:,:) !< Stage coefficients [s+1, s].
   integer(I4P)                                 :: i, k    !< Counters.

   associate(s=>so%stages)
   allocate(bb(s+1,s))
   bb = 0._R_P
   do i=2, s+1
      do k=1, i-1
         bb(i,:) = bb(i,:) + so%alpha(i,k) * bb(k,:)
         bb(i,k) = bb(i,k) + so%beta(i,k)
      enddo
   enddo
   bt%stages = s
   ! allocated before the assignments, as ADAM does: nvfortran 26.9 leaves an unallocated allocatable component
   ! unallocated when assigned the result of sum(..., dim=2)
   allocate(bt%a(s,s), bt%b(s), bt%c(s))
   bt%a = bb(1:s,:)
   bt%b = bb(s+1,:)
   bt%c = sum(bt%a, dim=2)
   endassociate
   if (present(order)) bt%order = order
   bt%ssp_coefficient = ssp_coefficient(so)
   endfunction butcher_from_shu_osher

   pure function shu_osher_from_2r(r) result(so)
   !< Shu–Osher form of a two-register form: alpha(i+1,1) = a_i, alpha(i+1,i) += b_i, beta(i+1,i) = c_i.
   type(shu_osher_2r_form), intent(in) :: r  !< Two-register form.
   type(shu_osher_form)                :: so !< Shu–Osher form.
   integer(I4P)                        :: i  !< Counter.

   associate(s=>r%stages)
   so%stages = s
   allocate(so%alpha(s+1,s), so%beta(s+1,s))
   so%alpha = 0._R_P
   so%beta = 0._R_P
   do i=1, s
      so%alpha(i+1,1) = so%alpha(i+1,1) + r%a(i)
      so%alpha(i+1,i) = so%alpha(i+1,i) + r%b(i)
      so%beta(i+1,i) = r%c(i)
   enddo
   endassociate
   endfunction shu_osher_from_2r

   pure function butcher_from_shu_osher_2r(r, order) result(bt)
   !< Butcher tableau of a two-register Shu–Osher form.
   type(shu_osher_2r_form), intent(in)           :: r     !< Two-register form.
   integer(I4P),            intent(in), optional :: order !< Declared order.
   type(butcher_tableau)                         :: bt    !< Butcher tableau.

   bt = butcher_from_shu_osher(shu_osher_from_2r(r), order)
   endfunction butcher_from_shu_osher_2r

   pure function butcher_from_williamson_2n(w, order) result(bt)
   !< Butcher tableau of a Williamson 2N form.
   !<
   !< \(\Delta Q_i = \sum_j D_{ij} \Delta t F(Y_j)\) with \(D_i = A_i D_{i-1} + e_i\); the stage rows follow
   !< \(Y_{i+1} = Y_i + B_i D_i\). Rows 1..s give A, row s+1 gives b, and c = A·1.
   type(williamson_2n_form), intent(in)           :: w       !< Williamson 2N form.
   integer(I4P),             intent(in), optional :: order   !< Declared order.
   type(butcher_tableau)                          :: bt      !< Butcher tableau.
   real(R_P), allocatable                         :: yy(:,:) !< Stage coefficients [s+1, s].
   real(R_P), allocatable                         :: dd(:)   !< Coefficients of Delta Q [s].
   integer(I4P)                                   :: i       !< Counter.

   associate(s=>w%stages)
   allocate(yy(s+1,s), dd(s))
   yy = 0._R_P
   dd = 0._R_P
   do i=1, s
      dd = w%a(i) * dd
      dd(i) = dd(i) + 1._R_P
      yy(i+1,:) = yy(i,:) + w%b(i) * dd
   enddo
   bt%stages = s
   allocate(bt%a(s,s), bt%b(s), bt%c(s)) ! before the assignments (see butcher_from_shu_osher)
   bt%a = yy(1:s,:)
   bt%b = yy(s+1,:)
   bt%c = sum(bt%a, dim=2)
   endassociate
   if (present(order)) bt%order = order
   endfunction butcher_from_williamson_2n

   pure function is_explicit(bt) result(yes)
   !< True when the tableau is strictly lower triangular.
   type(butcher_tableau), intent(in) :: bt  !< Butcher tableau.
   logical                           :: yes !< Explicit scheme.
   integer(I4P)                      :: i   !< Counter.

   yes = .true.
   do i=1, bt%stages
      if (any(bt%a(i,i:) /= 0._R_P)) then
         yes = .false.
         return
      endif
   enddo
   endfunction is_explicit

   pure function ssp_coefficient(so) result(ssp)
   !< SSP coefficient of a Shu–Osher form: \(\min \alpha_{ik} / \beta_{ik}\) over \(\beta_{ik} > 0\).
   !<
   !< 0 when a coefficient is negative, or when \(\beta_{ik} > 0\) with \(\alpha_{ik} = 0\): the form is then not a
   !< convex combination of forward-Euler steps. It is the coefficient of this representation, a lower bound of the
   !< scheme's optimal one; the catalogue stores the optimal representations.
   type(shu_osher_form), intent(in) :: so  !< Shu–Osher form.
   real(R_P)                        :: ssp !< SSP coefficient.
   integer(I4P)                     :: i, k !< Counters.

   ssp = 0._R_P
   if (any(so%alpha < 0._R_P) .or. any(so%beta < 0._R_P)) return
   ssp = huge(1._R_P)
   do i=2, so%stages+1
      do k=1, i-1
         if (so%beta(i,k) > 0._R_P) then
            if (so%alpha(i,k) <= 0._R_P) then
               ssp = 0._R_P
               return
            endif
            ssp = min(ssp, so%alpha(i,k) / so%beta(i,k))
         endif
      enddo
   enddo
   if (ssp == huge(1._R_P)) ssp = 0._R_P
   endfunction ssp_coefficient
   pure function euler_step_fractions(so) result(f)
   !< Per stage k, the largest forward-Euler step its residual takes in a Shu–Osher form, as a fraction of the step:
   !< \(f_k = \max_{i > k,\ \beta_{ik} > 0} \beta_{ik} / \alpha_{ik}\).
   !<
   !< Stage i of the form is the convex combination \(Y_i = \sum_{k<i} \alpha_{ik} (Y_k + \Delta t\, (\beta_{ik} /
   !< \alpha_{ik}) F(Y_k))\): the residual of stage k enters forward-Euler steps of \(\Delta t\, \beta_{ik} /
   !< \alpha_{ik}\). A limiter that makes \(Y_k + \tau F(Y_k)\) admissible (positive, bounded) for \(\tau = f_k
   !< \Delta t\) makes it so for every \(0 \le \tau \le f_k \Delta t\) (a convex combination of \(Y_k\) and that
   !< state), hence every later stage and the step. \(\max_k f_k = 1 / C\), C the SSP coefficient of the form.
   !< Every \(f_k\) is 0 when the form is not a convex combination of forward-Euler steps (a negative coefficient, or
   !< \(\beta_{ik} > 0\) with \(\alpha_{ik} = 0\): `ssp_coefficient` is 0 then), and \(f_k = 0\) for a stage whose
   !< residual no later stage uses.
   type(shu_osher_form), intent(in) :: so             !< Shu–Osher form.
   real(R_P)                        :: f(so%stages)   !< Forward-Euler fractions [s].
   integer(I4P)                     :: i, k           !< Counters.

   f = 0._R_P
   if (ssp_coefficient(so) <= 0._R_P) return
   do k=1, so%stages
      do i=k+1, so%stages+1
         if (so%beta(i,k) > 0._R_P) f(k) = max(f(k), so%beta(i,k) / so%alpha(i,k))
      enddo
   enddo
   endfunction euler_step_fractions
   pure function stability_polynomial(bt) result(g)
   !< Coefficients of the stability polynomial of an explicit tableau, \(R(z) = \sum_{k=0}^{s} g_k z^k\):
   !< \(g_0 = 1\), \(g_k = b^T A^{k-1} \mathbf 1\).
   type(butcher_tableau), intent(in) :: bt                 !< Butcher tableau (explicit).
   real(R_P)                         :: g(0:bt%stages)     !< Coefficients.
   real(R_P)                         :: v(bt%stages)       !< A^{k-1} 1.
   integer(I4P)                      :: k                  !< Counter.

   g(0) = 1._R_P
   v = 1._R_P
   do k=1, bt%stages
      g(k) = sum(bt%b * v)
      v = matmul(bt%a, v)
   enddo
   endfunction stability_polynomial

   pure function imag_axis_bound(bt) result(beta)
   !< The imaginary-axis stability bound: the largest \(\beta\) with \(|R(iy)| \le 1\) for \(|y| \le \beta\), 0 when
   !< none (euler, ssp-22). The step of an oscillatory mode of frequency \(\omega\) (a gyration, a wave) is stable
   !< for \(\omega \Delta t \le \beta\).
   !<
   !< \(P(w) = |R(iy)|^2 - 1\) is a polynomial in \(w = y^2\); its coefficients of \(w^j\), \(2j \le p\), vanish
   !< for a scheme of order p (\(R\) agrees with \(e^z\) to \(z^p\)) and are set to 0, so that the rounding of the
   !< tableau's digits does not decide the sign near 0. \(\beta = 0\) when the lowest remaining coefficient is
   !< positive, else \(\sqrt{w^*}\), \(w^*\) the first point where \(P\) turns positive (a scan of \(y\) by
   !< 1e-3, then bisection to round-off).
   type(butcher_tableau), intent(in) :: bt          !< Butcher tableau (explicit).
   real(R_P)                         :: beta        !< Bound.
   real(R_P)                         :: g(0:bt%stages), pw(0:bt%stages) !< R's coefficients; P's in w.
   real(R_P)                         :: re(0:bt%stages), im(0:bt%stages) !< Re R(iy) and Im R(iy)/y as polynomials in w.
   real(R_P)                         :: lo, hi, mid !< Bisection.
   integer(I4P)                      :: j, k, s     !< Counters, stages.

   s = bt%stages
   g = stability_polynomial(bt)
   re = 0._R_P ; im = 0._R_P
   do k=0, s
      if (mod(k, 2_I4P) == 0_I4P) then
         re(k/2) = re(k/2) + (-1._R_P)**(k/2) * g(k)
      else
         im((k-1)/2) = im((k-1)/2) + (-1._R_P)**((k-1)/2) * g(k)
      endif
   enddo
   pw = 0._R_P
   do j=0, s
      do k=0, s - j
         pw(j+k) = pw(j+k) + re(j) * re(k)
         if (j + k + 1 <= s) pw(j+k+1) = pw(j+k+1) + im(j) * im(k)
      enddo
   enddo
   pw(0) = pw(0) - 1._R_P
   do j=0, s
      if (2 * j <= bt%order) pw(j) = 0._R_P
   enddo
   beta = 0._R_P
   do j=0, s
      if (pw(j) > 0._R_P) return
      if (pw(j) < 0._R_P) exit
   enddo
   if (all(pw == 0._R_P)) then
      beta = huge(1._R_P)
      return
   endif
   hi = 0._R_P
   do while (p_of(hi + 1.e-3_R_P) <= 0._R_P .and. hi < 100._R_P)
      hi = hi + 1.e-3_R_P
   enddo
   lo = hi ; hi = hi + 1.e-3_R_P
   do k=1, 100
      mid = 0.5_R_P * (lo + hi)
      if (p_of(mid) <= 0._R_P) then ; lo = mid ; else ; hi = mid ; endif
   enddo
   beta = lo
   contains
      pure function p_of(y) result(p)
      !< P(y^2), Horner.
      real(R_P), intent(in) :: y
      real(R_P)             :: p
      integer(I4P)          :: i

      p = pw(s)
      do i=s-1, 0, -1
         p = p * y**2 + pw(i)
      enddo
      endfunction p_of
   endfunction imag_axis_bound

   pure function real_axis_bound(bt) result(beta)
   !< The real-axis stability bound: the largest \(\beta\) with \(|R(-x)| \le 1\) for \(0 \le x \le \beta\) (a decay of
   !< rate \(\lambda\) is stable for \(\lambda \Delta t \le \beta\)); a scan by 1e-3, then bisection.
   type(butcher_tableau), intent(in) :: bt          !< Butcher tableau (explicit).
   real(R_P)                         :: beta        !< Bound.
   real(R_P)                         :: g(0:bt%stages) !< R's coefficients.
   real(R_P)                         :: lo, hi, mid !< Bisection.
   integer(I4P)                      :: k           !< Counter.

   g = stability_polynomial(bt)
   hi = 0._R_P
   do while (abs(r_of(-(hi + 1.e-3_R_P))) <= 1._R_P .and. hi < 100._R_P)
      hi = hi + 1.e-3_R_P
   enddo
   lo = hi ; hi = hi + 1.e-3_R_P
   do k=1, 100
      mid = 0.5_R_P * (lo + hi)
      if (abs(r_of(-mid)) <= 1._R_P) then ; lo = mid ; else ; hi = mid ; endif
   enddo
   beta = lo
   contains
      pure function r_of(x) result(r)
      !< R(x), Horner.
      real(R_P), intent(in) :: x
      real(R_P)             :: r
      integer(I4P)          :: i

      r = g(bt%stages)
      do i=bt%stages-1, 0, -1
         r = r * x + g(i)
      enddo
      endfunction r_of
   endfunction real_axis_bound
endmodule selah_tableaux
