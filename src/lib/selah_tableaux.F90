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
endmodule selah_tableaux
