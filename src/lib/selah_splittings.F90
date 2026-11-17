!< SELAH, L0: splitting (composition) schemes for u' = F_1(u) + F_2(u).
module selah_splittings
!< SELAH, L0: splitting (composition) schemes for u' = F_1(u) + F_2(u).
!<
!< A splitting scheme advances one step as the product of the sub-flows of the two partitions, each over a fraction of
!< the step: \(\varphi_{c_m \Delta t}^{[p_m]} \circ \dots \circ \varphi_{c_1 \Delta t}^{[p_1]}\), applied left to
!< right in the arrays (flow 1 first). With explicit Euler sub-steps on each partition, as ADAM's PRISM does for
!< Maxwell's B and D, a sub-flow is exact when that partition's residual does not depend on its own variables.
!<
!< Every scheme here is a symmetric (palindromic) composition, so its order is even. It is stored as the half
!< sequence ending with the middle flow; `splitting` expands it:
!<
!<| Scheme | Flows | Order | Source |
!<|--------|-------|-------|--------|
!<| `strang`          | 3  | 2 | Strang (1968) |
!<| `yoshida-4`       | 7  | 4 | Yoshida (1990), the triple jump |
!<| `blanes-moan-s6`  | 13 | 4 | Blanes & Moan (2002, J. Comput. Appl. Math. 142, 313), S6 |
!<| `blanes-moan-s10` | 21 | 6 | Blanes & Moan (2002), S10 |
!<
!< The Blanes–Moan digits are those of ADAM's `blanesmoan_object` (ADAM issue #128, D3, has the composition ADAM's loop
!< got wrong).
use penf, only : I4P
use selah_kind_parameters, only : R_P
implicit none
private
public :: splitting_composition
public :: SCHEME_STRANG, SCHEME_YOSHIDA_4, SCHEME_BLANES_MOAN_S6, SCHEME_BLANES_MOAN_S10
public :: SPLITTING_SCHEMES
public :: splitting
public :: palindrome

character(*), parameter :: SCHEME_STRANG          = 'strang'          !< Strang, order 2.
character(*), parameter :: SCHEME_YOSHIDA_4       = 'yoshida-4'       !< Yoshida triple jump, order 4.
character(*), parameter :: SCHEME_BLANES_MOAN_S6  = 'blanes-moan-s6'  !< Blanes–Moan S6, order 4.
character(*), parameter :: SCHEME_BLANES_MOAN_S10 = 'blanes-moan-s10' !< Blanes–Moan S10, order 6.
character(15), parameter :: SPLITTING_SCHEMES(4) = [character(15) :: SCHEME_STRANG, SCHEME_YOSHIDA_4,          &
                                                    SCHEME_BLANES_MOAN_S6, SCHEME_BLANES_MOAN_S10] !< All of them.

type :: splitting_composition
   !< A composition of the sub-flows of two partitions.
   integer(I4P)              :: flows = 0_I4P !< Number of sub-flows.
   integer(I4P)              :: order = 0_I4P !< Declared order.
   integer(I4P), allocatable :: partition(:)  !< Partition of each sub-flow, 1 or 2 [flows].
   real(R_P),    allocatable :: coefficient(:) !< Fraction of the step of each sub-flow [flows].
endtype splitting_composition

contains
   pure function palindrome(first, half, order) result(sc)
   !< Symmetric composition from its half sequence: the flows alternate partitions from `first`, and the last entry of
   !< `half` is the middle flow.
   integer(I4P), intent(in)    :: first   !< Partition of the first flow, 1 or 2.
   real(R_P),    intent(in)    :: half(:) !< Coefficients up to and including the middle flow.
   integer(I4P), intent(in)    :: order   !< Declared order.
   type(splitting_composition) :: sc      !< Composition.
   integer(I4P)                :: m, h, i !< Flows, half length, counter.

   h = size(half, kind=I4P)
   m = 2_I4P * h - 1_I4P
   sc%flows = m
   sc%order = order
   allocate(sc%partition(m), sc%coefficient(m))
   do i=1, h
      sc%coefficient(i) = half(i)
      sc%coefficient(m+1-i) = half(i)
   enddo
   do i=1, m
      sc%partition(i) = merge(first, 3_I4P - first, mod(i, 2_I4P) == 1_I4P)
   enddo
   endfunction palindrome

   pure function splitting(scheme, first) result(sc)
   !< Splitting scheme by name; `first` (default 1) is the partition of the first flow.
   character(*), intent(in)           :: scheme !< Scheme name.
   integer(I4P), intent(in), optional :: first  !< Partition of the first flow.
   type(splitting_composition)        :: sc     !< Composition (0 flows if the name is unknown).
   integer(I4P)                       :: f      !< First partition.
   real(R_P)                          :: w0, w1, cbrt2 !< Yoshida coefficients.
   real(R_P)                          :: a(6), b(5)    !< Blanes–Moan coefficients.

   f = 1_I4P ; if (present(first)) f = first
   select case(trim(adjustl(scheme)))
   case(SCHEME_STRANG)
      sc = palindrome(f, [0.5_R_P, 1._R_P], 2_I4P)
   case(SCHEME_YOSHIDA_4)
      ! w1 = 1/(2 - 2^(1/3)), w0 = -2^(1/3) w1: the triple jump w1, w0, w1 of the Strang scheme
      cbrt2 = 2._R_P**(1._R_P/3._R_P)
      w1 = 1._R_P / (2._R_P - cbrt2)
      w0 = -cbrt2 * w1
      sc = palindrome(f, [0.5_R_P*w1, w1, 0.5_R_P*(w0 + w1), w0], 4_I4P)
   case(SCHEME_BLANES_MOAN_S6)
      a(1) =  0.0792036964311957_R_P
      a(2) =  0.3531729060497740_R_P
      a(3) = -0.0420650803577195_R_P
      a(4) =  1._R_P - 2._R_P*(a(1) + a(2) + a(3))
      b(1) =  0.209515106613362_R_P
      b(2) = -0.143851773179818_R_P
      b(3) =  0.5_R_P - (b(1) + b(2))
      sc = palindrome(f, [a(1), b(1), a(2), b(2), a(3), b(3), a(4)], 4_I4P)
   case(SCHEME_BLANES_MOAN_S10)
      a(1) =  0.0502627644003922_R_P
      a(2) =  0.4135143004283440_R_P
      a(3) =  0.0450798897943977_R_P
      a(4) = -0.1880548538195690_R_P
      a(5) =  0.5419606784507800_R_P
      a(6) =  1._R_P - 2._R_P*(a(1) + a(2) + a(3) + a(4) + a(5))
      b(1) =  0.148816447901042_R_P
      b(2) = -0.132385865767784_R_P
      b(3) =  0.067307604692185_R_P
      b(4) =  0.432666402578175_R_P
      b(5) =  0.5_R_P - (b(1) + b(2) + b(3) + b(4))
      sc = palindrome(f, [a(1), b(1), a(2), b(2), a(3), b(3), a(4), b(4), a(5), b(5), a(6)], 6_I4P)
   endselect
   endfunction splitting
endmodule selah_splittings
