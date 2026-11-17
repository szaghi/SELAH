!< SELAH, the working real kind.
module selah_kind_parameters
!< SELAH, the working real kind.
!<
!< Every operator computes in PENF's R_P, chosen when SELAH (and PENF, compiled from its sources by the same build mode)
!< is built: -D_R_P_IS_R4P, -D_R_P_IS_R8P, or nothing (R8P, PENF's default). One precision per build, no generic
!< interfaces: the kernels are compiled once, for that kind (notes/PLAN.md, principle P5).
!<
!< A compile-time check ties the macro to PENF's R_P: a build that defines it for SELAH but not for PENF (or the
!< reverse) stops with a division by zero in a constant expression.
use penf, only : R_P, R4P, R8P
implicit none
private
public :: R_P

#if defined _R_P_IS_R4P
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R4P) !< Stops the build if PENF's R_P is not R4P.
#elif defined _R_P_IS_R16P
#error "SELAH: R16P is not supported (no device code, no caller)"
#else
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R8P) !< Stops the build if PENF's R_P is not R8P.
#endif
endmodule selah_kind_parameters
