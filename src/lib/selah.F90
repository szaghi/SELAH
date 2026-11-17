!< SELAH, Stage-wise Evolution Library for Accelerators and HPC: the library's entry point.
module selah
!< SELAH, Stage-wise Evolution Library for Accelerators and HPC: the library's entry point.
!<
!< `use selah` gives the working kinds (R_P, the operators' kind, chosen at build time: selah_kind_parameters). The
!< scheme tableaux (L0), the kernels (L1) and the stepper (L2) join this module as they land (notes/PLAN.md).
use penf, only : I4P, I8P, R4P, R8P
use selah_kind_parameters
implicit none
private
public :: I4P, I8P, R4P, R8P, R_P
endmodule selah
