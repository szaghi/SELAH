!< SELAH, Stage-wise Evolution Library for Accelerators and HPC: the library's entry point.
module selah
!< SELAH, Stage-wise Evolution Library for Accelerators and HPC: the library's entry point.
!<
!< `use selah` gives the working kinds (R_P, the operators' kind, chosen at build time: selah_kind_parameters) and L0,
!< the coefficients: the coefficient forms and their conversions (selah_tableaux), the Runge–Kutta catalogue
!< (selah_schemes), the splitting schemes (selah_splittings) and the buffers each storage form needs
!< (selah_storage); and L1, the kernels: the box (selah_box) and the stage operations over a box or a flat array
!< (selah_kernels); and L2, the stepper a caller drives stage by stage (selah_stepper).
use penf, only : I4P, I8P, R4P, R8P
use selah_kind_parameters
use selah_tableaux
use selah_schemes
use selah_splittings
use selah_storage
use selah_box
use selah_kernels
use selah_stepper
implicit none
private
! kinds
public :: I4P, I8P, R4P, R8P, R_P
! coefficient forms and conversions
public :: butcher_tableau, shu_osher_form, shu_osher_2r_form, williamson_2n_form
public :: butcher_from_shu_osher, butcher_from_shu_osher_2r, butcher_from_williamson_2n, shu_osher_from_2r
public :: is_explicit, ssp_coefficient
! Runge–Kutta catalogue
public :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54, SCHEME_SSP_104
public :: SCHEME_RK4, SCHEME_LSRK_33, SCHEME_LSRK_54, RK_SCHEMES
public :: butcher, shu_osher, shu_osher_2r, williamson_2n
! splittings
public :: splitting_composition, palindrome, splitting
public :: SCHEME_STRANG, SCHEME_YOSHIDA_4, SCHEME_BLANES_MOAN_S6, SCHEME_BLANES_MOAN_S10, SPLITTING_SCHEMES
! storage
public :: FORM_BUTCHER, FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_KETCHESON_2R, FORM_SPLITTING
public :: registers_needed
! L1: box and kernels
public :: box, LAYOUT_V_FIRST, LAYOUT_V_LAST, box_interior, box_range, box_components, box_dim, component_dim
public :: selah_copy, selah_axpy, selah_lincomb, selah_stage_state, selah_step_update
public :: selah_shu_osher_2r_stage, selah_williamson_2n_stage
! L2: stepper
public :: stepper
endmodule selah
