   ! One variant of the box kernel shu_osher_2r_stage of SELAH (selah_kernels_box) as a procedure of its own: SELAH_PROC
   ! is its name (undefined at the end), SELAH_ACTIVE and SELAH_MASK_ARG select the variant (see the private procedures
   ! of selah_kernels_box). The dummy arguments are those of shu_osher_2r_stage_box.
   subroutine SELAH_PROC(bx, dt, a, b, c, q0, r, q, mask, mask_comp)
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: dt                  !< Time step.
   real(R_P),    intent(in)           :: a, b, c             !< Stage coefficients.
   real(R_P),    intent(in)           :: q0(SELAH_BOUNDS)    !< Solution at the start of the step.
   real(R_P),    intent(in)           :: r(SELAH_BOUNDS)     !< Residual of the current stage.
   real(R_P),    intent(inout)        :: q(SELAH_BOUNDS)     !< Current stage, advanced.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
#include "selah_kernel_shu_osher_2r_stage.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
