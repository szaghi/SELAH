   ! One variant of the box kernel stage_state of SELAH (selah_kernels_box) as a procedure of its own: SELAH_PROC is its
   ! name (undefined at the end), SELAH_ACTIVE, SELAH_FILL and SELAH_MASK_ARG select the variant (see the private
   ! procedures of selah_kernels_box). The dummy arguments are those of stage_state_box, without `masked_q` (it selects
   ! SELAH_FILL, here at compile time); n > 0 (the dispatcher handles n = 0).
   subroutine SELAH_PROC(bx, n, dt, coef, q, regs, u, mask, mask_comp)
   type(box),    intent(in)           :: bx                   !< Box.
   integer(I4P), intent(in)           :: n                    !< Registers summed.
   real(R_P),    intent(in)           :: dt                   !< Time step.
   real(R_P),    intent(in)           :: coef(1:)             !< Coefficients [n or more].
   real(R_P),    intent(in)           :: q(SELAH_BOUNDS)      !< Solution.
   real(R_P),    intent(in)           :: regs(SELAH_BOUNDS,1:) !< Residual registers [..., n or more].
   real(R_P),    intent(inout)        :: u(SELAH_BOUNDS)      !< Stage state.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)   !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp            !< Mask component (default: the last).
   real(R_P)                          :: d                    !< Private sum.
   integer(I4P)                       :: k                    !< Counter.
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
#include "selah_kernel_stage_state.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
