   ! One variant of the box kernel step_update of SELAH (selah_kernels_box) as a procedure of its own: SELAH_PROC is its
   ! name (undefined at the end), SELAH_ACTIVE, SELAH_MASK_ARG, SELAH_WITH_DQ and SELAH_DQ_ARG select the variant (see
   ! the private procedures of selah_kernels_box). The dummy arguments are those of step_update_box; `dq` is a dummy of
   ! the variants with SELAH_WITH_DQ only.
#ifdef SELAH_WITH_DQ
   subroutine SELAH_PROC(bx, n, dt, coef, regs, q, dq, mask, mask_comp)
#else
   subroutine SELAH_PROC(bx, n, dt, coef, regs, q, mask, mask_comp)
#endif
   type(box),    intent(in)              :: bx                    !< Box.
   integer(I4P), intent(in)              :: n                     !< Registers summed.
   real(R_P),    intent(in)              :: dt                    !< Time step.
   real(R_P),    intent(in)              :: coef(1:)              !< Weights [n or more].
   real(R_P),    intent(in)              :: regs(SELAH_BOUNDS,1:) !< Residual registers [..., n or more].
   real(R_P),    intent(inout)           :: q(SELAH_BOUNDS)       !< Solution, advanced.
#ifdef SELAH_WITH_DQ
   real(R_P),    intent(inout)           :: dq(SELAH_BOUNDS)      !< Step residual.
#endif
   real(R_P),    intent(in),    optional :: mask(SELAH_BOUNDS)    !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in),    optional :: mask_comp             !< Mask component (default: the last).
   real(R_P)                             :: d                     !< Private sum.
   integer(I4P)                          :: k                     !< Counter.
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
#include "selah_kernel_step_update.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
