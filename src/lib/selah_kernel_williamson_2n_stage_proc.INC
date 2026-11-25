   ! One variant of the box kernel williamson_2n_stage of SELAH (selah_kernels_box) as a procedure of its own: SELAH_PROC
   ! is its name (undefined at the end), SELAH_ACTIVE, SELAH_MASK_ARG and SELAH_A_ZERO select the variant (see the private
   ! procedures of selah_kernels_box). The dummy arguments are those of williamson_2n_stage_box; `a` is a dummy of the
   ! variants without SELAH_A_ZERO only (with a = 0 the kernel never reads it).
#ifdef SELAH_A_ZERO
   subroutine SELAH_PROC(bx, dt, b, r, dqr, q, mask, mask_comp)
#else
   subroutine SELAH_PROC(bx, dt, a, b, r, dqr, q, mask, mask_comp)
#endif
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: dt                  !< Time step.
#ifndef SELAH_A_ZERO
   real(R_P),    intent(in)           :: a                   !< Stage coefficient A_i.
#endif
   real(R_P),    intent(in)           :: b                   !< Stage coefficient B_i.
   real(R_P),    intent(in)           :: r(SELAH_BOUNDS)     !< Residual of the current stage.
   real(R_P),    intent(inout)        :: dqr(SELAH_BOUNDS)   !< The Delta Q register.
   real(R_P),    intent(inout)        :: q(SELAH_BOUNDS)     !< Current stage, advanced.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
#include "selah_kernel_williamson_2n_stage.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
