   ! One variant of the box kernel lincomb of SELAH (selah_kernels_box) as a procedure of its own: SELAH_PROC is its name
   ! (undefined at the end), SELAH_ACTIVE and SELAH_MASK_ARG select the variant (see the private procedures of
   ! selah_kernels_box). The dummy arguments are those of lincomb_box.
   subroutine SELAH_PROC(bx, alpha, x, beta, y, mask, mask_comp)
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: alpha               !< Coefficient of x.
   real(R_P),    intent(in)           :: x(SELAH_BOUNDS)     !< First field.
   real(R_P),    intent(in)           :: beta                !< Coefficient of y.
   real(R_P),    intent(inout)        :: y(SELAH_BOUNDS)     !< Second field, overwritten.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
#include "selah_kernel_lincomb.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
