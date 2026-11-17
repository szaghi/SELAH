!< SELAH, L1: the box kernels, the stage operations over the box of a rank-5 field.

#include "fundal.H"

! the dummy arrays carry the lower bounds of the allocation (ghosts included): no rebasing to 1
#define SELAH_BOUNDS bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):
! the mask joins the data clauses of the directives only on the masked loops; dq only when it is computed
#define SELAH_MASK_ARG
#define SELAH_DQ_ARG

module selah_kernels_box
!< SELAH, L1: the box kernels, the stage operations over the box of a rank-5 field.
!<
!< Every kernel is a pointwise update of the cells of the box (selah_box): never a ghost cell, never a pool slot past
!< the active blocks. One source serves the host (plain loops, OpenMP with `-fopenmp`), OpenACC and OpenMP offload, by
!< the backend macros of FUNDAL (`fundal.H`) in the directive style of ADAM's FNL kernels; in a device build every array
!< argument, the coefficients included, must be device resident.
!<
!< **Mask.** With `mask` present, a cell is active where `mask(..., mask_comp) < 0` (ADAM's immersed-boundary distance,
!< negative in the fluid; `mask_comp` defaults to the last component, ADAM's all-solids summary). The mask has the
!< field's shape with its own number of components. Inactive cells are left untouched, except the step residual `dq`,
!< which is 0 there.
!<
!< **Summation order.** The sums over registers form a private sum and store once, as ADAM's FNL kernels do, so the
!< results are bitwise those of ADAM's device path (src/tests: the parity test).
!<
!< **Aliasing.** The output arrays must not alias the inputs (Fortran's rule on dummy arguments).
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_box, only : box, LAYOUT_V_FIRST, component_dim
implicit none
private
public :: copy_box
public :: axpy_box
public :: lincomb_box
public :: stage_state_box
public :: step_update_box
public :: shu_osher_2r_stage_box
public :: williamson_2n_stage_box

contains
   subroutine copy_box(bx, x, y, mask, mask_comp)
   !< y = x over the box.
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: x(SELAH_BOUNDS)     !< Source.
   real(R_P),    intent(inout)        :: y(SELAH_BOUNDS)     !< Destination.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_copy.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_copy.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_copy.INC"
   endif
   endsubroutine copy_box

   subroutine axpy_box(bx, alpha, x, y, mask, mask_comp)
   !< y = y + alpha x over the box (a partition update of a splitting scheme, a forward-Euler step).
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: alpha               !< Coefficient.
   real(R_P),    intent(in)           :: x(SELAH_BOUNDS)     !< Increment.
   real(R_P),    intent(inout)        :: y(SELAH_BOUNDS)     !< Updated field.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_axpy.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_axpy.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_axpy.INC"
   endif
   endsubroutine axpy_box

   subroutine lincomb_box(bx, alpha, x, beta, y, mask, mask_comp)
   !< y = alpha x + beta y over the box.
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: alpha               !< Coefficient of x.
   real(R_P),    intent(in)           :: x(SELAH_BOUNDS)     !< First field.
   real(R_P),    intent(in)           :: beta                !< Coefficient of y.
   real(R_P),    intent(inout)        :: y(SELAH_BOUNDS)     !< Second field, overwritten.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_lincomb.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_lincomb.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_lincomb.INC"
   endif
   endsubroutine lincomb_box

   subroutine stage_state_box(bx, n, dt, coef, q, regs, u, mask, mask_comp)
   !< u = q + sum_{k=1..n} dt coef(k) regs(k) over the box: the state of a stage of a Butcher form, from the solution and
   !< the residuals of the previous stages (coef = row s of A, n = s - 1). n = 0 gives u = q, through the copy kernel:
   !< the summing loop with nothing to sum ran at 0.57 of the copy's speed on the A30 (M5, job 1368098).
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

   if (n == 0_I4P) then
      call copy_box(bx, q, u, mask, mask_comp)
      return
   endif
#include "selah_kernel_setup.INC"
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_stage_state.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_stage_state.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_stage_state.INC"
   endif
   endsubroutine stage_state_box

   subroutine step_update_box(bx, n, dt, coef, regs, q, dq, mask, mask_comp)
   !< q = q + dt sum_{k=1..n} coef(k) regs(k) over the box: the end of a step of a Butcher form (coef = b, n = s), fused
   !< with the step residual dq = sum_k coef(k) regs(k) when dq is present (0 where the mask is inactive).
   type(box),    intent(in)              :: bx                    !< Box.
   integer(I4P), intent(in)              :: n                     !< Registers summed.
   real(R_P),    intent(in)              :: dt                    !< Time step.
   real(R_P),    intent(in)              :: coef(1:)              !< Weights [n or more].
   real(R_P),    intent(in)              :: regs(SELAH_BOUNDS,1:) !< Residual registers [..., n or more].
   real(R_P),    intent(inout)           :: q(SELAH_BOUNDS)       !< Solution, advanced.
   real(R_P),    intent(inout), optional :: dq(SELAH_BOUNDS)      !< Step residual.
   real(R_P),    intent(in),    optional :: mask(SELAH_BOUNDS)    !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in),    optional :: mask_comp             !< Mask component (default: the last).
   real(R_P)                             :: d                     !< Private sum.
   integer(I4P)                          :: k                     !< Counter.
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
   if (present(dq)) then
#define SELAH_WITH_DQ
#undef SELAH_DQ_ARG
#define SELAH_DQ_ARG , dq
      if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
         if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_step_update.INC"
#undef SELAH_ACTIVE
         else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_step_update.INC"
#undef SELAH_ACTIVE
         endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
      else
#include "selah_kernel_step_update.INC"
      endif
#undef SELAH_WITH_DQ
#undef SELAH_DQ_ARG
#define SELAH_DQ_ARG
   else
      if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
         if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_step_update.INC"
#undef SELAH_ACTIVE
         else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_step_update.INC"
#undef SELAH_ACTIVE
         endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
      else
#include "selah_kernel_step_update.INC"
      endif
   endif
   endsubroutine step_update_box

   subroutine shu_osher_2r_stage_box(bx, dt, a, b, c, q0, r, q, mask, mask_comp)
   !< q = a q0 + b q + dt c r over the box: a stage of a two-register Shu–Osher form (q0 = u^n, r = F(q)).
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
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_shu_osher_2r_stage.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_shu_osher_2r_stage.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_shu_osher_2r_stage.INC"
   endif
   endsubroutine shu_osher_2r_stage_box

   subroutine williamson_2n_stage_box(bx, dt, a, b, r, dqr, q, mask, mask_comp)
   !< dqr = a dqr + dt r; q = q + b dqr over the box: a stage of a Williamson 2N form (r = F(q)). With a = 0 dqr is not
   !< read (dqr = dt r): the first stage needs no prior content of the register.
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: dt                  !< Time step.
   real(R_P),    intent(in)           :: a, b                !< Stage coefficients A_i, B_i.
   real(R_P),    intent(in)           :: r(SELAH_BOUNDS)     !< Residual of the current stage.
   real(R_P),    intent(inout)        :: dqr(SELAH_BOUNDS)   !< The Delta Q register.
   real(R_P),    intent(inout)        :: q(SELAH_BOUNDS)     !< Current stage, advanced.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).
#include "selah_kernel_decl.INC"

#include "selah_kernel_setup.INC"
   if (a == 0._R_P) then
#define SELAH_A_ZERO
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_williamson_2n_stage.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_williamson_2n_stage.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_williamson_2n_stage.INC"
   endif
#undef SELAH_A_ZERO
   else
   if (present(mask)) then
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
      if (bx%layout == LAYOUT_V_FIRST) then
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#include "selah_kernel_williamson_2n_stage.INC"
#undef SELAH_ACTIVE
      else
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#include "selah_kernel_williamson_2n_stage.INC"
#undef SELAH_ACTIVE
      endif
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
   else
#include "selah_kernel_williamson_2n_stage.INC"
   endif
   endif
   endsubroutine williamson_2n_stage_box
endmodule selah_kernels_box
