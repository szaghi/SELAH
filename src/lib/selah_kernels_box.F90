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
!< which is 0 there, and the stage state with `masked_q` (`stage_state_box`), which is q there.
!<
!< **Summation order.** The sums over registers form a private sum and store once, as ADAM's FNL kernels do, so the
!< results are bitwise those of ADAM's device path (src/tests: the parity test).
!<
!< **Aliasing.** The output arrays must not alias the inputs (Fortran's rule on dummy arguments).
!<
!< **One procedure per variant.** The loop of a kernel is an include file (`selah_kernel_<name>.INC`), compiled once per
!< variant: masked or not, once per layout, with the step residual or not, with the stage-state fill, with `a = 0`. Each
!< variant is a private procedure of its own, made from the wrapper `selah_kernel_<name>_proc.INC` under the name given
!< by the macro `SELAH_PROC` (the block of private procedures at the end of the module), and the public routines only
!< dispatch to them. The reason is a defect of AMD's flang (AFAR 23.2.1 and 24.4.0): it pairs the host launch of an
!< OpenMP target region with the device kernel of ANOTHER region when several regions of one procedure come from the
!< same line of an included file, so that on an AMD GPU a call ran the kernel of an arbitrary variant: the time
!< integrator did not advance the state (ADAM issue #155). The same included line in different procedures is
!< not affected. Do not include one kernel file more than once in a procedure.
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
public :: stage_current_box

contains
   ! ---- the kernels: each dispatches to the private procedure of its variant (the end of the module)
   subroutine copy_box(bx, x, y, mask, mask_comp)
   !< y = x over the box.
   type(box),    intent(in)           :: bx                  !< Box.
   real(R_P),    intent(in)           :: x(SELAH_BOUNDS)     !< Source.
   real(R_P),    intent(inout)        :: y(SELAH_BOUNDS)     !< Destination.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)  !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp           !< Mask component (default: the last).

   if (present(mask)) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call copy_box_masked_first(bx, x, y, mask, mask_comp)
      else
         call copy_box_masked_last(bx, x, y, mask, mask_comp)
      endif
   else
      call copy_box_plain(bx, x, y)
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

   if (present(mask)) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call axpy_box_masked_first(bx, alpha, x, y, mask, mask_comp)
      else
         call axpy_box_masked_last(bx, alpha, x, y, mask, mask_comp)
      endif
   else
      call axpy_box_plain(bx, alpha, x, y)
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

   if (present(mask)) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call lincomb_box_masked_first(bx, alpha, x, beta, y, mask, mask_comp)
      else
         call lincomb_box_masked_last(bx, alpha, x, beta, y, mask, mask_comp)
      endif
   else
      call lincomb_box_plain(bx, alpha, x, beta, y)
   endif
   endsubroutine lincomb_box

   subroutine stage_state_box(bx, n, dt, coef, q, regs, u, mask, mask_comp, masked_q)
   !< u = q + sum_{k=1..n} dt coef(k) regs(k) over the box: the state of a stage of a Butcher form, from the solution and
   !< the residuals of the previous stages (coef = row s of A, n = s - 1). n = 0 gives u = q, through the copy kernel:
   !< the summing loop with nothing to sum ran at 0.57 of the copy's speed on the A30 (M5, job 1368098). With
   !< `masked_q`, the inactive cells of u take q (the stepper's mask policy) in the same pass, instead of being left
   !< untouched.
   type(box),    intent(in)           :: bx                   !< Box.
   integer(I4P), intent(in)           :: n                    !< Registers summed.
   real(R_P),    intent(in)           :: dt                   !< Time step.
   real(R_P),    intent(in)           :: coef(1:)             !< Coefficients [n or more].
   real(R_P),    intent(in)           :: q(SELAH_BOUNDS)      !< Solution.
   real(R_P),    intent(in)           :: regs(SELAH_BOUNDS,1:) !< Residual registers [..., n or more].
   real(R_P),    intent(inout)        :: u(SELAH_BOUNDS)      !< Stage state.
   real(R_P),    intent(in), optional :: mask(SELAH_BOUNDS)   !< Active where mask(..., mask_comp) < 0.
   integer(I4P), intent(in), optional :: mask_comp            !< Mask component (default: the last).
   logical,      intent(in), optional :: masked_q             !< The inactive cells of u take q (default .false.).
   logical                            :: fill                 !< The inactive cells take q.

   fill = .false.
   if (present(masked_q)) fill = masked_q
   if (n == 0_I4P) then
      if (fill) then ! every cell of the box takes q
         call copy_box(bx, q, u)
      else
         call copy_box(bx, q, u, mask, mask_comp)
      endif
      return
   endif
   if (present(mask) .and. fill) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call stage_state_box_fill_first(bx, n, dt, coef, q, regs, u, mask, mask_comp)
      else
         call stage_state_box_fill_last(bx, n, dt, coef, q, regs, u, mask, mask_comp)
      endif
   elseif (present(mask)) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call stage_state_box_masked_first(bx, n, dt, coef, q, regs, u, mask, mask_comp)
      else
         call stage_state_box_masked_last(bx, n, dt, coef, q, regs, u, mask, mask_comp)
      endif
   else
      call stage_state_box_plain(bx, n, dt, coef, q, regs, u)
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

   if (present(dq)) then
      if (present(mask)) then
         if (bx%layout == LAYOUT_V_FIRST) then
            call step_update_box_dq_masked_first(bx, n, dt, coef, regs, q, dq, mask, mask_comp)
         else
            call step_update_box_dq_masked_last(bx, n, dt, coef, regs, q, dq, mask, mask_comp)
         endif
      else
         call step_update_box_dq_plain(bx, n, dt, coef, regs, q, dq)
      endif
   else
      if (present(mask)) then
         if (bx%layout == LAYOUT_V_FIRST) then
            call step_update_box_masked_first(bx, n, dt, coef, regs, q, mask, mask_comp)
         else
            call step_update_box_masked_last(bx, n, dt, coef, regs, q, mask, mask_comp)
         endif
      else
         call step_update_box_plain(bx, n, dt, coef, regs, q)
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

   if (present(mask)) then
      if (bx%layout == LAYOUT_V_FIRST) then
         call shu_osher_2r_stage_box_masked_first(bx, dt, a, b, c, q0, r, q, mask, mask_comp)
      else
         call shu_osher_2r_stage_box_masked_last(bx, dt, a, b, c, q0, r, q, mask, mask_comp)
      endif
   else
      call shu_osher_2r_stage_box_plain(bx, dt, a, b, c, q0, r, q)
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

   if (a == 0._R_P) then
      if (present(mask)) then
         if (bx%layout == LAYOUT_V_FIRST) then
            call williamson_2n_stage_box_a0_masked_first(bx, dt, b, r, dqr, q, mask, mask_comp)
         else
            call williamson_2n_stage_box_a0_masked_last(bx, dt, b, r, dqr, q, mask, mask_comp)
         endif
      else
         call williamson_2n_stage_box_a0_plain(bx, dt, b, r, dqr, q)
      endif
   else
      if (present(mask)) then
         if (bx%layout == LAYOUT_V_FIRST) then
            call williamson_2n_stage_box_masked_first(bx, dt, a, b, r, dqr, q, mask, mask_comp)
         else
            call williamson_2n_stage_box_masked_last(bx, dt, a, b, r, dqr, q, mask, mask_comp)
         endif
      else
         call williamson_2n_stage_box_plain(bx, dt, a, b, r, dqr, q)
      endif
   endif
   endsubroutine williamson_2n_stage_box

   subroutine stage_current_box(bx, n, coef, h, regs, y)
   !< y = (h - sum_{r=1..n-1} coef(r) regs(r)) / coef(n) over the box: the stage current that makes the weighted sum of
   !< the stage currents equal a given one (coef = row n+1 of A, or b at the last stage; selah_stepper `stage_current`).
   !< One variant, so one loop: it needs no procedure of its own.
   type(box),    intent(in)    :: bx                    !< Box.
   integer(I4P), intent(in)    :: n                     !< Stage (coef(n) /= 0).
   real(R_P),    intent(in)    :: coef(1:)              !< Coefficients [n or more].
   real(R_P),    intent(in)    :: h(SELAH_BOUNDS)       !< The weighted sum wanted.
   real(R_P),    intent(in)    :: regs(SELAH_BOUNDS,1:) !< The previous stage currents [..., n-1 or more].
   real(R_P),    intent(inout) :: y(SELAH_BOUNDS)       !< Stage current n.
   real(R_P)                   :: d                     !< Private sum.
   integer(I4P)                :: k                     !< Counter.
#include "selah_kernel_decl.INC"

   lo1 = bx%lo(1) ; lo2 = bx%lo(2) ; lo3 = bx%lo(3) ; lo4 = bx%lo(4) ; lo5 = bx%lo(5) ! no mask: the setup's range only
   hi1 = bx%hi(1) ; hi2 = bx%hi(2) ; hi3 = bx%hi(3) ; hi4 = bx%hi(4) ; hi5 = bx%hi(5)
#include "selah_kernel_stage_current.INC"
   endsubroutine stage_current_box

   ! ---- the private procedures: one per variant of a kernel (see the module documentation, "One procedure per
   ! variant"). SELAH_PROC names the procedure; SELAH_ACTIVE is the mask condition of a masked variant per layout (first
   ! = LAYOUT_V_FIRST, last = LAYOUT_V_LAST); SELAH_MASK_ARG adds the mask to the data clauses; SELAH_FILL, SELAH_WITH_DQ
   ! and SELAH_A_ZERO select the other variants of their kernels.

   ! copy: masked per layout, unmasked
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC copy_box_masked_first
#include "selah_kernel_copy_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC copy_box_masked_last
#include "selah_kernel_copy_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC copy_box_plain
#include "selah_kernel_copy_proc.INC"

   ! axpy: masked per layout, unmasked
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC axpy_box_masked_first
#include "selah_kernel_axpy_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC axpy_box_masked_last
#include "selah_kernel_axpy_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC axpy_box_plain
#include "selah_kernel_axpy_proc.INC"

   ! lincomb: masked per layout, unmasked
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC lincomb_box_masked_first
#include "selah_kernel_lincomb_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC lincomb_box_masked_last
#include "selah_kernel_lincomb_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC lincomb_box_plain
#include "selah_kernel_lincomb_proc.INC"

   ! stage_state: masked with the fill of the inactive cells per layout, masked per layout, unmasked
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_FILL
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC stage_state_box_fill_first
#include "selah_kernel_stage_state_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC stage_state_box_fill_last
#include "selah_kernel_stage_state_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_FILL
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC stage_state_box_masked_first
#include "selah_kernel_stage_state_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC stage_state_box_masked_last
#include "selah_kernel_stage_state_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC stage_state_box_plain
#include "selah_kernel_stage_state_proc.INC"

   ! step_update: with the step residual dq (masked per layout, unmasked), without it (masked per layout, unmasked)
#define SELAH_WITH_DQ
#undef SELAH_DQ_ARG
#define SELAH_DQ_ARG , dq
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC step_update_box_dq_masked_first
#include "selah_kernel_step_update_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC step_update_box_dq_masked_last
#include "selah_kernel_step_update_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC step_update_box_dq_plain
#include "selah_kernel_step_update_proc.INC"
#undef SELAH_WITH_DQ
#undef SELAH_DQ_ARG
#define SELAH_DQ_ARG
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC step_update_box_masked_first
#include "selah_kernel_step_update_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC step_update_box_masked_last
#include "selah_kernel_step_update_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC step_update_box_plain
#include "selah_kernel_step_update_proc.INC"

   ! shu_osher_2r_stage: masked per layout, unmasked
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC shu_osher_2r_stage_box_masked_first
#include "selah_kernel_shu_osher_2r_stage_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC shu_osher_2r_stage_box_masked_last
#include "selah_kernel_shu_osher_2r_stage_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC shu_osher_2r_stage_box_plain
#include "selah_kernel_shu_osher_2r_stage_proc.INC"

   ! williamson_2n_stage: a = 0 (masked per layout, unmasked), a /= 0 (masked per layout, unmasked)
#define SELAH_A_ZERO
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC williamson_2n_stage_box_a0_masked_first
#include "selah_kernel_williamson_2n_stage_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC williamson_2n_stage_box_a0_masked_last
#include "selah_kernel_williamson_2n_stage_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC williamson_2n_stage_box_a0_plain
#include "selah_kernel_williamson_2n_stage_proc.INC"
#undef SELAH_A_ZERO
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG , mask
#define SELAH_ACTIVE mask(mc,i2,i3,i4,i5) < 0._R_P
#define SELAH_PROC williamson_2n_stage_box_masked_first
#include "selah_kernel_williamson_2n_stage_proc.INC"
#undef SELAH_ACTIVE
#define SELAH_ACTIVE mask(i1,i2,i3,i4,mc) < 0._R_P
#define SELAH_PROC williamson_2n_stage_box_masked_last
#include "selah_kernel_williamson_2n_stage_proc.INC"
#undef SELAH_ACTIVE
#undef SELAH_MASK_ARG
#define SELAH_MASK_ARG
#define SELAH_PROC williamson_2n_stage_box_plain
#include "selah_kernel_williamson_2n_stage_proc.INC"
endmodule selah_kernels_box
