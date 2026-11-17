!< SELAH, L1: the flat kernels, the stage operations over a whole contiguous rank-1 array.

#include "fundal.H"

module selah_kernels_flat
!< SELAH, L1: the flat kernels, the stage operations over a whole contiguous rank-1 array.
!<
!< For state that has no ghost cells and no mask: particles (position, velocity), perfectly-matched-layer slabs, any
!< array the caller flattens. Every element is updated; the operations and their summation orders are those of the
!< box kernels (selah_kernels_box). In a device build every array argument must be device resident.
use penf, only : I4P
use selah_kind_parameters, only : R_P
implicit none
private
public :: copy_flat
public :: axpy_flat
public :: lincomb_flat
public :: stage_state_flat
public :: step_update_flat
public :: shu_osher_2r_stage_flat
public :: williamson_2n_stage_flat

contains
   subroutine copy_flat(x, y)
   !< y = x.
   real(R_P), intent(in)    :: x(1:) !< Source.
   real(R_P), intent(inout) :: y(1:) !< Destination.
   integer(I4P)             :: i     !< Counter.

   !$acc parallel loop independent DEVICEVAR(x, y)
   !$omp OMPLOOP DEVICEPTR(x, y)
   do i=1, size(y, kind=I4P)
      y(i) = x(i)
   enddo
   endsubroutine copy_flat

   subroutine axpy_flat(alpha, x, y)
   !< y = y + alpha x.
   real(R_P), intent(in)    :: alpha !< Coefficient.
   real(R_P), intent(in)    :: x(1:) !< Increment.
   real(R_P), intent(inout) :: y(1:) !< Updated array.
   integer(I4P)             :: i     !< Counter.

   !$acc parallel loop independent DEVICEVAR(x, y)
   !$omp OMPLOOP DEVICEPTR(x, y)
   do i=1, size(y, kind=I4P)
      y(i) = y(i) + alpha * x(i)
   enddo
   endsubroutine axpy_flat

   subroutine lincomb_flat(alpha, x, beta, y)
   !< y = alpha x + beta y.
   real(R_P), intent(in)    :: alpha !< Coefficient of x.
   real(R_P), intent(in)    :: x(1:) !< First array.
   real(R_P), intent(in)    :: beta  !< Coefficient of y.
   real(R_P), intent(inout) :: y(1:) !< Second array, overwritten.
   integer(I4P)             :: i     !< Counter.

   !$acc parallel loop independent DEVICEVAR(x, y)
   !$omp OMPLOOP DEVICEPTR(x, y)
   do i=1, size(y, kind=I4P)
      y(i) = alpha * x(i) + beta * y(i)
   enddo
   endsubroutine lincomb_flat

   subroutine stage_state_flat(n, dt, coef, q, regs, u)
   !< u = q + sum_{k=1..n} dt coef(k) regs(:,k).
   integer(I4P), intent(in)    :: n          !< Registers summed.
   real(R_P),    intent(in)    :: dt         !< Time step.
   real(R_P),    intent(in)    :: coef(1:)   !< Coefficients [n or more].
   real(R_P),    intent(in)    :: q(1:)      !< Solution.
   real(R_P),    intent(in)    :: regs(1:,1:) !< Residual registers [size(q), n or more].
   real(R_P),    intent(inout) :: u(1:)      !< Stage state.
   real(R_P)                   :: d          !< Private sum.
   integer(I4P)                :: i, k       !< Counters.

   !$acc parallel loop independent DEVICEVAR(coef, q, regs, u) private(d)
   !$omp OMPLOOP DEVICEPTR(coef, q, regs, u) private(d)
   do i=1, size(u, kind=I4P)
      d = 0._R_P
      !$acc loop seq
      do k=1, n
         d = d + dt * coef(k) * regs(i,k)
      enddo
      u(i) = q(i) + d
   enddo
   endsubroutine stage_state_flat

   subroutine step_update_flat(n, dt, coef, regs, q, dq)
   !< q = q + dt sum_{k=1..n} coef(k) regs(:,k), fused with dq = sum_k coef(k) regs(:,k) when dq is present.
   integer(I4P), intent(in)              :: n           !< Registers summed.
   real(R_P),    intent(in)              :: dt          !< Time step.
   real(R_P),    intent(in)              :: coef(1:)    !< Weights [n or more].
   real(R_P),    intent(in)              :: regs(1:,1:) !< Residual registers [size(q), n or more].
   real(R_P),    intent(inout)           :: q(1:)       !< Solution, advanced.
   real(R_P),    intent(inout), optional :: dq(1:)      !< Step residual.
   real(R_P)                             :: d           !< Private sum.
   integer(I4P)                          :: i, k        !< Counters.

   if (present(dq)) then
      !$acc parallel loop independent DEVICEVAR(coef, regs, q, dq) private(d)
      !$omp OMPLOOP DEVICEPTR(coef, regs, q, dq) private(d)
      do i=1, size(q, kind=I4P)
         d = 0._R_P
         !$acc loop seq
         do k=1, n
            d = d + coef(k) * regs(i,k)
         enddo
         q(i) = q(i) + dt * d
         dq(i) = d
      enddo
   else
      !$acc parallel loop independent DEVICEVAR(coef, regs, q) private(d)
      !$omp OMPLOOP DEVICEPTR(coef, regs, q) private(d)
      do i=1, size(q, kind=I4P)
         d = 0._R_P
         !$acc loop seq
         do k=1, n
            d = d + coef(k) * regs(i,k)
         enddo
         q(i) = q(i) + dt * d
      enddo
   endif
   endsubroutine step_update_flat

   subroutine shu_osher_2r_stage_flat(dt, a, b, c, q0, r, q)
   !< q = a q0 + b q + dt c r: a stage of a two-register Shu–Osher form.
   real(R_P), intent(in)    :: dt      !< Time step.
   real(R_P), intent(in)    :: a, b, c !< Stage coefficients.
   real(R_P), intent(in)    :: q0(1:)  !< Solution at the start of the step.
   real(R_P), intent(in)    :: r(1:)   !< Residual of the current stage.
   real(R_P), intent(inout) :: q(1:)   !< Current stage, advanced.
   integer(I4P)             :: i       !< Counter.

   !$acc parallel loop independent DEVICEVAR(q0, r, q)
   !$omp OMPLOOP DEVICEPTR(q0, r, q)
   do i=1, size(q, kind=I4P)
      q(i) = a * q0(i) + b * q(i) + dt * c * r(i)
   enddo
   endsubroutine shu_osher_2r_stage_flat

   subroutine williamson_2n_stage_flat(dt, a, b, r, dqr, q)
   !< dqr = a dqr + dt r; q = q + b dqr: a stage of a Williamson 2N form.
   real(R_P), intent(in)    :: dt     !< Time step.
   real(R_P), intent(in)    :: a, b   !< Stage coefficients A_i, B_i.
   real(R_P), intent(in)    :: r(1:)  !< Residual of the current stage.
   real(R_P), intent(inout) :: dqr(1:) !< The Delta Q register.
   real(R_P), intent(inout) :: q(1:)  !< Current stage, advanced.
   integer(I4P)             :: i      !< Counter.

   !$acc parallel loop independent DEVICEVAR(r, dqr, q)
   !$omp OMPLOOP DEVICEPTR(r, dqr, q)
   do i=1, size(q, kind=I4P)
      dqr(i) = a * dqr(i) + dt * r(i)
      q(i) = q(i) + b * dqr(i)
   enddo
   endsubroutine williamson_2n_stage_flat
endmodule selah_kernels_flat
