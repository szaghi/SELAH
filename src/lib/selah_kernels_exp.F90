!< SELAH, L1: the exponential box kernels, a Butcher stage and step with an exactly integrated linear decay.

#include "fundal.H"

! the dummy arrays carry the lower bounds of the allocation (ghosts included): no rebasing to 1
#define SELAH_BOUNDS bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):

module selah_kernels_exp
!< SELAH, L1: the exponential box kernels, a Butcher stage and step with an exactly integrated linear decay.
!<
!< For a state whose equation is psi' = g - lambda psi, lambda >= 0 per element and constant in time (the auxiliary
!< fields of a perfectly matched layer), the caller's residual is the source g only and the decay is integrated
!< exactly: with z = -lambda dt, the state of stage i is e^(c_i z) psi^n + phi_1(c_i z) dt sum_j a_ij g_j and the step
!< e^z psi^n + phi_1(z) dt sum_j b_j g_j, phi_1(z) = (e^z - 1)/z (exponential Euler on the Runge–Kutta-averaged
!< source: `notes/companions-pushers-pml.md` §5, variant F). Exact for a source frozen over the step, any lambda dt;
!< second order with a varying source; where lambda = 0, bitwise the explicit kernels (selah_kernels_box). No mask.
!< phi_1 is evaluated without cancellation: e^(z/2) sinh(z/2)/(z/2), and (e^z - 1)/z where e^z is negligible.
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_box, only : box
implicit none
private
public :: stage_state_exp_box
public :: step_update_exp_box

contains
   subroutine stage_state_exp_box(bx, n, dt, coef, c, q, regs, lam, u)
   !< u = e^(c z) q + phi_1(c z) sum_{k=1..n} dt coef(k) regs(k), z = -lam dt, over the box (u = q + the sum where
   !< lam = 0, as stage_state_box).
   type(box),    intent(in)    :: bx                    !< Box.
   integer(I4P), intent(in)    :: n                     !< Registers summed.
   real(R_P),    intent(in)    :: dt                    !< Time step.
   real(R_P),    intent(in)    :: coef(1:)              !< Coefficients [n or more] (row of A).
   real(R_P),    intent(in)    :: c                     !< The stage's abscissa.
   real(R_P),    intent(in)    :: q(SELAH_BOUNDS)       !< Solution.
   real(R_P),    intent(in)    :: regs(SELAH_BOUNDS,1:) !< Source registers [..., n or more].
   real(R_P),    intent(in)    :: lam(SELAH_BOUNDS)     !< Decay rate, >= 0.
   real(R_P),    intent(inout) :: u(SELAH_BOUNDS)       !< Stage state.
   real(R_P)                   :: d, z, x, ph           !< Private sum, c z, z/2, phi_1.
   integer(I4P)                :: k                     !< Counter.
#include "selah_kernel_decl.INC"

   lo1 = bx%lo(1) ; lo2 = bx%lo(2) ; lo3 = bx%lo(3) ; lo4 = bx%lo(4) ; lo5 = bx%lo(5)
   hi1 = bx%hi(1) ; hi2 = bx%hi(2) ; hi3 = bx%hi(3) ; hi4 = bx%hi(4) ; hi5 = bx%hi(5)
   !$acc parallel loop independent collapse(5) DEVICEVAR(coef, q, regs, lam, u) private(d, z, x, ph)
   !$omp OMPLOOP collapse(5) DEVICEPTR(coef, q, regs, lam, u) private(d, z, x, ph)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
                  d = 0._R_P
                  !$acc loop seq
                  do k=1, n
                     d = d + dt * coef(k) * regs(i1,i2,i3,i4,i5,k)
                  enddo
                  z = -c * lam(i1,i2,i3,i4,i5) * dt
                  if (z == 0._R_P) then
                     u(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5) + d
                  else
                     x = 0.5_R_P * z
                     if (x < -20._R_P) then
                        ph = (exp(z) - 1._R_P) / z
                     else
                        ph = exp(x) * sinh(x) / x
                     endif
                     u(i1,i2,i3,i4,i5) = exp(z) * q(i1,i2,i3,i4,i5) + ph * d
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine stage_state_exp_box

   subroutine step_update_exp_box(bx, n, dt, coef, regs, lam, q)
   !< q = e^z q + phi_1(z) dt sum_{k=1..n} coef(k) regs(k), z = -lam dt, over the box (q = q + dt times the sum where
   !< lam = 0, as step_update_box).
   type(box),    intent(in)    :: bx                    !< Box.
   integer(I4P), intent(in)    :: n                     !< Registers summed.
   real(R_P),    intent(in)    :: dt                    !< Time step.
   real(R_P),    intent(in)    :: coef(1:)              !< Weights [n or more] (b).
   real(R_P),    intent(in)    :: regs(SELAH_BOUNDS,1:) !< Source registers [..., n or more].
   real(R_P),    intent(in)    :: lam(SELAH_BOUNDS)     !< Decay rate, >= 0.
   real(R_P),    intent(inout) :: q(SELAH_BOUNDS)       !< Solution, advanced.
   real(R_P)                   :: d, z, x, ph           !< Private sum, z, z/2, phi_1.
   integer(I4P)                :: k                     !< Counter.
#include "selah_kernel_decl.INC"

   lo1 = bx%lo(1) ; lo2 = bx%lo(2) ; lo3 = bx%lo(3) ; lo4 = bx%lo(4) ; lo5 = bx%lo(5)
   hi1 = bx%hi(1) ; hi2 = bx%hi(2) ; hi3 = bx%hi(3) ; hi4 = bx%hi(4) ; hi5 = bx%hi(5)
   !$acc parallel loop independent collapse(5) DEVICEVAR(coef, regs, lam, q) private(d, z, x, ph)
   !$omp OMPLOOP collapse(5) DEVICEPTR(coef, regs, lam, q) private(d, z, x, ph)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
                  d = 0._R_P
                  !$acc loop seq
                  do k=1, n
                     d = d + coef(k) * regs(i1,i2,i3,i4,i5,k)
                  enddo
                  z = -lam(i1,i2,i3,i4,i5) * dt
                  if (z == 0._R_P) then
                     q(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5) + dt * d
                  else
                     x = 0.5_R_P * z
                     if (x < -20._R_P) then
                        ph = (exp(z) - 1._R_P) / z
                     else
                        ph = exp(x) * sinh(x) / x
                     endif
                     q(i1,i2,i3,i4,i5) = exp(z) * q(i1,i2,i3,i4,i5) + ph * (dt * d)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine step_update_exp_box
endmodule selah_kernels_exp
