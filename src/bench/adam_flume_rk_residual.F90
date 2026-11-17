!< ADAM's FLUME FNL step residual kernel, the oracle of the M5 benchmark.
!<
!< The subroutine is a verbatim copy of `compute_rk_ssp_residual_dev`, ADAM `src/app/flume/fnl/adam_flume_fnl_kernels.F90`
!< lines 578-614 (file sha256 974e532aa828a5a92e452870b2c31a0ea55cc39becb570f08e57bd731b288001, read 2026-11-17).
!< The module around it is SELAH's. FLUME FNL calls it after `rk_update_q_dev`, which returns no step residual
!< (ADAM#128 D4): it rereads every stage register.

#include "fundal.H"

module adam_flume_rk_residual
!< ADAM's FLUME FNL step residual kernel, the oracle of the M5 benchmark.
use :: fundal
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_rk_ssp_residual_dev

contains
   subroutine compute_rk_ssp_residual_dev(ni, nj, nk, ngc, nv, blocks_number, nrk, beta_gpu, q_rk_gpu, dq_gpu)
   !< Compute the residual of a strong stability preserving step, `dq = sum_s beta_s dq_s`, from the stored stages.
   !<
   !< Device twin of the residual that the host `rk_object%update_q` returns: the residuals history of the two backends
   !< then reports the same quantity (the effective step residual, not the last stage one).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                      !< Variables number.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: nrk                                     !< Runge-Kutta stages number.
   real(R8P),    intent(in)    :: beta_gpu(1:)                            !< Runge-Kutta beta coefficients.
   real(R8P),    intent(in)    :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Stored stage residuals.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Step residual.
   real(R8P)                   :: dq                                      !< Private step residual of one cell.
   integer(I4P)                :: b, i, j, k, v, s                        !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(beta_gpu,q_rk_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,nv,blocks_number,nrk) private(dq)
   !$omp OMPLOOP collapse(4) DEVICEPTR(beta_gpu,q_rk_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,nv,blocks_number,nrk) private(dq)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, nv
         dq = 0._R8P
         !$acc loop seq
         do s=1, nrk
            dq = dq + beta_gpu(s) * q_rk_gpu(b,i,j,k,v,s)
         enddo
         dq_gpu(b,i,j,k,v) = dq
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_rk_ssp_residual_dev
endmodule adam_flume_rk_residual
