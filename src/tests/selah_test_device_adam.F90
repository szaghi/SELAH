!< SELAH test: the kernels against ADAM's FNL kernels, both on device arrays, bitwise.

#include "fundal.H"

program selah_test_device_adam
!< SELAH test: the kernels against ADAM's FNL kernels, both on device arrays, bitwise.
!<
!< The device twin of the FNL half of selah_test_kernels_adam: ADAM's adam_fnl_rk_kernels (verbatim, the oracle) and
!< SELAH's kernels run on the same device fields (FUNDAL; layout V_LAST, ADAM's device layout) an SSP(5,4) step and a
!< two-register SSP(3,3) step, unmasked and masked; every stage state, the advanced solution and every two-register
!< stage, copied back, are bitwise equal (==). Coefficients are device resident, as the kernels require: ADAM's whole
!< alph matrix and beta, SELAH's rows. R8P only (ADAM's kernels are R8P).
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device, dev_assign_to_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use penf, only : I4P, I8P, R8P
use selah
use adam_fnl_rk_kernels
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: NV=5, NI=6, NJ=5, NK=4, NGC=3, NB=4, BN=3, NS=2
real(R8P),    parameter :: DT=0.0137_R8P
type(butcher_tableau)   :: bt
type(shu_osher_2r_form) :: r

if (R_P /= R8P) error stop 'selah_test_device_adam: R8P builds only'
#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build: both kernel sets run on the device'
#else
print '(A)', 'host build: FUNDAL host memory'
#endif
bt = butcher(SCHEME_SSP_54)
r = shu_osher_2r(SCHEME_SSP_33)
call butcher_step(.false.) ; call butcher_step(.true.)
call two_register_step(.false.) ; call two_register_step(.true.)

if (failures > 0) error stop 'selah_test_device_adam: failed'
print '(A)', 'selah_test_device_adam: ok'

contains
   subroutine field(nc, a, seed)
   !< A host field of the pool in layout V_LAST (blocks first), nc components, values in [-1, 1].
   integer(I4P),           intent(in)    :: nc, seed
   real(R8P), allocatable, intent(inout) :: a(:,:,:,:,:)
   integer(I8P)                          :: s
   integer(I4P)                          :: i1, i2, i3, i4, i5

   if (allocated(a)) deallocate(a)
   allocate(a(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:nc))
   s = int(seed, I8P) * 7919_I8P + 1_I8P
   do i5=lbound(a,5), ubound(a,5) ; do i4=lbound(a,4), ubound(a,4) ; do i3=lbound(a,3), ubound(a,3)
   do i2=lbound(a,2), ubound(a,2) ; do i1=lbound(a,1), ubound(a,1)
      s = modulo(s * 1103515245_I8P + 12345_I8P, 2147483648_I8P)
      a(i1,i2,i3,i4,i5) = real(s, R8P) / 1073741824._R8P - 1._R8P
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine field

   subroutine dalloc5(a, a_d)
   !< A device array with the bounds of a host one, filled from it.
   real(R8P),          intent(in)    :: a(:,:,:,:,:)
   real(R8P), pointer, intent(inout) :: a_d(:,:,:,:,:)
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=[NB, NI+NGC, NJ+NGC, NK+NGC, size(a,5,I4P)], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P], ierr=ierr)
   if (ierr /= 0) error stop 'selah_test_device_adam: dev_alloc failed'
   call dev_memcpy_to_device(dst=a_d, src=a)
   endsubroutine dalloc5

   subroutine dalloc6(n, a_d)
   !< Zeroed device stage buffers of the pool, ghosts included.
   integer(I4P),       intent(in)    :: n
   real(R8P), pointer, intent(inout) :: a_d(:,:,:,:,:,:)
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=[NB, NI+NGC, NJ+NGC, NK+NGC, NV, n], &
                  lbounds=[1_I4P, 1_I4P-NGC, 1_I4P-NGC, 1_I4P-NGC, 1_I4P, 1_I4P], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0) error stop 'selah_test_device_adam: dev_alloc failed'
   endsubroutine dalloc6

   subroutine compare(phi, masked, mine_d, theirs_d, label)
   !< Copy both device results back; the cells the kernels update (interior of the active blocks, fluid if masked)
   !< must be bitwise equal.
   real(R8P),          intent(in) :: phi(1:,1-NGC:,1-NGC:,1-NGC:,1:)
   logical,            intent(in) :: masked
   real(R8P), target,  intent(in) :: mine_d(:,:,:,:,:), theirs_d(:,:,:,:,:) !< Device results (any section).
   character(*),       intent(in) :: label
   real(R8P), allocatable         :: mine(:,:,:,:,:), theirs(:,:,:,:,:)
   integer(I4P)                   :: b, i, j, k, v, n, differ

   allocate(mine(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NV), theirs(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, &
            1-NGC:NK+NGC, 1:NV))
   call dev_memcpy_from_device(dst=mine, src=mine_d) ; call dev_memcpy_from_device(dst=theirs, src=theirs_d)
   n = 0 ; differ = 0
   do v=1, NV ; do k=1, NK ; do j=1, NJ ; do i=1, NI ; do b=1, BN
      if (masked) then
         if (phi(b,i,j,k,NS) >= 0._R8P) cycle
      endif
      n = n + 1
      if (mine(b,i,j,k,v) /= theirs(b,i,j,k,v)) differ = differ + 1
   enddo ; enddo ; enddo ; enddo ; enddo
   call check(n > 0 .and. differ == 0, label//': bitwise equal')
   endsubroutine compare

   subroutine butcher_step(masked)
   !< SSP(5,4): ADAM's initialize/compute/assign/update against SELAH's stage_state/copy/step_update, on the device.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:), row(:)
   real(R8P), pointer     :: q_d(:,:,:,:,:)=>null(), qa_d(:,:,:,:,:)=>null(), qs_d(:,:,:,:,:)=>null()
   real(R8P), pointer     :: u_d(:,:,:,:,:)=>null(), phi_d(:,:,:,:,:)=>null(), res_d(:,:,:,:,:)=>null()
   real(R8P), pointer     :: q_rk_d(:,:,:,:,:,:)=>null(), regs_d(:,:,:,:,:,:)=>null()
   real(R8P), pointer     :: alph_d(:,:)=>null(), beta_d(:)=>null(), row_d(:)=>null()
   type(box)              :: bx
   character(:), allocatable :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, BN)
   call field(NV, q, 1_I4P) ; call field(NS, phi, 2_I4P)
   call dalloc5(q, q_d) ; call dalloc5(q, qa_d) ; call dalloc5(q, qs_d) ; call dalloc5(q, u_d) ; call dalloc5(phi, phi_d)
   call dalloc6(bt%stages, q_rk_d) ; call dalloc6(bt%stages, regs_d)
   call dev_assign_to_device(dst=alph_d, src=bt%a) ; call dev_assign_to_device(dst=beta_d, src=bt%b)
   call rk_initialize_stages_dev(ni=NI, nj=NJ, nk=NK, ngc=NGC, nv=NV, blocks_number=BN, q_gpu=q_d, q_rk_gpu=q_rk_d)
   do s=1, bt%stages
      row = bt%a(s,:) ! a contiguous copy of the row of A
      call dev_assign_to_device(dst=row_d, src=row)
      if (masked) then
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, BN, s, DT, alph_d, phi_gpu=phi_d, q_rk_gpu=q_rk_d)
         call selah_stage_state(bx, s-1, DT, row_d, q_d, regs_d, u_d, mask=phi_d)
      else
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, BN, s, DT, alph_d, q_rk_gpu=q_rk_d)
         call selah_stage_state(bx, s-1, DT, row_d, q_d, regs_d, u_d)
      endif
      call compare(phi, masked, u_d, q_rk_d(:,:,:,:,:,s), 'device ssp-54 '//trim(tag)//' stage '//achar(48+s)//' state')
      call field(NV, res, 10_I4P + s)
      call dalloc5(res, res_d)
      if (masked) then
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, BN, s, phi_gpu=phi_d, q_gpu=res_d, q_rk_gpu=q_rk_d)
         call selah_copy(bx, res_d, regs_d(:,:,:,:,:,s), mask=phi_d)
      else
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, BN, s, q_gpu=res_d, q_rk_gpu=q_rk_d)
         call selah_copy(bx, res_d, regs_d(:,:,:,:,:,s))
      endif
      call dev_free(res_d) ; call dev_free(row_d)
   enddo
   if (masked) then
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, beta_d, phi_gpu=phi_d, q_rk_gpu=q_rk_d, q_gpu=qa_d)
      call selah_step_update(bx, bt%stages, DT, beta_d, regs_d, qs_d, mask=phi_d)
   else
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, beta_d, q_rk_gpu=q_rk_d, q_gpu=qa_d)
      call selah_step_update(bx, bt%stages, DT, beta_d, regs_d, qs_d)
   endif
   call compare(phi, masked, qs_d, qa_d, 'device ssp-54 '//trim(tag)//' advanced solution')
   call dev_free(q_d) ; call dev_free(qa_d) ; call dev_free(qs_d) ; call dev_free(u_d) ; call dev_free(phi_d)
   call dev_free(q_rk_d) ; call dev_free(regs_d) ; call dev_free(alph_d) ; call dev_free(beta_d)
   endsubroutine butcher_step

   subroutine two_register_step(masked)
   !< SSP(3,3) two-register: ADAM's compute_stage_ls against SELAH's shu_osher_2r_stage, on the device.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:)
   real(R8P), pointer     :: q_d(:,:,:,:,:)=>null(), qa_d(:,:,:,:,:)=>null(), qs_d(:,:,:,:,:)=>null()
   real(R8P), pointer     :: phi_d(:,:,:,:,:)=>null(), res_d(:,:,:,:,:)=>null(), q_rk_d(:,:,:,:,:,:)=>null()
   type(box)              :: bx
   character(:), allocatable :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, BN)
   call field(NV, q, 3_I4P) ; call field(NS, phi, 4_I4P)
   call dalloc5(q, q_d) ; call dalloc5(q, qa_d) ; call dalloc5(q, qs_d) ; call dalloc5(phi, phi_d)
   call dalloc6(1_I4P, q_rk_d)
   call rk_initialize_stages_dev(ni=NI, nj=NJ, nk=NK, ngc=NGC, nv=NV, blocks_number=BN, q_gpu=q_d, q_rk_gpu=q_rk_d)
   do s=1, r%stages
      call field(NV, res, 20_I4P + s)
      call dalloc5(res, res_d)
      if (masked) then
         call rk_compute_stage_ls_dev(NI, NJ, NK, NGC, NV, BN, DT, r%a(s), r%b(s), r%c(s), phi_gpu=phi_d, &
                                      q_n_gpu=q_rk_d(:,:,:,:,:,1), dq_gpu=res_d, q_rk_gpu=qa_d)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q_d, res_d, qs_d, mask=phi_d)
      else
         call rk_compute_stage_ls_dev(NI, NJ, NK, NGC, NV, BN, DT, r%a(s), r%b(s), r%c(s), &
                                      q_n_gpu=q_rk_d(:,:,:,:,:,1), dq_gpu=res_d, q_rk_gpu=qa_d)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q_d, res_d, qs_d)
      endif
      call compare(phi, masked, qs_d, qa_d, 'device ssp-33 two-register '//trim(tag)//' stage '//achar(48+s))
      call dev_free(res_d)
   enddo
   call dev_free(q_d) ; call dev_free(qa_d) ; call dev_free(qs_d) ; call dev_free(phi_d) ; call dev_free(q_rk_d)
   endsubroutine two_register_step
endprogram selah_test_device_adam
