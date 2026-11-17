!< SELAH test: the kernels against ADAM's, bitwise on the device path, within a measured envelope on the CPU path.
program selah_test_kernels_adam
!< SELAH test: the kernels against ADAM's, bitwise on the device path, within a measured envelope on the CPU path.
!<
!< **FNL path** (layout V_LAST): ADAM's `adam_fnl_rk_kernels` verbatim (src/tests/adam_fnl_rk_kernels.F90, sha256
!< a693ca35...) and SELAH's kernels run the same SSP(5,4) step (initialize, five stages, update) and the same
!< two-register SSP(3,3) step, unmasked and masked by an immersed-boundary distance: every stage state, the advanced
!< solution and every low-storage stage are **bitwise** equal (`==`) in the cells the kernels update.
!<
!< **CPU path** (layout V_FIRST): ADAM's `rk_object` loops, transcribed (selah_adam_cpu_reference), sum the stages
!< in a different order (in place, term by term; ADAM issue #128, D6). The first stage (a copy) and the low-storage
!< stages (the same expression in both codes) are bitwise equal; the later stages, the step residual and the advanced
!< solution differ in the last bits of the summands, measured (in eps, at the O(1) scale of the summands) and bounded
!< here. Even where the summation orders coincide (stage 2: one term) the results may differ when the compiler
!< contracts to FMA, because ADAM's in-place update contracts into one rounding and SELAH's private sum rounds twice
!< (nvfortran 26.9 does this; gfortran 16.2 and ifx 2026.2 do not here). Host only, R8P only (ADAM's kernels are R8P).
use penf, only : I4P, R8P
use selah
use adam_fnl_rk_kernels
use selah_adam_cpu_reference
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: NV=5, NI=6, NJ=5, NK=4, NGC=3, NB=4, BN=3, NS=2 !< Field, pool, active blocks, solids.
real(R8P),    parameter :: DT=0.0137_R8P                                   !< Time step.
type(butcher_tableau)   :: bt                                              !< SSP(5,4).
type(shu_osher_2r_form) :: r                                               !< SSP(3,3), two-register.
real(R8P)               :: worst_ulps                                      !< Largest CPU-path difference [eps].

if (R_P /= R8P) error stop 'selah_test_kernels_adam: R8P builds only'
bt = butcher(SCHEME_SSP_54)
r = shu_osher_2r(SCHEME_SSP_33)
call fnl_butcher(.false.) ; call fnl_butcher(.true.)
call fnl_two_register(.false.) ; call fnl_two_register(.true.)
worst_ulps = 0._R8P
call cpu_butcher(.false.) ; call cpu_butcher(.true.)
call cpu_two_register(.false.) ; call cpu_two_register(.true.)
print '(A,F6.2,A)', '        CPU path: largest difference of stages 2-5, dq and the solution: ', worst_ulps, ' eps'
call check(worst_ulps <= 4._R8P, 'CPU path: stages 2-5, dq and solution within 4 eps of ADAM''s CPU loops')

if (failures > 0) error stop 'selah_test_kernels_adam: failed'
print '(A)', 'selah_test_kernels_adam: ok'

contains
   subroutine field(layout, nc, a, seed)
   !< A field of the pool with ghosts, nc components, values in [-1, 1] (a linear congruential sequence).
   integer(I4P),           intent(in)    :: layout, nc, seed
   real(R8P), allocatable, intent(inout) :: a(:,:,:,:,:)
   integer(I8P)                          :: s
   integer(I4P)                          :: i1, i2, i3, i4, i5

   if (allocated(a)) deallocate(a)
   if (layout == LAYOUT_V_LAST) then
      allocate(a(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:nc))
   else
      allocate(a(1:nc, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NB))
   endif
   s = int(seed, I8P) * 7919_I8P + 1_I8P
   do i5=lbound(a,5), ubound(a,5) ; do i4=lbound(a,4), ubound(a,4) ; do i3=lbound(a,3), ubound(a,3)
   do i2=lbound(a,2), ubound(a,2) ; do i1=lbound(a,1), ubound(a,1)
      s = modulo(s * 1103515245_I8P + 12345_I8P, 2147483648_I8P)
      a(i1,i2,i3,i4,i5) = real(s, R8P) / 1073741824._R8P - 1._R8P
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine field

   subroutine stages(layout, n, q_rk)
   !< Stage buffers of the pool, ghosts included, zeroed (ADAM allocates them zeroed).
   integer(I4P),           intent(in)    :: layout, n
   real(R8P), allocatable, intent(inout) :: q_rk(:,:,:,:,:,:)

   if (allocated(q_rk)) deallocate(q_rk)
   if (layout == LAYOUT_V_LAST) then
      allocate(q_rk(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NV, 1:n))
   else
      allocate(q_rk(1:NV, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NB, 1:n))
   endif
   q_rk = 0._R8P
   endsubroutine stages

   logical function updated(bx, phi, masked, i)
   !< The kernels update this cell: in the box and, when masked, in the fluid (phi of the last solid < 0).
   type(box),    intent(in) :: bx
   real(R8P),    intent(in) :: phi(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   logical,      intent(in) :: masked
   integer(I4P), intent(in) :: i(5)

   updated = all(i >= bx%lo) .and. all(i <= bx%hi)
   if (updated .and. masked) then
      if (bx%layout == LAYOUT_V_LAST) then
         updated = phi(i(1),i(2),i(3),i(4),NS) < 0._R8P
      else
         updated = phi(NS,i(2),i(3),i(4),i(5)) < 0._R8P
      endif
   endif
   endfunction updated

   subroutine compare(bx, phi, masked, mine, theirs, label, bitwise, ulps)
   !< Compare the updated cells: bitwise (==), or measure the largest difference in ulps of the value.
   type(box),    intent(in)              :: bx
   real(R8P),    intent(in)              :: phi(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   logical,      intent(in)              :: masked
   real(R8P),    intent(in)              :: mine(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R8P),    intent(in)              :: theirs(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   character(*), intent(in)              :: label
   logical,      intent(in)              :: bitwise
   real(R8P),    intent(inout), optional :: ulps
   integer(I4P)                          :: i1, i2, i3, i4, i5, n, differ

   n = 0 ; differ = 0
   do i5=bx%lo(5), bx%hi(5) ; do i4=bx%lo(4), bx%hi(4) ; do i3=bx%lo(3), bx%hi(3)
   do i2=bx%lo(2), bx%hi(2) ; do i1=bx%lo(1), bx%hi(1)
      if (.not. updated(bx, phi, masked, [i1,i2,i3,i4,i5])) cycle
      n = n + 1
      if (mine(i1,i2,i3,i4,i5) /= theirs(i1,i2,i3,i4,i5)) then
         differ = differ + 1
         ! in units of eps at the scale of the summands (|q| <= 1, small dt terms): an ulp of the result itself would
         ! inflate a last-bit difference of the O(1) summands wherever the result cancels toward 0
         if (present(ulps)) ulps = max(ulps, abs(mine(i1,i2,i3,i4,i5) - theirs(i1,i2,i3,i4,i5)) / epsilon(1._R8P))
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   if (bitwise) call check(n > 0 .and. differ == 0, label//': bitwise equal')
   endsubroutine compare

   subroutine fnl_butcher(masked)
   !< SSP(5,4) step, FNL layout: ADAM's kernels against SELAH's, stage by stage.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), qa(:,:,:,:,:), qs(:,:,:,:,:), u(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:)
   real(R8P), allocatable :: q_rk(:,:,:,:,:,:), regs(:,:,:,:,:,:)
   type(box)              :: bx
   character(8)           :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, BN)
   call field(LAYOUT_V_LAST, NV, q, 1_I4P) ; call field(LAYOUT_V_LAST, NS, phi, 2_I4P)
   qa = q ; qs = q ; u = q
   call stages(LAYOUT_V_LAST, bt%stages, q_rk) ; call stages(LAYOUT_V_LAST, bt%stages, regs)
   call rk_initialize_stages_dev(ni=NI, nj=NJ, nk=NK, ngc=NGC, nv=NV, blocks_number=BN, q_gpu=q, q_rk_gpu=q_rk)
   do s=1, bt%stages
      if (masked) then
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, BN, s, DT, bt%a, phi_gpu=phi, q_rk_gpu=q_rk)
         call selah_stage_state(bx, s-1, DT, bt%a(s,:), q, regs, u, mask=phi)
      else
         call rk_compute_stage_dev(NI, NJ, NK, NGC, NV, BN, s, DT, bt%a, q_rk_gpu=q_rk)
         call selah_stage_state(bx, s-1, DT, bt%a(s,:), q, regs, u)
      endif
      call compare(bx, phi, masked, u, q_rk(:,:,:,:,:,s), 'FNL ssp-54 '//trim(tag)//' stage '//achar(48+s)//' state', .true.)
      call field(LAYOUT_V_LAST, NV, res, 10_I4P + s) ! the residual of the stage (any field will do)
      if (masked) then
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, BN, s, phi_gpu=phi, q_gpu=res, q_rk_gpu=q_rk)
         call selah_copy(bx, res, regs(:,:,:,:,:,s), mask=phi)
      else
         call rk_assign_stage_dev(NI, NJ, NK, NGC, NV, BN, s, q_gpu=res, q_rk_gpu=q_rk)
         call selah_copy(bx, res, regs(:,:,:,:,:,s))
      endif
   enddo
   if (masked) then
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, bt%b, phi_gpu=phi, q_rk_gpu=q_rk, q_gpu=qa)
      call selah_step_update(bx, bt%stages, DT, bt%b, regs, qs, mask=phi)
   else
      call rk_update_q_dev(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, bt%b, q_rk_gpu=q_rk, q_gpu=qa)
      call selah_step_update(bx, bt%stages, DT, bt%b, regs, qs)
   endif
   call compare(bx, phi, masked, qs, qa, 'FNL ssp-54 '//trim(tag)//' advanced solution', .true.)
   endsubroutine fnl_butcher

   subroutine fnl_two_register(masked)
   !< SSP(3,3) two-register step, FNL layout: ADAM's compute_stage_ls against SELAH's shu_osher_2r_stage.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), qa(:,:,:,:,:), qs(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:)
   real(R8P), allocatable :: q_rk(:,:,:,:,:,:)
   type(box)              :: bx
   character(8)           :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, BN)
   call field(LAYOUT_V_LAST, NV, q, 3_I4P) ; call field(LAYOUT_V_LAST, NS, phi, 4_I4P)
   qa = q ; qs = q
   call stages(LAYOUT_V_LAST, 1_I4P, q_rk)
   call rk_initialize_stages_dev(ni=NI, nj=NJ, nk=NK, ngc=NGC, nv=NV, blocks_number=BN, q_gpu=q, q_rk_gpu=q_rk)
   do s=1, r%stages
      call field(LAYOUT_V_LAST, NV, res, 20_I4P + s)
      if (masked) then
         call rk_compute_stage_ls_dev(NI, NJ, NK, NGC, NV, BN, DT, r%a(s), r%b(s), r%c(s), phi_gpu=phi, &
                                      q_n_gpu=q_rk(:,:,:,:,:,1), dq_gpu=res, q_rk_gpu=qa)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q, res, qs, mask=phi)
      else
         call rk_compute_stage_ls_dev(NI, NJ, NK, NGC, NV, BN, DT, r%a(s), r%b(s), r%c(s), &
                                      q_n_gpu=q_rk(:,:,:,:,:,1), dq_gpu=res, q_rk_gpu=qa)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q, res, qs)
      endif
      call compare(bx, phi, masked, qs, qa, 'FNL ssp-33 two-register '//trim(tag)//' stage '//achar(48+s), .true.)
   enddo
   endsubroutine fnl_two_register

   subroutine cpu_butcher(masked)
   !< SSP(5,4) step, CPU layout: ADAM's transcribed CPU loops against SELAH's kernels.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), qa(:,:,:,:,:), qs(:,:,:,:,:), u(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:)
   real(R8P), allocatable :: dqa(:,:,:,:,:), dqs(:,:,:,:,:), q_rk(:,:,:,:,:,:), regs(:,:,:,:,:,:)
   type(box)              :: bx
   character(8)           :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_FIRST, NV, NI, NJ, NK, NGC, BN)
   call field(LAYOUT_V_FIRST, NV, q, 5_I4P) ; call field(LAYOUT_V_FIRST, NS, phi, 6_I4P)
   qa = q ; qs = q ; u = q ; dqa = 0._R8P * q ; dqs = dqa
   call stages(LAYOUT_V_FIRST, bt%stages, q_rk) ; call stages(LAYOUT_V_FIRST, bt%stages, regs)
   call adam_cpu_initialize_stages(NI, NJ, NK, NGC, NV, BN, q, q_rk)
   do s=1, bt%stages
      if (masked) then
         call adam_cpu_compute_stage(NI, NJ, NK, NGC, NV, BN, s, DT, bt%a, q_rk, phi)
         call selah_stage_state(bx, s-1, DT, bt%a(s,:), q, regs, u, mask=phi)
      else
         call adam_cpu_compute_stage(NI, NJ, NK, NGC, NV, BN, s, DT, bt%a, q_rk)
         call selah_stage_state(bx, s-1, DT, bt%a(s,:), q, regs, u)
      endif
      if (s == 1) then ! a copy, no arithmetic: bitwise on every compiler
         call compare(bx, phi, masked, u, q_rk(:,:,:,:,:,s), 'CPU ssp-54 '//trim(tag)//' stage '//achar(48+s)//' state', &
                      .true.)
      else
         call compare(bx, phi, masked, u, q_rk(:,:,:,:,:,s), '', .false., worst_ulps)
      endif
      call field(LAYOUT_V_FIRST, NV, res, 30_I4P + s)
      if (masked) then
         call adam_cpu_assign_stage(NI, NJ, NK, NGC, NV, BN, s, res, q_rk, phi)
         call selah_copy(bx, res, regs(:,:,:,:,:,s), mask=phi)
      else
         call adam_cpu_assign_stage(NI, NJ, NK, NGC, NV, BN, s, res, q_rk)
         call selah_copy(bx, res, regs(:,:,:,:,:,s))
      endif
   enddo
   if (masked) then
      call adam_cpu_update_q(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, bt%b, q_rk, qa, phi=phi)
      call selah_step_update(bx, bt%stages, DT, bt%b, regs, qs, mask=phi)
   else
      call adam_cpu_update_q(NI, NJ, NK, NGC, NV, BN, bt%stages, DT, bt%b, q_rk, qa, dq=dqa)
      call selah_step_update(bx, bt%stages, DT, bt%b, regs, qs, dq=dqs)
      call compare(bx, phi, masked, dqs, dqa, '', .false., worst_ulps)
   endif
   call compare(bx, phi, masked, qs, qa, '', .false., worst_ulps)
   endsubroutine cpu_butcher

   subroutine cpu_two_register(masked)
   !< SSP(3,3) two-register step, CPU layout: ADAM's transcribed compute_stage_ls against shu_osher_2r_stage.
   logical, intent(in)    :: masked
   real(R8P), allocatable :: q(:,:,:,:,:), qa(:,:,:,:,:), qs(:,:,:,:,:), phi(:,:,:,:,:), res(:,:,:,:,:)
   real(R8P), allocatable :: q_rk(:,:,:,:,:,:)
   type(box)              :: bx
   character(8)           :: tag
   integer(I4P)           :: s

   tag = merge('masked  ', 'unmasked', masked)
   bx = box_interior(LAYOUT_V_FIRST, NV, NI, NJ, NK, NGC, BN)
   call field(LAYOUT_V_FIRST, NV, q, 7_I4P) ; call field(LAYOUT_V_FIRST, NS, phi, 8_I4P)
   qa = q ; qs = q
   call stages(LAYOUT_V_FIRST, 1_I4P, q_rk)
   call adam_cpu_initialize_stages(NI, NJ, NK, NGC, NV, BN, q, q_rk)
   do s=1, r%stages
      call field(LAYOUT_V_FIRST, NV, res, 40_I4P + s)
      if (masked) then
         call adam_cpu_compute_stage_ls(NI, NJ, NK, NGC, NV, BN, s, DT, r%a, r%b, r%c, q_rk, res, qa, phi)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q, res, qs, mask=phi)
      else
         call adam_cpu_compute_stage_ls(NI, NJ, NK, NGC, NV, BN, s, DT, r%a, r%b, r%c, q_rk, res, qa)
         call selah_shu_osher_2r_stage(bx, DT, r%a(s), r%b(s), r%c(s), q, res, qs)
      endif
      call compare(bx, phi, masked, qs, qa, 'CPU ssp-33 two-register '//trim(tag)//' stage '//achar(48+s), .true.)
   enddo
   endsubroutine cpu_two_register
endprogram selah_test_kernels_adam
