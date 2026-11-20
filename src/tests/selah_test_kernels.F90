!< SELAH test: the box and flat kernels, value by value and cell by cell.
program selah_test_kernels
!< SELAH test: the box and flat kernels, value by value and cell by cell.
!<
!< Fields of 3 components, 4x3x2 interior cells with 2 ghosts, in a pool of 6 blocks, in both layouts; every kernel
!< runs on the interior of the first 3 and then the first 5 blocks (the active prefix grows), unmasked and masked
!< (the last mask component and the first, both named), and on a component range (a partition). The default mask
!< component (mask_comp absent) is the last one, bitwise. The 2N stage with a = 0 does not read its register (NaN
!< in, finite out). Each cell of the box must
!< hold the reference value (the same expression in a plain loop, within 16 eps: compilers may contract to FMA
!< differently), and **every other cell must be bitwise untouched**: ghosts, pool slots past the active blocks,
!< components outside the range, masked-out cells; dq is exactly 0 where the mask is inactive. The flat kernels are
!< checked against plain loops. Run with -fopenmp too (tests-gnu-omp) for the parallel loops.
use selah
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: NV=3, NI=4, NJ=3, NK=2, NGC=2, NB=6, NM=2, NR=3 !< Components, cells, ghosts, pool, mask
                                                                            !< components, registers.
real(R_P),    parameter :: TOL = 16._R_P * epsilon(1._R_P)                  !< Value tolerance.
real(R_P),    parameter :: DT = 0.37_R_P                                    !< Time step.
real(R_P),    parameter :: COEF(NR) = [0.25_R_P, -0.5_R_P, 1.125_R_P]       !< Coefficients.
integer(I4P) :: layout, blocks
character(7) :: lname

do layout=LAYOUT_V_FIRST, LAYOUT_V_LAST
   lname = merge('v-first', 'v-last ', layout == LAYOUT_V_FIRST)
   do blocks=3, 5, 2
      print '(A,I0,A)', trim(lname)//', ', blocks, ' active blocks:'
      call test_box(box_interior(layout, NV, NI, NJ, NK, NGC, blocks), trim(lname))
   enddo
   call test_partition(box_interior(layout, NV, NI, NJ, NK, NGC, 4_I4P), trim(lname))
   call test_default_mask(box_interior(layout, NV, NI, NJ, NK, NGC, 4_I4P), trim(lname))
   call test_2n_first_stage(box_interior(layout, NV, NI, NJ, NK, NGC, 4_I4P), trim(lname))
enddo
call test_flat

if (failures > 0) error stop 'selah_test_kernels: failed'
print '(A)', 'selah_test_kernels: ok'

contains
   subroutine allocate_field(bx, a, nc, seed)
   !< Allocate a field of the pool (ghosts included, nc components) and fill it with values in [-1, 1].
   type(box),              intent(in)    :: bx    !< Box (layout, lower bounds).
   real(R_P), allocatable, intent(inout) :: a(:,:,:,:,:) !< Field.
   integer(I4P),           intent(in)    :: nc    !< Components.
   integer(I4P),           intent(in)    :: seed  !< Seed of the values.
   integer(I4P)                          :: ub(5) !< Upper bounds.

   if (bx%layout == LAYOUT_V_FIRST) then
      ub = [nc, NI+NGC, NJ+NGC, NK+NGC, NB]
   else
      ub = [NB, NI+NGC, NJ+NGC, NK+NGC, nc]
   endif
   if (allocated(a)) deallocate(a)
   allocate(a(bx%lb(1):ub(1), bx%lb(2):ub(2), bx%lb(3):ub(3), bx%lb(4):ub(4), bx%lb(5):ub(5)))
   call fill(a, seed)
   endsubroutine allocate_field

   subroutine fill(a, seed)
   !< Deterministic values in [-1, 1] (a linear congruential sequence: no compiler-dependent random numbers).
   real(R_P),    intent(inout) :: a(:,:,:,:,:) !< Field.
   integer(I4P), intent(in)    :: seed         !< Seed.
   integer(I8P)                :: s            !< State.
   integer(I4P)                :: i1, i2, i3, i4, i5

   s = int(seed, I8P) * 7919_I8P + 1_I8P
   do i5=1, size(a,5) ; do i4=1, size(a,4) ; do i3=1, size(a,3) ; do i2=1, size(a,2) ; do i1=1, size(a,1)
      s = modulo(s * 1103515245_I8P + 12345_I8P, 2147483648_I8P)
      a(i1,i2,i3,i4,i5) = real(s, R_P) / 1073741824._R_P - 1._R_P
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine fill

   pure logical function inside(bx, i)
   !< The cell (i1..i5) is in the loop range of the box.
   type(box),    intent(in) :: bx   !< Box.
   integer(I4P), intent(in) :: i(5) !< Indices.

   inside = all(i >= bx%lo) .and. all(i <= bx%hi)
   endfunction inside

   pure logical function active(bx, mask, mc, i)
   !< The mask is negative at the cell (component mc of the mask's component dimension).
   type(box),    intent(in) :: bx                !< Box.
   real(R_P),    intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Mask.
   integer(I4P), intent(in) :: mc                !< Mask component.
   integer(I4P), intent(in) :: i(5)              !< Indices.

   if (bx%layout == LAYOUT_V_FIRST) then
      active = mask(mc,i(2),i(3),i(4),i(5)) < 0._R_P
   else
      active = mask(i(1),i(2),i(3),i(4),mc) < 0._R_P
   endif
   endfunction active

   subroutine compare(bx, got, before, expected, mode, mask, mc, label)
   !< Updated cells (in the box, and active if mode > 0) hold `expected` within TOL; all other cells hold `before`
   !< bitwise.
   type(box),    intent(in) :: bx       !< Box.
   real(R_P),    intent(in) :: got(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)      !< Kernel result.
   real(R_P),    intent(in) :: before(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)   !< Before the kernel.
   real(R_P),    intent(in) :: expected(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Reference.
   integer(I4P), intent(in) :: mode     !< 0: unmasked; 1: masked.
   real(R_P),    intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)     !< Mask.
   integer(I4P), intent(in) :: mc       !< Mask component.
   character(*), intent(in) :: label    !< What is checked.
   logical                  :: values_ok, untouched_ok, upd
   integer(I4P)             :: i1, i2, i3, i4, i5, updated

   values_ok = .true. ; untouched_ok = .true. ; updated = 0
   do i5=lbound(got,5), ubound(got,5) ; do i4=lbound(got,4), ubound(got,4) ; do i3=lbound(got,3), ubound(got,3)
   do i2=lbound(got,2), ubound(got,2) ; do i1=lbound(got,1), ubound(got,1)
      upd = inside(bx, [i1,i2,i3,i4,i5])
      if (upd .and. mode > 0) upd = active(bx, mask, mc, [i1,i2,i3,i4,i5])
      if (upd) then
         updated = updated + 1
         if (abs(got(i1,i2,i3,i4,i5) - expected(i1,i2,i3,i4,i5)) > TOL) values_ok = .false.
      else
         if (got(i1,i2,i3,i4,i5) /= before(i1,i2,i3,i4,i5)) untouched_ok = .false.
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   call check(values_ok .and. updated > 0, label//': updated cells hold the reference value')
   call check(untouched_ok, label//': every other cell bitwise untouched')
   endsubroutine compare

   function masked_reference(bx, e, q, mask, mc) result(m)
   !< e in the active cells of the box, q in its inactive ones, e elsewhere (compare checks the box cells only).
   type(box),    intent(in) :: bx
   real(R_P),    intent(in) :: e(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),    intent(in) :: q(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),    intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   integer(I4P), intent(in) :: mc
   real(R_P), allocatable   :: m(:,:,:,:,:)
   integer(I4P)             :: i1, i2, i3, i4, i5

   m = e
   do i5=lbound(e,5), ubound(e,5) ; do i4=lbound(e,4), ubound(e,4) ; do i3=lbound(e,3), ubound(e,3)
   do i2=lbound(e,2), ubound(e,2) ; do i1=lbound(e,1), ubound(e,1)
      if (inside(bx, [i1,i2,i3,i4,i5])) then
         if (.not. active(bx, mask, mc, [i1,i2,i3,i4,i5])) m(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5)
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   endfunction masked_reference

   subroutine test_box(bx, lname)
   !< Every box kernel, unmasked and masked.
   type(box),    intent(in) :: bx    !< Box.
   character(*), intent(in) :: lname !< Layout name.
   real(R_P), allocatable   :: x(:,:,:,:,:), y(:,:,:,:,:), y0(:,:,:,:,:), e(:,:,:,:,:), mask(:,:,:,:,:)
   real(R_P), allocatable   :: q(:,:,:,:,:), u(:,:,:,:,:), dq(:,:,:,:,:), dq0(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P), allocatable   :: regs(:,:,:,:,:,:)
   integer(I4P)             :: mode, mc, k, n

   call allocate_field(bx, x, NV, 1_I4P) ; call allocate_field(bx, y0, NV, 2_I4P) ; call allocate_field(bx, q, NV, 3_I4P)
   call allocate_field(bx, r, NV, 4_I4P) ; call allocate_field(bx, mask, NM, 5_I4P)
   call allocate_field(bx, dq0, NV, 9_I4P)
   allocate(regs(lbound(x,1):ubound(x,1), lbound(x,2):ubound(x,2), lbound(x,3):ubound(x,3), lbound(x,4):ubound(x,4), &
                 lbound(x,5):ubound(x,5), NR))
   do k=1, NR
      call fill(regs(:,:,:,:,:,k), 10_I4P + k)
   enddo
   do mode=0, 2 ! 0 unmasked, 1 masked (component NM, the last), 2 masked (component 1)
      mc = merge(1_I4P, NM, mode == 2)
      associate(tag=>lname//merge(' unmasked        ', ' masked          ', mode == 0)//merge(' comp 1', '       ', mode == 2))
      ! copy
      y = y0 ; e = x
      if (mode == 0) then ; call selah_copy(bx, x, y) ; else ; call selah_copy(bx, x, y, mask, mc) ; endif
      call compare(bx, y, y0, e, min(mode, 1_I4P), mask, mc, 'copy '//tag)
      ! axpy
      y = y0 ; e = y0 + 0.7_R_P * x
      if (mode == 0) then ; call selah_axpy(bx, 0.7_R_P, x, y) ; else ; call selah_axpy(bx, 0.7_R_P, x, y, mask, mc) ; endif
      call compare(bx, y, y0, e, min(mode, 1_I4P), mask, mc, 'axpy '//tag)
      ! lincomb
      y = y0 ; e = 0.6_R_P * x - 0.4_R_P * y0
      if (mode == 0) then
         call selah_lincomb(bx, 0.6_R_P, x, -0.4_R_P, y)
      else
         call selah_lincomb(bx, 0.6_R_P, x, -0.4_R_P, y, mask, mc)
      endif
      call compare(bx, y, y0, e, min(mode, 1_I4P), mask, mc, 'lincomb '//tag)
      ! stage_state, 0 and 3 registers
      do n=0, NR, NR
         u = y0 ; e = q
         do k=1, n
            e = e + DT * COEF(k) * regs(:,:,:,:,:,k)
         enddo
         if (mode == 0) then
            call selah_stage_state(bx, n, DT, COEF, q, regs, u)
         else
            call selah_stage_state(bx, n, DT, COEF, q, regs, u, mask, mc)
         endif
         call compare(bx, u, y0, e, min(mode, 1_I4P), mask, mc, 'stage_state n='//achar(48+n)//' '//tag)
         if (mode > 0) then ! the stepper's mask policy: the inactive cells of the box take q
            u = y0
            call selah_stage_state(bx, n, DT, COEF, q, regs, u, mask, mc, masked_q=.true.)
            call compare(bx, u, y0, masked_reference(bx, e, q, mask, mc), 0_I4P, mask, mc, &
                         'stage_state masked_q n='//achar(48+n)//' '//tag)
         endif
      enddo
      ! step_update, with the step residual
      y = y0 ; dq = dq0 ; e = 0._R_P
      do k=1, NR
         e = e + COEF(k) * regs(:,:,:,:,:,k)
      enddo
      if (mode == 0) then
         call selah_step_update(bx, NR, DT, COEF, regs, y, dq)
      else
         call selah_step_update(bx, NR, DT, COEF, regs, y, dq, mask, mc)
      endif
      call compare(bx, y, y0, y0 + DT * e, min(mode, 1_I4P), mask, mc, 'step_update '//tag)
      if (mode == 0) then
         call compare(bx, dq, dq0, e, 0_I4P, mask, mc, 'step_update dq '//tag)
      else
         call check(dq_masked_ok(bx, dq, dq0, e, mask, mc), 'step_update dq '//tag//': d active, 0 inactive, else untouched')
      endif
      ! step_update, without
      y = y0
      if (mode == 0) then
         call selah_step_update(bx, NR, DT, COEF, regs, y)
      else
         call selah_step_update(bx, NR, DT, COEF, regs, y, mask=mask, mask_comp=mc)
      endif
      call compare(bx, y, y0, y0 + DT * e, min(mode, 1_I4P), mask, mc, 'step_update (no dq) '//tag)
      ! shu_osher_2r_stage
      y = y0 ; e = 0.75_R_P * q + 0.25_R_P * y0 + DT * 0.25_R_P * r
      if (mode == 0) then
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y)
      else
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y, mask, mc)
      endif
      call compare(bx, y, y0, e, min(mode, 1_I4P), mask, mc, 'shu_osher_2r_stage '//tag)
      ! williamson_2n_stage: both outputs
      y = y0 ; dq = dq0
      if (mode == 0) then
         call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r, dq, y)
      else
         call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r, dq, y, mask, mc)
      endif
      e = -0.5_R_P * dq0 + DT * r
      call compare(bx, dq, dq0, e, min(mode, 1_I4P), mask, mc, 'williamson_2n_stage dQ '//tag)
      call compare(bx, y, y0, y0 + 0.8_R_P * e, min(mode, 1_I4P), mask, mc, 'williamson_2n_stage q '//tag)
      endassociate
   enddo
   endsubroutine test_box

   logical function dq_masked_ok(bx, dq, dq0, e, mask, mc) result(ok)
   !< dq: e where in the box and active, exactly 0 where in the box and inactive, untouched outside the box.
   type(box),    intent(in) :: bx
   real(R_P),    intent(in) :: dq(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),    intent(in) :: dq0(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),    intent(in) :: e(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),    intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   integer(I4P), intent(in) :: mc
   integer(I4P)             :: i1, i2, i3, i4, i5, nact, ninact

   ok = .true. ; nact = 0 ; ninact = 0
   do i5=lbound(dq,5), ubound(dq,5) ; do i4=lbound(dq,4), ubound(dq,4) ; do i3=lbound(dq,3), ubound(dq,3)
   do i2=lbound(dq,2), ubound(dq,2) ; do i1=lbound(dq,1), ubound(dq,1)
      if (.not. inside(bx, [i1,i2,i3,i4,i5])) then
         if (dq(i1,i2,i3,i4,i5) /= dq0(i1,i2,i3,i4,i5)) ok = .false.
      elseif (active(bx, mask, mc, [i1,i2,i3,i4,i5])) then
         nact = nact + 1
         if (abs(dq(i1,i2,i3,i4,i5) - e(i1,i2,i3,i4,i5)) > TOL) ok = .false.
      else
         ninact = ninact + 1
         if (dq(i1,i2,i3,i4,i5) /= 0._R_P) ok = .false.
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   ok = ok .and. nact > 0 .and. ninact > 0
   endfunction dq_masked_ok

   subroutine test_partition(bx, lname)
   !< A component range (a partition of a splitting scheme): every kernel updates components 2..3 only.
   type(box),    intent(in) :: bx    !< Box of all components.
   character(*), intent(in) :: lname !< Layout name.
   real(R_P), allocatable   :: x(:,:,:,:,:), y(:,:,:,:,:), y0(:,:,:,:,:), mask(:,:,:,:,:), q(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P), allocatable   :: dq(:,:,:,:,:), dq0(:,:,:,:,:), e(:,:,:,:,:), regs(:,:,:,:,:,:)
   type(box)                :: sub
   integer(I4P)             :: k

   call allocate_field(bx, x, NV, 21_I4P) ; call allocate_field(bx, y0, NV, 22_I4P) ; call allocate_field(bx, mask, NM, 23_I4P)
   call allocate_field(bx, q, NV, 24_I4P) ; call allocate_field(bx, r, NV, 25_I4P) ; call allocate_field(bx, dq0, NV, 26_I4P)
   allocate(regs(lbound(x,1):ubound(x,1), lbound(x,2):ubound(x,2), lbound(x,3):ubound(x,3), lbound(x,4):ubound(x,4), &
                 lbound(x,5):ubound(x,5), NR))
   do k=1, NR
      call fill(regs(:,:,:,:,:,k), 30_I4P + k)
   enddo
   sub = box_components(bx, 2_I4P, 3_I4P)
   call check(sub%lo(component_dim(bx)) == 2 .and. sub%hi(component_dim(bx)) == 3, lname//' partition: component range set')
   associate(tag=>lname//' partition, components 2..3')
   y = y0 ; call selah_copy(sub, x, y) ; call compare(sub, y, y0, x, 0_I4P, mask, 1_I4P, 'copy '//tag)
   y = y0 ; call selah_axpy(sub, 1.5_R_P, x, y)
   call compare(sub, y, y0, y0 + 1.5_R_P * x, 0_I4P, mask, 1_I4P, 'axpy '//tag)
   y = y0 ; call selah_lincomb(sub, 0.6_R_P, x, -0.4_R_P, y)
   call compare(sub, y, y0, 0.6_R_P * x - 0.4_R_P * y0, 0_I4P, mask, 1_I4P, 'lincomb '//tag)
   e = q
   do k=1, NR
      e = e + DT * COEF(k) * regs(:,:,:,:,:,k)
   enddo
   y = y0 ; call selah_stage_state(sub, NR, DT, COEF, q, regs, y) ; call compare(sub, y, y0, e, 0_I4P, mask, 1_I4P, &
                                                                                 'stage_state '//tag)
   e = 0._R_P
   do k=1, NR
      e = e + COEF(k) * regs(:,:,:,:,:,k)
   enddo
   y = y0 ; dq = dq0 ; call selah_step_update(sub, NR, DT, COEF, regs, y, dq)
   call compare(sub, y, y0, y0 + DT * e, 0_I4P, mask, 1_I4P, 'step_update '//tag)
   call compare(sub, dq, dq0, e, 0_I4P, mask, 1_I4P, 'step_update dq '//tag)
   y = y0 ; call selah_shu_osher_2r_stage(sub, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y)
   call compare(sub, y, y0, 0.75_R_P * q + 0.25_R_P * y0 + DT * 0.25_R_P * r, 0_I4P, mask, 1_I4P, 'shu_osher_2r '//tag)
   y = y0 ; dq = dq0 ; call selah_williamson_2n_stage(sub, DT, -0.5_R_P, 0.8_R_P, r, dq, y)
   call compare(sub, dq, dq0, -0.5_R_P * dq0 + DT * r, 0_I4P, mask, 1_I4P, 'williamson_2n dQ '//tag)
   call compare(sub, y, y0, y0 + 0.8_R_P * (-0.5_R_P * dq0 + DT * r), 0_I4P, mask, 1_I4P, 'williamson_2n q '//tag)
   endassociate
   endsubroutine test_partition

   subroutine test_default_mask(bx, lname)
   !< mask_comp absent: every masked kernel uses the last mask component, bitwise the result with mask_comp = NM.
   type(box),    intent(in) :: bx    !< Box.
   character(*), intent(in) :: lname !< Layout name.
   real(R_P), allocatable   :: x(:,:,:,:,:), y0(:,:,:,:,:), mask(:,:,:,:,:), q(:,:,:,:,:), r(:,:,:,:,:), dq0(:,:,:,:,:)
   real(R_P), allocatable   :: ya(:,:,:,:,:), yb(:,:,:,:,:), da(:,:,:,:,:), db(:,:,:,:,:), regs(:,:,:,:,:,:)
   integer(I4P)             :: k
   logical                  :: ok

   call allocate_field(bx, x, NV, 41_I4P) ; call allocate_field(bx, y0, NV, 42_I4P) ; call allocate_field(bx, mask, NM, 43_I4P)
   call allocate_field(bx, q, NV, 44_I4P) ; call allocate_field(bx, r, NV, 45_I4P) ; call allocate_field(bx, dq0, NV, 46_I4P)
   allocate(regs(lbound(x,1):ubound(x,1), lbound(x,2):ubound(x,2), lbound(x,3):ubound(x,3), lbound(x,4):ubound(x,4), &
                 lbound(x,5):ubound(x,5), NR))
   do k=1, NR
      call fill(regs(:,:,:,:,:,k), 50_I4P + k)
   enddo
   ok = .true.
   ya = y0 ; yb = y0 ; call selah_copy(bx, x, ya, mask) ; call selah_copy(bx, x, yb, mask, NM) ; ok = ok .and. all(ya == yb)
   ya = y0 ; yb = y0 ; call selah_axpy(bx, 0.7_R_P, x, ya, mask) ; call selah_axpy(bx, 0.7_R_P, x, yb, mask, NM)
   ok = ok .and. all(ya == yb)
   ya = y0 ; yb = y0 ; call selah_lincomb(bx, 0.6_R_P, x, -0.4_R_P, ya, mask)
   call selah_lincomb(bx, 0.6_R_P, x, -0.4_R_P, yb, mask, NM) ; ok = ok .and. all(ya == yb)
   ya = y0 ; yb = y0 ; call selah_stage_state(bx, NR, DT, COEF, q, regs, ya, mask)
   call selah_stage_state(bx, NR, DT, COEF, q, regs, yb, mask, NM) ; ok = ok .and. all(ya == yb)
   ya = y0 ; yb = y0 ; da = dq0 ; db = dq0 ; call selah_step_update(bx, NR, DT, COEF, regs, ya, da, mask)
   call selah_step_update(bx, NR, DT, COEF, regs, yb, db, mask, NM) ; ok = ok .and. all(ya == yb) .and. all(da == db)
   ya = y0 ; yb = y0 ; call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, ya, mask)
   call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, yb, mask, NM) ; ok = ok .and. all(ya == yb)
   ya = y0 ; yb = y0 ; da = dq0 ; db = dq0 ; call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r, da, ya, mask)
   call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r, db, yb, mask, NM)
   ok = ok .and. all(ya == yb) .and. all(da == db)
   call check(ok, lname//' every masked kernel: mask_comp absent = the last component, bitwise')
   endsubroutine test_default_mask

   subroutine test_2n_first_stage(bx, lname)
   !< The 2N stage with a = 0 never reads its register: a NaN register gives dQ = dt r and q + b dt r in the box.
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   type(box),    intent(in) :: bx    !< Box.
   character(*), intent(in) :: lname !< Layout name.
   real(R_P), allocatable   :: y(:,:,:,:,:), y0(:,:,:,:,:), r(:,:,:,:,:), dq(:,:,:,:,:), e(:,:,:,:,:), mask(:,:,:,:,:)
   real(R_P), allocatable   :: dq0(:,:,:,:,:)

   call allocate_field(bx, y0, NV, 61_I4P) ; call allocate_field(bx, r, NV, 62_I4P) ; call allocate_field(bx, mask, NM, 63_I4P)
   call allocate_field(bx, dq0, NV, 64_I4P)
   ! finite register: the values of the general expression's case a = 0, and every other cell untouched
   y = y0 ; dq = dq0 ; call selah_williamson_2n_stage(bx, DT, 0._R_P, 0.8_R_P, r, dq, y)
   e = DT * r
   call compare(bx, dq, dq0, e, 0_I4P, mask, 1_I4P, lname//' williamson_2n a=0 dQ')
   call compare(bx, y, y0, y0 + 0.8_R_P * e, 0_I4P, mask, 1_I4P, lname//' williamson_2n a=0 q')
   ! NaN register: the box comes out finite (the register is written, never read)
   y = y0 ; dq = y0 ; dq = ieee_value(1._R_P, ieee_quiet_nan) ; call selah_williamson_2n_stage(bx, DT, 0._R_P, 0.8_R_P, r, dq, y)
   call check(all(ieee_is_finite(y(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),bx%lo(4):bx%hi(4), &
                                   bx%lo(5):bx%hi(5)))), lname//' williamson_2n a=0: a NaN register is not read')
   endsubroutine test_2n_first_stage

   subroutine test_flat
   !< The flat kernels against plain loops.
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   integer(I4P), parameter :: N = 37
   real(R_P)               :: x(N), y(N), y0(N), q(N), r(N), dq(N), dqr(N), dqr0(N), regs(N,NR), e(N)
   integer(I4P)            :: k

   x = [(sin(real(k, R_P)), k=1, N)] ; y0 = [(cos(real(k, R_P)), k=1, N)] ; q = [(0.1_R_P*k, k=1, N)]
   r = [(1._R_P/k, k=1, N)] ; dqr0 = [(-0.01_R_P*k, k=1, N)]
   do k=1, NR
      regs(:,k) = x * real(k, R_P) - y0
   enddo
   print '(A)', 'flat:'
   y = y0 ; call selah_copy(x, y) ; call check(all(y == x), 'flat copy')
   y = y0 ; call selah_axpy(0.7_R_P, x, y) ; call check(all(abs(y - (y0 + 0.7_R_P*x)) <= TOL), 'flat axpy')
   y = y0 ; call selah_lincomb(0.6_R_P, x, -0.4_R_P, y)
   call check(all(abs(y - (0.6_R_P*x - 0.4_R_P*y0)) <= TOL), 'flat lincomb')
   e = q
   do k=1, NR
      e = e + DT * COEF(k) * regs(:,k)
   enddo
   y = y0 ; call selah_stage_state(NR, DT, COEF, q, regs, y) ; call check(all(abs(y - e) <= TOL), 'flat stage_state')
   y = y0 ; call selah_stage_state(0_I4P, DT, COEF, q, regs, y) ; call check(all(y == q), 'flat stage_state n=0 (a copy)')
   e = 0._R_P
   do k=1, NR
      e = e + COEF(k) * regs(:,k)
   enddo
   y = y0 ; call selah_step_update(NR, DT, COEF, regs, y, dq)
   call check(all(abs(y - (y0 + DT*e)) <= TOL) .and. all(abs(dq - e) <= TOL), 'flat step_update with dq')
   y = y0 ; call selah_step_update(NR, DT, COEF, regs, y)
   call check(all(abs(y - (y0 + DT*e)) <= TOL), 'flat step_update without dq')
   y = y0 ; call selah_shu_osher_2r_stage(DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y)
   call check(all(abs(y - (0.75_R_P*q + 0.25_R_P*y0 + DT*0.25_R_P*r)) <= TOL), 'flat shu_osher_2r_stage')
   y = y0 ; dqr = dqr0 ; call selah_williamson_2n_stage(DT, -0.5_R_P, 0.8_R_P, r, dqr, y)
   call check(all(abs(dqr - (-0.5_R_P*dqr0 + DT*r)) <= TOL) .and. all(abs(y - (y0 + 0.8_R_P*dqr)) <= TOL), &
              'flat williamson_2n_stage')
   y = y0 ; dqr = ieee_value(1._R_P, ieee_quiet_nan) ! 0 * NaN would be NaN, and trap in a debug build
   call selah_williamson_2n_stage(DT, 0._R_P, 0.8_R_P, r, dqr, y)
   call check(all(dqr == DT*r) .and. all(abs(y - (y0 + 0.8_R_P*DT*r)) <= TOL), 'flat williamson_2n_stage a=0: register not read')
   endsubroutine test_flat
endprogram selah_test_kernels
