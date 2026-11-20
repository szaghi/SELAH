!< SELAH test: the kernels on device arrays (FUNDAL), value by value and cell by cell.

#include "fundal.H"

program selah_test_device_kernels
!< SELAH test: the kernels on device arrays (FUNDAL), value by value and cell by cell.
!<
!< The arrays are allocated and filled through FUNDAL (device memory in an OpenACC build, host memory otherwise), every
!< kernel runs on them, and the results are copied back: each cell of the box holds the host reference (16 eps: the
!< device may contract to FMA), every other cell is bitwise untouched (ghosts, pool slots past the active blocks,
!< masked-out cells), dq is exactly 0 where the mask is inactive. Both layouts, unmasked and masked, 3 of 6 blocks
!< active; then the stage state with no register (a copy), the step update without dq, the 2N stage with a = 0, a
!< component range (box_components) and an explicit mask component; the stage state with masked_q (inactive cells
!< of the box = q); every flat kernel too. In a device build
!< `dev_init(require_device=.true.)` forbids a silent host fallback.
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device, dev_assign_to_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: NV=3, NI=5, NJ=4, NK=3, NGC=2, NB=6, BN=3, NM=2, NR=3
real(R_P),    parameter :: TOL = 16._R_P * epsilon(1._R_P)
real(R_P),    parameter :: DT = 0.37_R_P
real(R_P),    parameter :: COEF(NR) = [0.25_R_P, -0.5_R_P, 1.125_R_P]
real(R_P), pointer      :: coef_d(:) => null()
real(R_P)               :: coef_h(NR) !< COEF as a variable, the source of the device copy.
integer(I4P)            :: layout

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build: kernels run on the device'
#else
print '(A)', 'host build: FUNDAL host memory'
#endif
coef_h = COEF ! a variable: passing the named constant itself down FUNDAL to a target dummy left nvfortran 26.9 -Mbounds
               ! reading garbage bounds for COEF afterwards (debug builds only; release, gfortran and ifx are fine)
call dev_assign_to_device(dst=coef_d, src=coef_h)
do layout=LAYOUT_V_FIRST, LAYOUT_V_LAST
   call test_layout(box_interior(layout, NV, NI, NJ, NK, NGC, BN), merge('v-first', 'v-last ', layout == LAYOUT_V_FIRST))
enddo
call test_flat
call dev_free(coef_d)

if (failures > 0) error stop 'selah_test_device_kernels: failed'
print '(A)', 'selah_test_device_kernels: ok'

contains
   subroutine host_field(bx, a, nc, seed)
   !< A host field of the pool (ghosts included, nc components), values in [-1, 1].
   type(box),              intent(in)    :: bx
   real(R_P), allocatable, intent(inout) :: a(:,:,:,:,:)
   integer(I4P),           intent(in)    :: nc, seed
   integer(I4P)                          :: ub(5), i1, i2, i3, i4, i5
   integer(I8P)                          :: s

   ub = merge([nc, NI+NGC, NJ+NGC, NK+NGC, NB], [NB, NI+NGC, NJ+NGC, NK+NGC, nc], bx%layout == LAYOUT_V_FIRST)
   if (allocated(a)) deallocate(a)
   allocate(a(bx%lb(1):ub(1), bx%lb(2):ub(2), bx%lb(3):ub(3), bx%lb(4):ub(4), bx%lb(5):ub(5)))
   s = int(seed, I8P) * 7919_I8P + 1_I8P
   do i5=lbound(a,5), ubound(a,5) ; do i4=lbound(a,4), ubound(a,4) ; do i3=lbound(a,3), ubound(a,3)
   do i2=lbound(a,2), ubound(a,2) ; do i1=lbound(a,1), ubound(a,1)
      s = modulo(s * 1103515245_I8P + 12345_I8P, 2147483648_I8P)
      a(i1,i2,i3,i4,i5) = real(s, R_P) / 1073741824._R_P - 1._R_P
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine host_field

   subroutine to_device(a, a_d)
   !< A device copy of a host field, same bounds.
   real(R_P),          intent(in)    :: a(:,:,:,:,:)
   real(R_P), pointer, intent(inout) :: a_d(:,:,:,:,:)
   integer(I4P)                      :: ierr

   call dev_alloc(fptr_dev=a_d, ubounds=int(ubound(a), I4P), lbounds=int(lbound(a), I4P), ierr=ierr)
   if (ierr /= 0) error stop 'selah_test_device_kernels: dev_alloc failed'
   call dev_memcpy_to_device(dst=a_d, src=a)
   endsubroutine to_device

   logical function updated(bx, mask, masked, i)
   !< The kernel updates this cell: in the box and, when masked, active (last mask component < 0).
   type(box),    intent(in) :: bx
   real(R_P),    intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   logical,      intent(in) :: masked
   integer(I4P), intent(in) :: i(5)

   updated = all(i >= bx%lo) .and. all(i <= bx%hi)
   if (updated .and. masked) then
      if (bx%layout == LAYOUT_V_FIRST) then
         updated = mask(NM,i(2),i(3),i(4),i(5)) < 0._R_P
      else
         updated = mask(i(1),i(2),i(3),i(4),NM) < 0._R_P
      endif
   endif
   endfunction updated

   subroutine compare(bx, mask, masked, got_d, before, expected, label, zero_inactive)
   !< Copy the device result back: updated cells within TOL of expected, the others bitwise as before (or exactly 0
   !< for inactive cells of the box when zero_inactive).
   type(box),          intent(in) :: bx
   real(R_P),          intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   logical,            intent(in) :: masked
   real(R_P), pointer, intent(in) :: got_d(:,:,:,:,:)
   real(R_P),          intent(in) :: before(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),          intent(in) :: expected(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   character(*),       intent(in) :: label
   logical,  optional, intent(in) :: zero_inactive
   real(R_P), allocatable         :: got(:,:,:,:,:)
   logical                        :: values_ok, others_ok, in_box
   integer(I4P)                   :: i1, i2, i3, i4, i5, n

   allocate(got, mold=before)
   call dev_memcpy_from_device(dst=got, src=got_d)
   values_ok = .true. ; others_ok = .true. ; n = 0
   do i5=lbound(got,5), ubound(got,5) ; do i4=lbound(got,4), ubound(got,4) ; do i3=lbound(got,3), ubound(got,3)
   do i2=lbound(got,2), ubound(got,2) ; do i1=lbound(got,1), ubound(got,1)
      if (updated(bx, mask, masked, [i1,i2,i3,i4,i5])) then
         n = n + 1
         if (abs(got(i1,i2,i3,i4,i5) - expected(i1,i2,i3,i4,i5)) > TOL) values_ok = .false.
      else
         in_box = all([i1,i2,i3,i4,i5] >= bx%lo) .and. all([i1,i2,i3,i4,i5] <= bx%hi)
         if (present(zero_inactive) .and. in_box) then
            if (got(i1,i2,i3,i4,i5) /= 0._R_P) others_ok = .false.
         elseif (got(i1,i2,i3,i4,i5) /= before(i1,i2,i3,i4,i5)) then
            others_ok = .false.
         endif
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   call check(values_ok .and. n > 0 .and. others_ok, label)
   endsubroutine compare

   subroutine test_layout(bx, lname)
   !< Every box kernel on device arrays, unmasked and masked.
   type(box),    intent(in) :: bx
   character(*), intent(in) :: lname
   real(R_P), allocatable   :: x(:,:,:,:,:), y0(:,:,:,:,:), q(:,:,:,:,:), r(:,:,:,:,:), mask(:,:,:,:,:), e(:,:,:,:,:)
   real(R_P), allocatable   :: regs(:,:,:,:,:,:), sum_b(:,:,:,:,:)
   real(R_P), allocatable   :: em(:,:,:,:,:)
   integer(I4P)             :: i1, i2, i3, i4, i5
   real(R_P), pointer       :: x_d(:,:,:,:,:)=>null(), y_d(:,:,:,:,:)=>null(), q_d(:,:,:,:,:)=>null()
   real(R_P), pointer       :: r_d(:,:,:,:,:)=>null(), m_d(:,:,:,:,:)=>null(), dq_d(:,:,:,:,:)=>null()
   real(R_P), pointer       :: regs_d(:,:,:,:,:,:)=>null()
   integer(I4P)             :: k, ierr
   logical                  :: masked
   character(:), allocatable :: tag

   call host_field(bx, x, NV, 1_I4P) ; call host_field(bx, y0, NV, 2_I4P) ; call host_field(bx, q, NV, 3_I4P)
   call host_field(bx, r, NV, 4_I4P) ; call host_field(bx, mask, NM, 5_I4P)
   allocate(regs(lbound(x,1):ubound(x,1), lbound(x,2):ubound(x,2), lbound(x,3):ubound(x,3), lbound(x,4):ubound(x,4), &
                 lbound(x,5):ubound(x,5), NR))
   do k=1, NR
      call host_field(bx, e, NV, 10_I4P + k) ; regs(:,:,:,:,:,k) = e
   enddo
   sum_b = 0._R_P * x
   do k=1, NR
      sum_b = sum_b + COEF(k) * regs(:,:,:,:,:,k)
   enddo
   call to_device(x, x_d) ; call to_device(q, q_d) ; call to_device(r, r_d) ; call to_device(mask, m_d)
   call to_device(y0, y_d) ; call to_device(y0, dq_d)
   call dev_alloc(fptr_dev=regs_d, ubounds=int(ubound(regs), I4P), lbounds=int(lbound(regs), I4P), ierr=ierr)
   call dev_memcpy_to_device(dst=regs_d, src=regs)
   do k=0, 1
      masked = k == 1
      tag = lname//merge(' masked  ', ' unmasked', masked)//': '
      call dev_memcpy_to_device(dst=y_d, src=y0)
      if (masked) then ; call selah_copy(bx, x_d, y_d, m_d) ; else ; call selah_copy(bx, x_d, y_d) ; endif
      call compare(bx, mask, masked, y_d, y0, x, tag//'copy')
      call dev_memcpy_to_device(dst=y_d, src=y0)
      if (masked) then ; call selah_axpy(bx, 0.7_R_P, x_d, y_d, m_d) ; else ; call selah_axpy(bx, 0.7_R_P, x_d, y_d) ; endif
      call compare(bx, mask, masked, y_d, y0, y0 + 0.7_R_P * x, tag//'axpy')
      call dev_memcpy_to_device(dst=y_d, src=y0)
      if (masked) then
         call selah_lincomb(bx, 0.6_R_P, x_d, -0.4_R_P, y_d, m_d)
      else
         call selah_lincomb(bx, 0.6_R_P, x_d, -0.4_R_P, y_d)
      endif
      call compare(bx, mask, masked, y_d, y0, 0.6_R_P * x - 0.4_R_P * y0, tag//'lincomb')
      call dev_memcpy_to_device(dst=y_d, src=y0)
      if (masked) then
         call selah_stage_state(bx, NR, DT, coef_d, q_d, regs_d, y_d, m_d)
      else
         call selah_stage_state(bx, NR, DT, coef_d, q_d, regs_d, y_d)
      endif
      e = q
      do ierr=1, NR
         e = e + DT * COEF(ierr) * regs(:,:,:,:,:,ierr)
      enddo
      call compare(bx, mask, masked, y_d, y0, e, tag//'stage_state')
      if (masked) then ! the mask policy: the inactive cells of the box take q
         call dev_memcpy_to_device(dst=y_d, src=y0)
         call selah_stage_state(bx, NR, DT, coef_d, q_d, regs_d, y_d, m_d, masked_q=.true.)
         em = e
         do i5=bx%lo(5), bx%hi(5) ; do i4=bx%lo(4), bx%hi(4) ; do i3=bx%lo(3), bx%hi(3)
         do i2=bx%lo(2), bx%hi(2) ; do i1=bx%lo(1), bx%hi(1)
            if (.not. updated(bx, mask, .true., [i1,i2,i3,i4,i5])) em(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5)
         enddo ; enddo ; enddo ; enddo ; enddo
         call compare(bx, mask, .false., y_d, y0, em, tag//'stage_state masked_q (inactive cells = q)')
      endif
      call dev_memcpy_to_device(dst=y_d, src=y0) ; call dev_memcpy_to_device(dst=dq_d, src=y0)
      if (masked) then
         call selah_step_update(bx, NR, DT, coef_d, regs_d, y_d, dq_d, m_d)
      else
         call selah_step_update(bx, NR, DT, coef_d, regs_d, y_d, dq_d)
      endif
      call compare(bx, mask, masked, y_d, y0, y0 + DT * sum_b, tag//'step_update q')
      if (masked) then
         call compare(bx, mask, masked, dq_d, y0, sum_b, tag//'step_update dq (0 on inactive cells)', zero_inactive=.true.)
      else
         call compare(bx, mask, masked, dq_d, y0, sum_b, tag//'step_update dq')
      endif
      call dev_memcpy_to_device(dst=y_d, src=y0)
      if (masked) then
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q_d, r_d, y_d, m_d)
      else
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q_d, r_d, y_d)
      endif
      call compare(bx, mask, masked, y_d, y0, 0.75_R_P * q + 0.25_R_P * y0 + DT * 0.25_R_P * r, tag//'shu_osher_2r_stage')
      call dev_memcpy_to_device(dst=y_d, src=y0) ; call dev_memcpy_to_device(dst=dq_d, src=x)
      if (masked) then
         call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r_d, dq_d, y_d, m_d)
      else
         call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.8_R_P, r_d, dq_d, y_d)
      endif
      call compare(bx, mask, masked, dq_d, x, -0.5_R_P * x + DT * r, tag//'williamson_2n_stage dQ')
      call compare(bx, mask, masked, y_d, y0, y0 + 0.8_R_P * (-0.5_R_P * x + DT * r), tag//'williamson_2n_stage q')
   enddo
   tag = lname//' unmasked: '
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_stage_state(bx, 0_I4P, DT, coef_d, q_d, regs_d, y_d)
   call compare(bx, mask, .false., y_d, y0, q, tag//'stage_state n=0 (a copy)')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_step_update(bx, NR, DT, coef_d, regs_d, y_d)
   call compare(bx, mask, .false., y_d, y0, y0 + DT * sum_b, tag//'step_update without dq')
   call dev_memcpy_to_device(dst=y_d, src=y0) ; call dev_memcpy_to_device(dst=dq_d, src=x)
   call selah_williamson_2n_stage(bx, DT, 0._R_P, 0.8_R_P, r_d, dq_d, y_d)
   call compare(bx, mask, .false., dq_d, x, DT * r, tag//'williamson_2n_stage a=0 dQ')
   call compare(bx, mask, .false., y_d, y0, y0 + 0.8_R_P * DT * r, tag//'williamson_2n_stage a=0 q')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_axpy(box_components(bx, 2_I4P, 3_I4P), 0.7_R_P, x_d, y_d)
   call compare(box_components(bx, 2_I4P, 3_I4P), mask, .false., y_d, y0, y0 + 0.7_R_P * x, &
                tag//'axpy on components 2..3 (the others untouched)')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_axpy(bx, 0.7_R_P, x_d, y_d, m_d, 1_I4P)
   call compare_comp1(bx, mask, y0, y_d, y0 + 0.7_R_P * x, lname//' masked, mask_comp=1: axpy')
   call dev_free(x_d) ; call dev_free(y_d) ; call dev_free(q_d) ; call dev_free(r_d) ; call dev_free(m_d)
   call dev_free(dq_d) ; call dev_free(regs_d)
   endsubroutine test_layout

   subroutine compare_comp1(bx, mask, y0, got_d, expected, label)
   !< Masked with component 1: updated where mask(comp 1) < 0 in the box, the other cells bitwise as before.
   type(box),          intent(in) :: bx
   real(R_P),          intent(in) :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P),          intent(in) :: y0(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   real(R_P), pointer, intent(in) :: got_d(:,:,:,:,:)
   real(R_P),          intent(in) :: expected(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)
   character(*),       intent(in) :: label
   real(R_P), allocatable         :: got(:,:,:,:,:)
   logical                        :: ok, upd
   integer(I4P)                   :: i1, i2, i3, i4, i5, n

   allocate(got, mold=y0)
   call dev_memcpy_from_device(dst=got, src=got_d)
   ok = .true. ; n = 0
   do i5=lbound(got,5), ubound(got,5) ; do i4=lbound(got,4), ubound(got,4) ; do i3=lbound(got,3), ubound(got,3)
   do i2=lbound(got,2), ubound(got,2) ; do i1=lbound(got,1), ubound(got,1)
      upd = all([i1,i2,i3,i4,i5] >= bx%lo) .and. all([i1,i2,i3,i4,i5] <= bx%hi)
      if (upd) then
         if (bx%layout == LAYOUT_V_FIRST) then
            upd = mask(1,i2,i3,i4,i5) < 0._R_P
         else
            upd = mask(i1,i2,i3,i4,1) < 0._R_P
         endif
      endif
      if (upd) then
         n = n + 1
         if (abs(got(i1,i2,i3,i4,i5) - expected(i1,i2,i3,i4,i5)) > TOL) ok = .false.
      elseif (got(i1,i2,i3,i4,i5) /= y0(i1,i2,i3,i4,i5)) then
         ok = .false.
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   call check(ok .and. n > 0, label)
   endsubroutine compare_comp1

   subroutine test_flat
   !< The flat kernels on device arrays.
   integer(I4P), parameter :: N = 1001
   real(R_P)               :: x(N), y0(N), q(N), r(N), regs(N,NR), e(N), y(N), dq(N)
   real(R_P), pointer      :: x_d(:)=>null(), y_d(:)=>null(), q_d(:)=>null(), r_d(:)=>null(), dq_d(:)=>null()
   real(R_P), pointer      :: regs_d(:,:)=>null()
   integer(I4P)            :: k

   x = [(sin(real(k, R_P)), k=1, N)] ; y0 = [(cos(real(k, R_P)), k=1, N)] ; q = [(0.001_R_P*k, k=1, N)]
   r = [(1._R_P/k, k=1, N)]
   do k=1, NR
      regs(:,k) = x * real(k, R_P) - y0
   enddo
   call dev_assign_to_device(dst=x_d, src=x) ; call dev_assign_to_device(dst=y_d, src=y0)
   call dev_assign_to_device(dst=q_d, src=q) ; call dev_assign_to_device(dst=r_d, src=r)
   call dev_assign_to_device(dst=dq_d, src=y0) ; call dev_assign_to_device(dst=regs_d, src=regs)
   call selah_copy(x_d, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(y == x), 'flat copy')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_lincomb(0.6_R_P, x_d, -0.4_R_P, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(abs(y - (0.6_R_P*x - 0.4_R_P*y0)) <= TOL), 'flat lincomb')
   e = q
   do k=1, NR
      e = e + DT * COEF(k) * regs(:,k)
   enddo
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_stage_state(NR, DT, coef_d, q_d, regs_d, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(abs(y - e) <= TOL), 'flat stage_state')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_stage_state(0_I4P, DT, coef_d, q_d, regs_d, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(y == q), 'flat stage_state n=0 (a copy)')
   call dev_memcpy_to_device(dst=y_d, src=y0) ; call dev_memcpy_to_device(dst=dq_d, src=x)
   call selah_williamson_2n_stage(DT, -0.5_R_P, 0.8_R_P, r_d, dq_d, y_d)
   call dev_memcpy_from_device(dst=y, src=y_d) ; call dev_memcpy_from_device(dst=dq, src=dq_d)
   call check(all(abs(dq - (-0.5_R_P*x + DT*r)) <= TOL) .and. all(abs(y - (y0 + 0.8_R_P*dq)) <= TOL), &
              'flat williamson_2n_stage')
   call dev_memcpy_to_device(dst=y_d, src=y0) ; call dev_memcpy_to_device(dst=dq_d, src=x)
   call selah_williamson_2n_stage(DT, 0._R_P, 0.8_R_P, r_d, dq_d, y_d)
   call dev_memcpy_from_device(dst=y, src=y_d) ; call dev_memcpy_from_device(dst=dq, src=dq_d)
   call check(all(abs(dq - DT*r) <= TOL) .and. all(abs(y - (y0 + 0.8_R_P*DT*r)) <= TOL), 'flat williamson_2n_stage a=0')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_axpy(0.7_R_P, x_d, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(abs(y - (y0 + 0.7_R_P*x)) <= TOL), 'flat axpy')
   e = 0._R_P
   do k=1, NR
      e = e + COEF(k) * regs(:,k)
   enddo
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_step_update(NR, DT, coef_d, regs_d, y_d, dq_d)
   call dev_memcpy_from_device(dst=y, src=y_d) ; call dev_memcpy_from_device(dst=dq, src=dq_d)
   call check(all(abs(y - (y0 + DT*e)) <= TOL) .and. all(abs(dq - e) <= TOL), 'flat step_update with dq')
   call dev_memcpy_to_device(dst=y_d, src=y0)
   call selah_shu_osher_2r_stage(DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q_d, r_d, y_d) ; call dev_memcpy_from_device(dst=y, src=y_d)
   call check(all(abs(y - (0.75_R_P*q + 0.25_R_P*y0 + DT*0.25_R_P*r)) <= TOL), 'flat shu_osher_2r_stage')
   call dev_free(x_d) ; call dev_free(y_d) ; call dev_free(q_d) ; call dev_free(r_d) ; call dev_free(dq_d)
   call dev_free(regs_d)
   endsubroutine test_flat
endprogram selah_test_device_kernels
