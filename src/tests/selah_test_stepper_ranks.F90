!< SELAH test: the stepper on fields of rank 1 to 4, against the same data as rank 5, on the host or the device.

#include "fundal.H"

module selah_test_ranks_support
!< Fields of rank 1 to 5 over one contiguous storage, and the calls of a step dispatched on their rank.
!<
!< A field is a rank-5 storage `s` (FUNDAL memory: the device in a device build) with the bounds of its rank-5 view,
!< and a pointer of the field's own rank remapped onto it: the rank-k pointer is what a caller passes to the stepper.
!< The residuals are pointwise (no space operator): y' = -y^2, a constant, the indicator of one stage, an oscillator.
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
use selah
implicit none
public

integer(I4P), parameter :: RES_RICCATI = 1_I4P   !< r = -u^2.
integer(I4P), parameter :: RES_ONE     = 2_I4P   !< r = 1.
integer(I4P), parameter :: RES_STAGE   = 3_I4P   !< r = 1 at stage j, 0 at the others.
integer(I4P), parameter :: RES_OSC     = 4_I4P   !< r(1) = u(2), r(2) = -u(1) (two components).

type :: field
   !< A field of rank 1 to 5 over a contiguous rank-5 storage.
   integer(I4P)                    :: rank = 5_I4P          !< Rank of the field.
   logical                         :: comps = .true.        !< Its first dimension holds components.
   type(box)                       :: bx                    !< Box of the field's rank.
   integer(I4P)                    :: lb5(5), ub5(5)        !< Bounds of the storage (the rank-5 view).
   real(R_P), pointer, contiguous  :: s(:,:,:,:,:) => null() !< Storage.
   real(R_P), pointer              :: p1(:) => null(), p2(:,:) => null(), p3(:,:,:) => null(), p4(:,:,:,:) => null()
endtype field

contains
   subroutine make(f, rank, comps, lb, ub, lo, hi, ncomp)
   !< A field of rank `rank` allocated lb:ub, loop range lo:hi; ncomp overrides the components' upper bound (a mask).
   type(field),  intent(inout)        :: f
   integer(I4P), intent(in)           :: rank
   logical,      intent(in)           :: comps
   integer(I4P), intent(in)           :: lb(:), ub(:), lo(:), hi(:)
   integer(I4P), intent(in), optional :: ncomp
   integer(I4P)                       :: u(rank), d, ierr

   u = ub
   if (present(ncomp)) u(1) = ncomp
   f%rank = rank ; f%comps = comps
   if (rank == 5_I4P) then
      f%bx = box_range(lb, lo, hi, .true.)
   else
      f%bx = box_range(lb, lo, hi, comps)
   endif
   f%lb5 = 1_I4P ; f%ub5 = 1_I4P
   do d=1_I4P, rank
      f%lb5(box_dim(f%bx, d)) = lb(d)
      f%ub5(box_dim(f%bx, d)) = u(d)
   enddo
   call dev_alloc(fptr_dev=f%s, ubounds=f%ub5, lbounds=f%lb5, init_value=0._R_P, ierr=ierr)
   if (ierr /= 0) error stop 'selah_test_stepper_ranks: allocation failed'
   select case(rank)
   case(1) ; f%p1(lb(1):u(1)) => f%s
   case(2) ; f%p2(lb(1):u(1),lb(2):u(2)) => f%s
   case(3) ; f%p3(lb(1):u(1),lb(2):u(2),lb(3):u(3)) => f%s
   case(4) ; f%p4(lb(1):u(1),lb(2):u(2),lb(3):u(3),lb(4):u(4)) => f%s
   endselect
   endsubroutine make

   subroutine unmake(f)
   type(field), intent(inout) :: f
   if (associated(f%s)) call dev_free(f%s)
   nullify(f%s, f%p1, f%p2, f%p3, f%p4)
   endsubroutine unmake

   subroutine set_host(f, h)
   !< The storage from a host copy.
   type(field), intent(inout) :: f
   real(R_P),   intent(in)    :: h(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):)
   call dev_memcpy_to_device(dst=f%s, src=h)
   endsubroutine set_host

   function get_host(f) result(h)
   !< A host copy of the storage.
   type(field), intent(in) :: f
   real(R_P), allocatable  :: h(:,:,:,:,:)
   allocate(h(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4),f%lb5(5):f%ub5(5)))
   call dev_memcpy_from_device(dst=h, src=f%s)
   endfunction get_host

   subroutine init_stepper(st, f, scheme, form, partitions)
   !< The stepper for the field, at its rank.
   type(stepper), intent(inout)        :: st
   type(field),   intent(in)           :: f
   character(*),  intent(in)           :: scheme, form
   integer(I4P),  intent(in), optional :: partitions(2,2)
   select case(f%rank)
   case(1) ; call st%initialize(scheme, form, lbound(f%p1), ubound(f%p1), partitions, components=f%comps)
   case(2) ; call st%initialize(scheme, form, lbound(f%p2), ubound(f%p2), partitions, components=f%comps)
   case(3) ; call st%initialize(scheme, form, lbound(f%p3), ubound(f%p3), partitions, components=f%comps)
   case(4) ; call st%initialize(scheme, form, lbound(f%p4), ubound(f%p4), partitions, components=f%comps)
   case(5) ; call st%initialize(scheme, form, f%lb5, f%ub5, partitions)
   endselect
   endsubroutine init_stepper

   subroutine residual(f, kind, k, j, u, r)
   !< The pointwise residual of stage k over the cells of the box (rank-5 explicit shape: any rank by sequence
   !< association, as the stepper's bridges).
   type(field),  intent(in)    :: f
   integer(I4P), intent(in)    :: kind, k, j
   real(R_P),    intent(in)    :: u(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4), &
                                    f%lb5(5):f%ub5(5))
   real(R_P),    intent(inout) :: r(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4), &
                                    f%lb5(5):f%ub5(5))
   integer(I4P)                :: i1, i2, i3, i4, i5, lo(5), hi(5)
   real(R_P)                   :: c

   lo = f%bx%lo ; hi = f%bx%hi
   c = merge(1._R_P, 0._R_P, kind == RES_ONE .or. (kind == RES_STAGE .and. k == j))
   !$acc parallel loop independent collapse(5) DEVICEVAR(u, r)
   !$omp OMPLOOP collapse(5) DEVICEPTR(u, r)
   do i5=lo(5), hi(5)
      do i4=lo(4), hi(4)
         do i3=lo(3), hi(3)
            do i2=lo(2), hi(2)
               do i1=lo(1), hi(1)
                  if (kind == RES_RICCATI) then
                     r(i1,i2,i3,i4,i5) = -u(i1,i2,i3,i4,i5)**2
                  elseif (kind == RES_OSC) then
                     r(i1,i2,i3,i4,i5) = merge(u(2,i2,i3,i4,i5), -u(1,i2,i3,i4,i5), i1 == 1)
                  else
                     r(i1,i2,i3,i4,i5) = c
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine residual

   subroutine step(st, f, dt, kind, j, m, window_ok, coef_err, c0)
   !< One step of the field at its rank. window_ok: the views are null outside each stage window and right inside;
   !< coef_err (with c0, the state at the start of the step): max |state of stage k - c0 - dt stage_time(k)|.
   type(stepper), intent(inout)           :: st
   type(field),   intent(inout)           :: f
   real(R_P),     intent(in)              :: dt
   integer(I4P),  intent(in)              :: kind, j
   type(field),   intent(in),    optional :: m
   logical,       intent(inout), optional :: window_ok
   real(R_P),     intent(inout), optional :: coef_err
   real(R_P),     intent(in),    optional :: c0
   integer(I4P)                           :: k
   real(R_P), pointer                     :: u1(:), r1(:), u2(:,:), r2(:,:), u3(:,:,:), r3(:,:,:), u4(:,:,:,:), &
                                             r4(:,:,:,:), u5(:,:,:,:,:), r5(:,:,:,:,:)
   real(R_P), allocatable                 :: h(:,:,:,:,:)
   logical                                :: on_q

   call begin_step_f(st, f, m)
   if (present(window_ok)) call check_null(st, f, window_ok)
   do k=1_I4P, st%stages()
      call begin_stage_f(st, f, k, dt, m)
      on_q = st%stage_state_is_q(k)
      select case(f%rank)
      case(1)
         call st%stage_state_view(u1) ; call st%rhs_register_view(r1)
         if (present(window_ok)) window_ok = window_ok .and. associated(r1) .and. (associated(u1) .neqv. on_q)
         if (on_q) u1 => f%p1
         call residual(f, kind, k, j, u1, r1)
         if (present(coef_err)) h = host_of(f, u1)
      case(2)
         call st%stage_state_view(u2) ; call st%rhs_register_view(r2)
         if (present(window_ok)) window_ok = window_ok .and. associated(r2) .and. (associated(u2) .neqv. on_q)
         if (on_q) u2 => f%p2
         call residual(f, kind, k, j, u2, r2)
         if (present(coef_err)) h = host_of(f, u2)
      case(3)
         call st%stage_state_view(u3) ; call st%rhs_register_view(r3)
         if (present(window_ok)) window_ok = window_ok .and. associated(r3) .and. (associated(u3) .neqv. on_q)
         if (on_q) u3 => f%p3
         call residual(f, kind, k, j, u3, r3)
         if (present(coef_err)) h = host_of(f, u3)
      case(4)
         call st%stage_state_view(u4) ; call st%rhs_register_view(r4)
         if (present(window_ok)) window_ok = window_ok .and. associated(r4) .and. (associated(u4) .neqv. on_q)
         if (on_q) u4 => f%p4
         call residual(f, kind, k, j, u4, r4)
         if (present(coef_err)) h = host_of(f, u4)
      case(5)
         call st%stage_state_view(u5) ; call st%rhs_register_view(r5)
         if (present(window_ok)) window_ok = window_ok .and. associated(r5) .and. (associated(u5) .neqv. on_q)
         if (on_q) u5 => f%s
         call residual(f, kind, k, j, u5, r5)
         if (present(coef_err)) h = host_of(f, u5)
      endselect
      if (present(coef_err)) coef_err = max(coef_err, interior_max(f, abs(h - c0 - dt * st%stage_time(k))))
      call end_stage_f(st, f, k, dt, m)
      if (present(window_ok)) call check_null(st, f, window_ok)
   enddo
   call finish_step_f(st, f, dt, m)
   endsubroutine step

   function host_of(f, a) result(h)
   !< A host copy of an array with the storage's shape (any rank, by sequence association).
   type(field), intent(in) :: f
   real(R_P),   intent(in) :: a(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4), &
                                f%lb5(5):f%ub5(5))
   real(R_P), allocatable  :: h(:,:,:,:,:)
   allocate(h(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4),f%lb5(5):f%ub5(5)))
   call dev_memcpy_from_device(dst=h, src=a)
   endfunction host_of

   function interior_max(f, h) result(e)
   !< Max of h over the cells of the box.
   type(field), intent(in) :: f
   real(R_P),   intent(in) :: h(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):)
   real(R_P)               :: e
   associate(lo=>f%bx%lo, hi=>f%bx%hi)
   e = maxval(h(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4),lo(5):hi(5)))
   endassociate
   endfunction interior_max

   logical function outside_is(f, h, v) result(ok)
   !< Every element of h outside the box equals v, bit for bit.
   type(field), intent(in) :: f
   real(R_P),   intent(in) :: h(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):)
   real(R_P),   intent(in) :: v
   integer(I4P)            :: i1, i2, i3, i4, i5
   ok = .true.
   do i5=lbound(h,5), ubound(h,5) ; do i4=lbound(h,4), ubound(h,4) ; do i3=lbound(h,3), ubound(h,3)
   do i2=lbound(h,2), ubound(h,2) ; do i1=lbound(h,1), ubound(h,1)
      if (all([i1,i2,i3,i4,i5] >= f%bx%lo) .and. all([i1,i2,i3,i4,i5] <= f%bx%hi)) cycle
      if (h(i1,i2,i3,i4,i5) /= v) ok = .false.
   enddo ; enddo ; enddo ; enddo ; enddo
   endfunction outside_is

   subroutine check_null(st, f, ok)
   !< Outside a stage window both views are null.
   type(stepper), intent(in)    :: st
   type(field),   intent(in)    :: f
   logical,       intent(inout) :: ok
   real(R_P), pointer           :: v1(:), v2(:,:), v3(:,:,:), v4(:,:,:,:), v5(:,:,:,:,:)
   select case(f%rank)
   case(1) ; call st%stage_state_view(v1) ; ok = ok .and. .not. associated(v1)
             call st%rhs_register_view(v1) ; ok = ok .and. .not. associated(v1)
   case(2) ; call st%stage_state_view(v2) ; ok = ok .and. .not. associated(v2)
             call st%rhs_register_view(v2) ; ok = ok .and. .not. associated(v2)
   case(3) ; call st%stage_state_view(v3) ; ok = ok .and. .not. associated(v3)
             call st%rhs_register_view(v3) ; ok = ok .and. .not. associated(v3)
   case(4) ; call st%stage_state_view(v4) ; ok = ok .and. .not. associated(v4)
             call st%rhs_register_view(v4) ; ok = ok .and. .not. associated(v4)
   case(5) ; call st%stage_state_view(v5) ; ok = ok .and. .not. associated(v5)
             call st%rhs_register_view(v5) ; ok = ok .and. .not. associated(v5)
   endselect
   endsubroutine check_null

   ! ---- the four step procedures, at the field's rank, with or without a mask
   subroutine begin_step_f(st, f, m)
   type(stepper), intent(inout)        :: st
   type(field),   intent(in)           :: f
   type(field),   intent(in), optional :: m
   if (present(m)) then
      select case(f%rank)
      case(1) ; call st%begin_step(f%bx, f%p1, m%p1)
      case(2) ; call st%begin_step(f%bx, f%p2, m%p2)
      case(3) ; call st%begin_step(f%bx, f%p3, m%p3)
      case(4) ; call st%begin_step(f%bx, f%p4, m%p4)
      case(5) ; call st%begin_step(f%bx, f%s, m%s)
      endselect
   else
      select case(f%rank)
      case(1) ; call st%begin_step(f%bx, f%p1)
      case(2) ; call st%begin_step(f%bx, f%p2)
      case(3) ; call st%begin_step(f%bx, f%p3)
      case(4) ; call st%begin_step(f%bx, f%p4)
      case(5) ; call st%begin_step(f%bx, f%s)
      endselect
   endif
   endsubroutine begin_step_f

   subroutine begin_stage_f(st, f, k, dt, m)
   type(stepper), intent(inout)        :: st
   type(field),   intent(in)           :: f
   integer(I4P),  intent(in)           :: k
   real(R_P),     intent(in)           :: dt
   type(field),   intent(in), optional :: m
   if (present(m)) then
      select case(f%rank)
      case(1) ; call st%begin_stage(k, f%bx, dt, f%p1, m%p1)
      case(2) ; call st%begin_stage(k, f%bx, dt, f%p2, m%p2)
      case(3) ; call st%begin_stage(k, f%bx, dt, f%p3, m%p3)
      case(4) ; call st%begin_stage(k, f%bx, dt, f%p4, m%p4)
      case(5) ; call st%begin_stage(k, f%bx, dt, f%s, m%s)
      endselect
   else
      select case(f%rank)
      case(1) ; call st%begin_stage(k, f%bx, dt, f%p1)
      case(2) ; call st%begin_stage(k, f%bx, dt, f%p2)
      case(3) ; call st%begin_stage(k, f%bx, dt, f%p3)
      case(4) ; call st%begin_stage(k, f%bx, dt, f%p4)
      case(5) ; call st%begin_stage(k, f%bx, dt, f%s)
      endselect
   endif
   endsubroutine begin_stage_f

   subroutine end_stage_f(st, f, k, dt, m)
   type(stepper), intent(inout)        :: st
   type(field),   intent(inout)        :: f
   integer(I4P),  intent(in)           :: k
   real(R_P),     intent(in)           :: dt
   type(field),   intent(in), optional :: m
   if (present(m)) then
      select case(f%rank)
      case(1) ; call st%end_stage(k, f%bx, dt, f%p1, m%p1)
      case(2) ; call st%end_stage(k, f%bx, dt, f%p2, m%p2)
      case(3) ; call st%end_stage(k, f%bx, dt, f%p3, m%p3)
      case(4) ; call st%end_stage(k, f%bx, dt, f%p4, m%p4)
      case(5) ; call st%end_stage(k, f%bx, dt, f%s, m%s)
      endselect
   else
      select case(f%rank)
      case(1) ; call st%end_stage(k, f%bx, dt, f%p1)
      case(2) ; call st%end_stage(k, f%bx, dt, f%p2)
      case(3) ; call st%end_stage(k, f%bx, dt, f%p3)
      case(4) ; call st%end_stage(k, f%bx, dt, f%p4)
      case(5) ; call st%end_stage(k, f%bx, dt, f%s)
      endselect
   endif
   endsubroutine end_stage_f

   subroutine finish_step_f(st, f, dt, m)
   type(stepper), intent(inout)        :: st
   type(field),   intent(inout)        :: f
   real(R_P),     intent(in)           :: dt
   type(field),   intent(in), optional :: m
   if (present(m)) then
      select case(f%rank)
      case(1) ; call st%finish_step(f%bx, dt, f%p1, mask=m%p1)
      case(2) ; call st%finish_step(f%bx, dt, f%p2, mask=m%p2)
      case(3) ; call st%finish_step(f%bx, dt, f%p3, mask=m%p3)
      case(4) ; call st%finish_step(f%bx, dt, f%p4, mask=m%p4)
      case(5) ; call st%finish_step(f%bx, dt, f%s, mask=m%s)
      endselect
   else
      select case(f%rank)
      case(1) ; call st%finish_step(f%bx, dt, f%p1)
      case(2) ; call st%finish_step(f%bx, dt, f%p2)
      case(3) ; call st%finish_step(f%bx, dt, f%p3)
      case(4) ; call st%finish_step(f%bx, dt, f%p4)
      case(5) ; call st%finish_step(f%bx, dt, f%s)
      endselect
   endif
   endsubroutine finish_step_f
endmodule selah_test_ranks_support

program selah_test_stepper_ranks
!< SELAH test: the stepper on fields of rank 1 to 4, against the same data as rank 5, on the host or the device.
!<
!< Seven shapes: scalar fields of rank 1 to 4 (`(i)`, `(i, j)`, `(i, j, k)`, `(i, j, k, b)`, the last with a box from
!< `box_range`) and fields with two components of rank 2 to 4 (`(v, i)`, `(v, i, j)`, `(v, i, j, k)`), one ghost
!< cell each side of every cell dimension. The rank-k field is a pointer of its rank over contiguous storage, as a
!< caller holds it. For each shape, every Runge–Kutta scheme in every form:
!<
!<+ **Rank k = rank 5.** Three steps of y' = -y^2 give bitwise the values of the same data stepped as a rank-5 field.
!<+ **Order.** y' = -y^2 cell by cell (exact 1/(1/y0 + t)) reaches the scheme's order (R8P builds).
!<+ **Contract.** The ghost cells of q and of every register keep their sentinel; both views are null outside every
!<  stage window, the residual view is associated inside it, the state view exactly when the stage is not on q.
!<+ **Stage coefficients.** With the residual 1 the state of stage k is q^n + dt stage_time(k); with the indicator of
!<  stage j the step adds dt stage_weight(j).
!<+ **Masks** (the shapes with components): a mask of q's rank with three components, half of the cells solid. The
!<  solid cells keep their values bit for bit, the fluid ones are bitwise those of the unmasked run.
!<+ **Splittings** (the shapes with components): the oscillator q' = p, p' = -q cell by cell reaches each splitting's
!<  order (R8P builds).
use selah_test_ranks_support
use selah_test_tools, only : check, failures
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
implicit none
integer(I4P), parameter :: NI = 3_I4P, NV = 2_I4P              !< Interior cells per dimension, components.
real(R_P),    parameter :: SENTINEL = -777._R_P                 !< Ghost value nobody must touch.
character(9), parameter :: RK_FORMS(5) = [character(9) :: 'butcher', 'butcher-q', 'two-reg', '2n', 'ket']
integer(I4P)            :: c, s, fm

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
do c=1_I4P, 7_I4P
   do s=1_I4P, size(RK_SCHEMES, kind=I4P)
      do fm=1_I4P, size(RK_FORMS, kind=I4P)
         if (registers_needed(trim(RK_SCHEMES(s)), form_of(trim(RK_FORMS(fm)))) < 0) cycle
         call rk_case(c, trim(RK_SCHEMES(s)), trim(RK_FORMS(fm)))
      enddo
   enddo
   if (shape_comps(c)) then
      do s=1_I4P, size(SPLITTING_SCHEMES, kind=I4P)
         call split_case(c, trim(SPLITTING_SCHEMES(s)))
      enddo
   endif
enddo
if (failures > 0) error stop 'selah_test_stepper_ranks: failed'
print '(A)', 'selah_test_stepper_ranks: ok'

contains
   function form_of(tag) result(form)
   character(*), intent(in)  :: tag
   character(:), allocatable :: form
   select case(tag)
   case('butcher', 'butcher-q') ; form = FORM_BUTCHER
   case('two-reg') ; form = FORM_SHU_OSHER_2R
   case('2n')      ; form = FORM_WILLIAMSON_2N
   case default    ; form = FORM_KETCHESON_2R
   endselect
   endfunction form_of

   logical function shape_comps(c)
   integer(I4P), intent(in) :: c
   shape_comps = c >= 5_I4P
   endfunction shape_comps

   function shape_tag(c) result(tag)
   integer(I4P), intent(in)  :: c
   character(:), allocatable :: tag
   character(16), parameter  :: TAGS(7) = [character(16) :: '(i)', '(i,j)', '(i,j,k)', '(i,j,k,b) range', '(v,i)', &
                                           '(v,i,j)', '(v,i,j,k)']
   tag = trim(TAGS(c))
   endfunction shape_tag

   subroutine make_shape(f, c, rank5, ncomp)
   !< Field of shape c (1-4 scalar, 5-7 with components), or its rank-5 twin; ncomp: a mask's components.
   type(field),  intent(inout)        :: f
   integer(I4P), intent(in)           :: c
   logical,      intent(in)           :: rank5
   integer(I4P), intent(in), optional :: ncomp
   integer(I4P)                       :: k, lb(4), ub(4), lo(4), hi(4), nc
   logical                            :: comps
   type(box)                          :: bx

   comps = shape_comps(c)
   k = merge(c - 3_I4P, c, comps)
   lb = 0_I4P ; ub = NI + 1_I4P ; lo = 1_I4P ; hi = NI
   if (comps) then
      nc = NV ; if (present(ncomp)) nc = ncomp
      lb(1) = 1_I4P ; ub(1) = nc ; lo(1) = 1_I4P ; hi(1) = nc
   endif
   if (c == 4_I4P) then ! (i, j, k, b): two blocks, no ghosts on the blocks (box_range)
      lb(4) = 1_I4P ; ub(4) = 2_I4P ; lo(4) = 1_I4P ; hi(4) = 2_I4P
   endif
   if (rank5) then ! the same storage as rank 5: the mapped bounds and range
      bx = box_range(lb(1:k), lo(1:k), hi(1:k), comps)
      call make(f, 5_I4P, .true., bx%lb, front(bx, ub(1:k)), bx%lo, bx%hi)
   else
      call make(f, k, comps, lb(1:k), ub(1:k), lo(1:k), hi(1:k))
      if (c /= 4_I4P .and. .not. present(ncomp)) then ! the interior constructor gives the same box
         bx = box_interior(lb(1:k), ub(1:k), 1_I4P, comps)
         if (any(bx%lo /= f%bx%lo) .or. any(bx%hi /= f%bx%hi) .or. any(bx%lb /= f%bx%lb) .or. bx%rank /= k) &
            call check(.false., shape_tag(c)//': box_interior = box_range of the interior')
      endif
   endif
   endsubroutine make_shape

   function front(bx, ub) result(ub5)
   type(box),    intent(in) :: bx
   integer(I4P), intent(in) :: ub(:)
   integer(I4P)             :: ub5(5), d
   ub5 = 1_I4P
   do d=1_I4P, size(ub, kind=I4P)
      ub5(box_dim(bx, d)) = ub(d)
   enddo
   endfunction front

   subroutine fill(f, h)
   !< Interior values 1 + 0.01 n (n the cell's position), the sentinel elsewhere.
   type(field),            intent(inout) :: f
   real(R_P), allocatable, intent(out)   :: h(:,:,:,:,:)
   integer(I4P)                          :: i1, i2, i3, i4, i5, n
   allocate(h(f%lb5(1):f%ub5(1),f%lb5(2):f%ub5(2),f%lb5(3):f%ub5(3),f%lb5(4):f%ub5(4),f%lb5(5):f%ub5(5)))
   h = SENTINEL ; n = 0_I4P
   do i5=f%bx%lo(5), f%bx%hi(5) ; do i4=f%bx%lo(4), f%bx%hi(4) ; do i3=f%bx%lo(3), f%bx%hi(3)
   do i2=f%bx%lo(2), f%bx%hi(2) ; do i1=f%bx%lo(1), f%bx%hi(1)
      n = n + 1_I4P ; h(i1,i2,i3,i4,i5) = 1._R_P + 0.01_R_P * n
   enddo ; enddo ; enddo ; enddo ; enddo
   call set_host(f, h)
   endsubroutine fill

   subroutine sentinel_registers(st)
   type(stepper), intent(inout) :: st
   real(R_P), allocatable       :: h(:,:,:,:,:,:)
   allocate(h, mold=st%regs) ; h = SENTINEL
   call dev_memcpy_to_device(dst=st%regs, src=h)
   endsubroutine sentinel_registers

   logical function registers_untouched(st, f) result(ok)
   type(stepper), intent(in) :: st
   type(field),   intent(in) :: f
   real(R_P), allocatable    :: h(:,:,:,:,:,:)
   integer(I4P)              :: k
   allocate(h(lbound(st%regs,1):ubound(st%regs,1), lbound(st%regs,2):ubound(st%regs,2), &
              lbound(st%regs,3):ubound(st%regs,3), lbound(st%regs,4):ubound(st%regs,4), &
              lbound(st%regs,5):ubound(st%regs,5), lbound(st%regs,6):ubound(st%regs,6)))
   call dev_memcpy_from_device(dst=h, src=st%regs)
   ok = .true.
   do k=lbound(h,6), ubound(h,6)
      ok = ok .and. outside_is(f, h(:,:,:,:,:,k), SENTINEL)
   enddo
   endfunction registers_untouched

   subroutine rk_case(c, name, ftag)
   !< Every check of one Runge–Kutta scheme in one form on shape c.
   integer(I4P), intent(in) :: c
   character(*), intent(in) :: name, ftag
   type(field)              :: f, g, m
   type(stepper)            :: st, st5
   real(R_P), allocatable   :: h(:,:,:,:,:), h5(:,:,:,:,:), h0(:,:,:,:,:), hm(:,:,:,:,:)
   real(R_P)                :: err(4), dt, ce, tol
   integer(I4P)             :: i, n, step_, j
   logical                  :: ok, window_ok, on_q
   type(butcher_tableau)    :: bt
   character(:), allocatable :: tag

   tag = shape_tag(c)//' '//name//' '//ftag
   on_q = ftag == 'butcher-q'
   bt = butcher(name)
   ! rank k = rank 5, bitwise; contract (sentinels, windows)
   call make_shape(f, c, .false.) ; call make_shape(g, c, .true.)
   call fill(f, h0) ; call set_host(g, h0)
   call init_stepper(st, f, name, form_of(ftag)) ; call init_stepper(st5, g, name, form_of(ftag))
   if (on_q) then
      call st%destroy ; call st5%destroy
      call init_stepper_q(st, f, name) ; call init_stepper_q(st5, g, name)
   endif
   call sentinel_registers(st)
   window_ok = .true.
   do step_=1_I4P, 3_I4P
      call step(st, f, 0.1_R_P, RES_RICCATI, 0_I4P, window_ok=window_ok)
      call step(st5, g, 0.1_R_P, RES_RICCATI, 0_I4P)
   enddo
   h = get_host(f) ; h5 = get_host(g)
   call check(all(h == h5), tag//': rank k bitwise rank 5')
   call check(outside_is(f, h, SENTINEL) .and. registers_untouched(st, f), tag//': ghost cells untouched (q, registers)')
   call check(window_ok, tag//': views null outside the stage windows, associated inside')
   ! stage coefficients: residual 1 (stage_time), indicator of stage j (stage_weight)
   dt = 0.1_R_P ; tol = 64._R_P * epsilon(1._R_P) * dt
   h = h0 ; h = 0._R_P ; call set_host(f, h)
   ce = 0._R_P
   call step(st, f, dt, RES_ONE, 0_I4P, coef_err=ce, c0=0._R_P)
   ok = ce <= tol
   do j=1_I4P, st%stages()
      call set_host(f, h)
      call step(st, f, dt, RES_STAGE, j)
      ok = ok .and. interior_max(f, abs(get_host(f) - dt * st%stage_weight(j))) <= tol
   enddo
   call check(ok .and. ce <= tol, tag//': stage k is stage k of the tableau (stage_time, stage_weight)')
   call st%destroy ; call st5%destroy
   ! masks: the shapes with components
   if (shape_comps(c)) then
      call make_shape(m, c, .false., ncomp=3_I4P)
      hm = get_host(m) ; hm = 5._R_P                     ! the mask's components 1, 2 carry nothing
      hm(3,:,:,:,:) = 1._R_P                               ! component 3 (the last, the default): solid
      hm(3,lbound(hm,2):lbound(hm,2)+1,:,:,:) = -1._R_P    ! fluid: the ghost and first layer of the first cell index
      call set_host(m, hm)
      call fill(f, h0) ; call set_host(g, h0)
      call init_stepper(st, f, name, form_of(ftag)) ; call init_stepper(st5, g, name, form_of(ftag))
      do step_=1_I4P, 2_I4P
         call step(st, f, 0.1_R_P, RES_RICCATI, 0_I4P, m=m)
         call step(st5, g, 0.1_R_P, RES_RICCATI, 0_I4P)
      enddo
      h = get_host(f) ; h5 = get_host(g)
      call check(mask_ok(f, h, h5, h0, hm), tag//': mask: solid cells untouched, fluid cells bitwise the unmasked run')
      call st%destroy ; call st5%destroy ; call unmake(m)
   endif
   call unmake(f) ; call unmake(g)
#if !defined _R_P_IS_R4P
   ! order on y' = -y^2, cell by cell
   n = merge(8_I4P, 16_I4P, bt%order >= 4)
   do i=1_I4P, 4_I4P
      call make_shape(f, c, .false.) ; call fill(f, h0)
      if (on_q) then ; call init_stepper_q(st, f, name) ; else ; call init_stepper(st, f, name, form_of(ftag)) ; endif
      dt = 0.5_R_P / (n * 2_I4P**(i-1_I4P))
      do step_=1_I4P, n * 2_I4P**(i-1_I4P)
         call step(st, f, dt, RES_RICCATI, 0_I4P)
      enddo
      h = get_host(f)
      err(i) = interior_max(f, abs(h - 1._R_P / (1._R_P / h0 + 0.5_R_P)))
      call st%destroy ; call unmake(f)
   enddo
   err(4) = max(err(4), tiny(1._R_P))
   call check(log(err(3) / err(4)) / log(2._R_P) >= real(bt%order, R_P) - 0.2_R_P, tag//': observed order')
#endif
   endsubroutine rk_case

   logical function mask_ok(f, h, h5, h0, hm) result(ok)
   !< In the cells of the box: fluid (hm(3, cell) < 0) bitwise the unmasked run h5, solid bitwise the start h0.
   type(field), intent(in) :: f
   real(R_P),   intent(in) :: h(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):)  !< Masked run.
   real(R_P),   intent(in) :: h5(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):) !< Unmasked run.
   real(R_P),   intent(in) :: h0(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):) !< Start.
   real(R_P),   intent(in) :: hm(f%lb5(1):,f%lb5(2):,f%lb5(3):,f%lb5(4):,f%lb5(5):) !< Mask.
   integer(I4P)            :: i1, i2, i3, i4, i5
   logical                 :: some_fluid, some_solid
   ok = .true. ; some_fluid = .false. ; some_solid = .false.
   do i5=f%bx%lo(5), f%bx%hi(5) ; do i4=f%bx%lo(4), f%bx%hi(4) ; do i3=f%bx%lo(3), f%bx%hi(3)
   do i2=f%bx%lo(2), f%bx%hi(2) ; do i1=f%bx%lo(1), f%bx%hi(1)
      if (hm(3,i2,i3,i4,i5) < 0._R_P) then
         ok = ok .and. h(i1,i2,i3,i4,i5) == h5(i1,i2,i3,i4,i5) ; some_fluid = .true.
      else
         ok = ok .and. h(i1,i2,i3,i4,i5) == h0(i1,i2,i3,i4,i5) ; some_solid = .true.
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   ok = ok .and. some_fluid .and. some_solid
   endfunction mask_ok

   subroutine init_stepper_q(st, f, name)
   !< The Butcher form with stages_on_q.
   type(stepper), intent(inout) :: st
   type(field),   intent(in)    :: f
   character(*),  intent(in)    :: name
   select case(f%rank)
   case(1) ; call st%initialize(name, FORM_BUTCHER, lbound(f%p1), ubound(f%p1), stages_on_q=.true., components=f%comps)
   case(2) ; call st%initialize(name, FORM_BUTCHER, lbound(f%p2), ubound(f%p2), stages_on_q=.true., components=f%comps)
   case(3) ; call st%initialize(name, FORM_BUTCHER, lbound(f%p3), ubound(f%p3), stages_on_q=.true., components=f%comps)
   case(4) ; call st%initialize(name, FORM_BUTCHER, lbound(f%p4), ubound(f%p4), stages_on_q=.true., components=f%comps)
   case(5) ; call st%initialize(name, FORM_BUTCHER, f%lb5, f%ub5, stages_on_q=.true.)
   endselect
   endsubroutine init_stepper_q

   subroutine split_case(c, name)
   !< A splitting on the oscillator q' = p, p' = -q cell by cell (shape c has two components): its order.
   integer(I4P), intent(in) :: c
   character(*), intent(in) :: name
   type(field)              :: f
   type(stepper)            :: st
   type(splitting_composition) :: sc
   real(R_P), allocatable   :: h(:,:,:,:,:), h0(:,:,:,:,:), hx(:,:,:,:,:)
   real(R_P)                :: err(3), dt
   integer(I4P)             :: i, n, step_, order
   logical                  :: ok

   sc = splitting(name)
   order = sc%order
   ok = .true.
   do i=1_I4P, 3_I4P
      call make_shape(f, c, .false.) ; call fill(f, h0)
      call init_stepper(st, f, name, FORM_SPLITTING, partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
      n = 8_I4P * 2_I4P**(i-1_I4P) ; dt = 1._R_P / n
      do step_=1_I4P, n
         call step(st, f, dt, RES_OSC, 0_I4P)
      enddo
      h = get_host(f)
      hx = h0 ! the exact rotation by t = 1
      hx(1,:,:,:,:) = h0(1,:,:,:,:) * cos(1._R_P) + h0(2,:,:,:,:) * sin(1._R_P)
      hx(2,:,:,:,:) = -h0(1,:,:,:,:) * sin(1._R_P) + h0(2,:,:,:,:) * cos(1._R_P)
      err(i) = interior_max(f, abs(h - hx))
      call st%destroy ; call unmake(f)
   enddo
#if !defined _R_P_IS_R4P
   err(3) = max(err(3), tiny(1._R_P))
   call check(log(err(2) / err(3)) / log(2._R_P) >= real(order, R_P) - 0.2_R_P .or. err(3) < 1e-13_R_P, &
              shape_tag(c)//' '//name//': observed order')
#else
   call check(all(err < 1e-2_R_P), shape_tag(c)//' '//name//': converges')
#endif
   endsubroutine split_case
endprogram selah_test_stepper_ranks
