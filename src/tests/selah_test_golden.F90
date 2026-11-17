!< SELAH regression C4: the values themselves, pinned in a golden file written once and never edited.
program selah_test_golden
!< SELAH regression C4: the values themselves, pinned in a golden file written once and never edited.
!<
!< The verification cases assert orders, properties and bounds; a change of a few ulp, or an error growing inside its
!< bound, passes them. This test pins the values: every Runge–Kutta scheme in every storage form it has (3 steps of a
!< linear system), every splitting (3 steps of a partitioned linear system), and every kernel (box, both layouts,
!< masked and unmasked; flat), on inputs exact in binary (dyadic values and coefficients: no libm, no random numbers),
!< so that the inputs are the same bits on every compiler and option.
!<
!< Usage: selah_test_golden            compare with the golden file of the build's precision
!<        selah_test_golden write FILE write a new golden file (refuses to overwrite an existing one)
!< The golden file is src/tests/golden/selah_golden_v1.txt (R8P) or selah_golden_r4_v1.txt (R4P), looked up from the
!< working directory, then from $SELAH_SRC. Tolerance: 1000 eps(R_P), relative to max(|value|, 1). A deliberate change
!< of results gets a new file (v2, ...) and a new default path here; an existing file is never edited.
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use selah
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter   :: N=32
real(R_P),    parameter   :: TOL = 1000._R_P * epsilon(1._R_P)
character(*), parameter   :: GOLDEN_R8 = 'src/tests/golden/selah_golden_v1.txt'
character(*), parameter   :: GOLDEN_R4 = 'src/tests/golden/selah_golden_r4_v1.txt'
integer,      parameter   :: MAXREC=200
character(64)             :: labels(MAXREC) ! fixed capacity: growing a character array by concatenation lost every
real(R_P)                 :: values(N,MAXREC) ! element but the last with nvfortran 26.9
integer                   :: nrec = 0
character(512)            :: arg, path
integer                   :: nlines

call compute_all
if (command_argument_count() >= 2) then
   call get_command_argument(1, arg) ; call get_command_argument(2, path)
   if (trim(arg) /= 'write') error stop 'usage: selah_test_golden [write FILE]'
   call write_golden(trim(path))
else
   if (R_P == R8P) then
      call compare_golden(locate(GOLDEN_R8))
   else
      call compare_golden(locate(GOLDEN_R4))
   endif
endif
if (failures > 0) error stop 'selah_test_golden: failed'
print '(A)', 'selah_test_golden: ok'

contains
   subroutine record(label, v)
   !< Append one line of N values.
   character(*), intent(in) :: label
   real(R_P),    intent(in) :: v(:)

   if (nrec == MAXREC) error stop 'selah_test_golden: raise MAXREC'
   nrec = nrec + 1
   labels(nrec) = label
   values(:, nrec) = 0._R_P
   values(1:size(v), nrec) = v
   endsubroutine record

   pure function dyadic(i, seed) result(v)
   !< A value exact in binary, in [-1, 1]: (((i * 37 + seed * 11) mod 61) - 30) / 32.
   integer, intent(in) :: i, seed
   real(R_P)           :: v
   v = real(modulo(i * 37 + seed * 11, 61) - 30, R_P) / 32._R_P
   endfunction dyadic

   subroutine compute_all
   !< Every scheme in every form, every splitting, every kernel.
   character(14), parameter :: FORMS(4) = [character(14) :: FORM_BUTCHER, FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, &
                                                             FORM_KETCHESON_2R]
   integer :: s, f
   do s=1, size(RK_SCHEMES)
      do f=1, size(FORMS)
         if (registers_needed(trim(RK_SCHEMES(s)), trim(FORMS(f))) < 0) cycle
         call record('rk '//trim(RK_SCHEMES(s))//' '//trim(FORMS(f)), rk_run(trim(RK_SCHEMES(s)), trim(FORMS(f))))
      enddo
   enddo
   do s=1, size(SPLITTING_SCHEMES)
      call record('split '//trim(SPLITTING_SCHEMES(s)), split_run(trim(SPLITTING_SCHEMES(s))))
   enddo
   call kernels
   endsubroutine compute_all

   function rk_run(name, form) result(v)
   !< 3 steps of y' = L y, L periodic tridiagonal (1/4, -1/2 + 1/8, 1/4), dt = 1/16, y0 dyadic.
   character(*), intent(in)       :: name, form
   real(R_P)                      :: v(N)
   real(R_P), allocatable, target :: y(:,:,:,:,:)
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   integer(I4P)                   :: step, k, i

   allocate(y(1:1, 0:N+1, 0:2, 0:2, 1:1)) ; y = 0._R_P
   y(1,1:N,1,1,1) = [(dyadic(i, 1), i=1, N)]
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, N, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(name, form, lb=lbound(y), ub=ubound(y))
   do step=1, 3
      call st%begin_step(bx, y)
      do k=1, st%stages()
         call st%begin_stage(k, bx, 0.0625_R_P, y)
         if (st%stage_state_is_q(k)) then ; u => y ; else ; u => st%stage_state() ; endif
         r => st%rhs_register()
         u(1,0,1,1,1) = u(1,N,1,1,1) ; u(1,N+1,1,1,1) = u(1,1,1,1,1)
         do i=1, N
            r(1,i,1,1,1) = 0.25_R_P * (u(1,i-1,1,1,1) + u(1,i+1,1,1,1)) - 0.375_R_P * u(1,i,1,1,1)
         enddo
         call st%end_stage(k, bx, 0.0625_R_P, y)
      enddo
      call st%finish_step(bx, 0.0625_R_P, y)
   enddo
   v = y(1,1:N,1,1,1)
   call st%destroy
   endfunction rk_run

   function split_run(name) result(v)
   !< 3 steps, dt = 1/8, of q' = p/2, p' = -q/4 on N/2 cells (q, p dyadic); v = (q, p).
   character(*), intent(in)       :: name
   real(R_P)                      :: v(N)
   real(R_P), allocatable, target :: y(:,:,:,:,:)
   real(R_P), pointer             :: r(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   integer(I4P)                   :: step, k, i

   allocate(y(1:2, 0:N/2+1, 0:2, 0:2, 1:1)) ; y = 0._R_P
   y(1,1:N/2,1,1,1) = [(dyadic(i, 2), i=1, N/2)] ; y(2,1:N/2,1,1,1) = [(dyadic(i, 3), i=1, N/2)]
   bx = box_interior(LAYOUT_V_FIRST, 2_I4P, N/2, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   call st%initialize(name, FORM_SPLITTING, lb=lbound(y), ub=ubound(y), partitions=reshape([1, 1, 2, 2], [2, 2]))
   do step=1, 3
      call st%begin_step(bx, y)
      do k=1, st%stages()
         call st%begin_stage(k, bx, 0.125_R_P, y)
         r => st%rhs_register()
         if (st%stage_partition(k) == 1) then
            r(1,1:N/2,1,1,1) = 0.5_R_P * y(2,1:N/2,1,1,1)
         else
            r(2,1:N/2,1,1,1) = -0.25_R_P * y(1,1:N/2,1,1,1)
         endif
         call st%end_stage(k, bx, 0.125_R_P, y)
      enddo
      call st%finish_step(bx, 0.125_R_P, y)
   enddo
   v = [y(1,1:N/2,1,1,1), y(2,1:N/2,1,1,1)]
   call st%destroy
   endfunction split_run

   subroutine kernels
   !< Every kernel on dyadic inputs: box (V_FIRST unmasked, V_LAST masked: 2 components x 4 cells x 4 blocks = 32
   !< values) and flat (32 values).
   real(R_P), parameter   :: COEF(3) = [0.25_R_P, -0.5_R_P, 0.75_R_P], DT=0.125_R_P
   real(R_P), allocatable :: x(:,:,:,:,:), y(:,:,:,:,:), q(:,:,:,:,:), r(:,:,:,:,:), m(:,:,:,:,:), regs(:,:,:,:,:,:)
   real(R_P)              :: fx(N), fy(N), fq(N), fr(N), fregs(N,3), fd(N), fm(N)
   type(box)              :: bx
   integer                :: ly, i, k
   character(:), allocatable :: t

   do i=1, N
      fx(i) = dyadic(i, 4) ; fy(i) = dyadic(i, 5) ; fq(i) = dyadic(i, 6) ; fr(i) = dyadic(i, 7)
      fm(i) = merge(1._R_P, -1._R_P, mod(i, 3) == 0)
      do k=1, 3
         fregs(i,k) = dyadic(i, 7 + k)
      enddo
   enddo
   do ly=LAYOUT_V_FIRST, LAYOUT_V_LAST
      bx = box_interior(ly, 2_I4P, 4_I4P, 1_I4P, 1_I4P, 1_I4P, 4_I4P)
      call field(x, fx, ly) ; call field(y, fy, ly) ; call field(q, fq, ly) ; call field(r, fr, ly)
      call field(m, fm, ly)
      allocate(regs(lbound(x,1):ubound(x,1), lbound(x,2):ubound(x,2), lbound(x,3):ubound(x,3), lbound(x,4):ubound(x,4), &
                    lbound(x,5):ubound(x,5), 3))
      do k=1, 3
         call field(q, fregs(:,k), ly)
         regs(:,:,:,:,:,k) = q
      enddo
      call field(q, fq, ly)
      if (ly == LAYOUT_V_FIRST) then
         t = 'box v-first '
         call selah_copy(bx, x, y)
         call record(t//'copy', cells(y, ly))
         call field(y, fy, ly)
         call selah_axpy(bx, 0.75_R_P, x, y)
         call record(t//'axpy', cells(y, ly))
         call field(y, fy, ly)
         call selah_lincomb(bx, 0.5_R_P, x, -0.25_R_P, y)
         call record(t//'lincomb', cells(y, ly))
         call field(y, fy, ly)
         call selah_stage_state(bx, 3_I4P, DT, COEF, q, regs, y)
         call record(t//'stage_state', cells(y, ly))
         call field(y, fy, ly)
         call selah_step_update(bx, 3_I4P, DT, COEF, regs, y, r)
         call record(t//'step_update q', cells(y, ly))
         call record(t//'step_update dq', cells(r, ly))
         call field(y, fy, ly) ; call field(r, fr, ly)
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y)
         call record(t//'shu_osher_2r_stage', cells(y, ly))
         call field(y, fy, ly)
         call selah_williamson_2n_stage(bx, DT, -0.5_R_P, 0.75_R_P, r, x, y)
         call record(t//'williamson_2n_stage q', cells(y, ly))
         call record(t//'williamson_2n_stage dq', cells(x, ly))
      else
         t = 'box v-last masked '
         call selah_copy(bx, x, y, m)
         call record(t//'copy', cells(y, ly))
         call field(y, fy, ly)
         call selah_axpy(bx, 0.75_R_P, x, y, m)
         call record(t//'axpy', cells(y, ly))
         call field(y, fy, ly)
         call selah_lincomb(bx, 0.5_R_P, x, -0.25_R_P, y, m)
         call record(t//'lincomb', cells(y, ly))
         call field(y, fy, ly)
         call selah_stage_state(bx, 3_I4P, DT, COEF, q, regs, y, m)
         call record(t//'stage_state', cells(y, ly))
         call field(y, fy, ly) ; call field(r, fr, ly)
         call selah_step_update(bx, 3_I4P, DT, COEF, regs, y, r, m)
         call record(t//'step_update q', cells(y, ly))
         call record(t//'step_update dq', cells(r, ly))
         call field(y, fy, ly) ; call field(r, fr, ly)
         call selah_shu_osher_2r_stage(bx, DT, 0.75_R_P, 0.25_R_P, 0.25_R_P, q, r, y, m)
         call record(t//'shu_osher_2r_stage', cells(y, ly))
         call field(y, fy, ly)
         call selah_williamson_2n_stage(bx, DT, 0._R_P, 0.75_R_P, r, x, y, m)
         call record(t//'williamson_2n_stage a=0 q', cells(y, ly))
      endif
      deallocate(regs, x, y, q, r, m)
   enddo
   fd = fy ; call selah_copy(fx, fd) ; call record('flat copy', fd)
   fd = fy ; call selah_axpy(0.75_R_P, fx, fd) ; call record('flat axpy', fd)
   fd = fy ; call selah_lincomb(0.5_R_P, fx, -0.25_R_P, fd) ; call record('flat lincomb', fd)
   fd = fy ; call selah_stage_state(3_I4P, DT, COEF, fq, fregs, fd) ; call record('flat stage_state', fd)
   fd = fy ; fx = fr ; call selah_step_update(3_I4P, DT, COEF, fregs, fd, fx) ; call record('flat step_update q', fd)
   call record('flat step_update dq', fx)
   endsubroutine kernels

   subroutine field(a, v, layout)
   !< The 32 values v into the box's cells (component, cell, block order), 0 elsewhere.
   real(R_P), allocatable, intent(inout) :: a(:,:,:,:,:)
   real(R_P),              intent(in)    :: v(:)
   integer,                intent(in)    :: layout
   integer                               :: c, i, b, n
   if (.not. allocated(a)) then
      if (layout == LAYOUT_V_FIRST) then
         allocate(a(1:2, 0:5, 0:2, 0:2, 1:4))
      else
         allocate(a(1:4, 0:5, 0:2, 0:2, 1:2))
      endif
   endif
   a = 0._R_P ; n = 0
   do b=1, 4 ; do i=1, 4 ; do c=1, 2
      n = n + 1
      if (layout == LAYOUT_V_FIRST) then ; a(c,i,1,1,b) = v(n) ; else ; a(b,i,1,1,c) = v(n) ; endif
   enddo ; enddo ; enddo
   endsubroutine field

   function cells(a, layout) result(v)
   real(R_P), intent(in) :: a(:,0:,0:,0:,:)
   integer,   intent(in) :: layout
   real(R_P)             :: v(N)
   integer               :: c, i, b, n
   n = 0
   do b=1, 4 ; do i=1, 4 ; do c=1, 2
      n = n + 1
      if (layout == LAYOUT_V_FIRST) then ; v(n) = a(c,i,1,1,b) ; else ; v(n) = a(b,i,1,1,c) ; endif
   enddo ; enddo ; enddo
   endfunction cells

   function locate(rel) result(p)
   !< The golden file: rel from the working directory, else from $SELAH_SRC.
   character(*), intent(in)  :: rel
   character(:), allocatable :: p
   character(512)            :: src
   integer                   :: l, st
   logical                   :: there

   p = rel
   inquire(file=p, exist=there)
   if (there) return
   call get_environment_variable('SELAH_SRC', src, l, st)
   if (st == 0 .and. l > 0) p = trim(src)//'/'//rel
   endfunction locate

   subroutine write_golden(file)
   character(*), intent(in) :: file
   integer                  :: u, l, ios

   open(newunit=u, file=file, status='new', action='write', iostat=ios)
   if (ios /= 0) error stop 'selah_test_golden: the golden file exists (never overwrite it) or cannot be created'
   write(u, '(A)') '# SELAH golden values v1: '//compiler_version()
   write(u, '(A)') '# options: '//compiler_options()
   write(u, '(A,I0,A)') '# one line per case: label, then ', N, ' values (ES26.17)'
   do l=1, nrec
      write(u, '(A64,*(ES26.17))') labels(l), values(:,l)
   enddo
   close(u)
   print '(A,I0,A)', 'selah_test_golden: wrote ', nrec, ' lines to '//file
   endsubroutine write_golden

   subroutine compare_golden(file)
   character(*), intent(in) :: file
   character(64)            :: lab
   character(4096)          :: line
   real(R_P)                :: ref(N), worst
   integer                  :: u, ios, l, nbit, nval, nmiss
   logical                  :: there

   inquire(file=file, exist=there)
   if (.not. there) error stop 'selah_test_golden: golden file not found (run from the source root, or set SELAH_SRC)'
   open(newunit=u, file=file, status='old', action='read')
   nlines = 0 ; nbit = 0 ; nval = 0 ; nmiss = 0 ; worst = 0._R_P
   do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      lab = line(1:64)
      read(line(65:), *) ref
      nlines = nlines + 1
      l = 0
      do ios=1, nrec ! not findloc on a character array (portability)
         if (labels(ios) == lab) then
            l = ios ; exit
         endif
      enddo
      ios = 0
      if (l == 0) then
         nmiss = nmiss + 1 ; print '(A)', '  FAIL  golden case no longer computed: '//trim(lab) ; cycle
      endif
      nval = nval + N
      nbit = nbit + count(values(:,l) == ref)
      worst = max(worst, maxval(abs(values(:,l) - ref) / max(abs(ref), 1._R_P)))
      call check(all(abs(values(:,l) - ref) <= TOL * max(abs(ref), 1._R_P)), 'golden '//trim(lab))
   enddo
   close(u)
   call check(nlines == nrec .and. nmiss == 0, 'every computed case is in the golden file, and the reverse')
   print '(A,I0,A,I0,A,ES9.2)', 'golden: ', nbit, ' of ', nval, ' values bitwise equal; worst relative difference ', worst
   endsubroutine compare_golden
endprogram selah_test_golden
