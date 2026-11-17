program box_tour
!< Tutorial, chapter 3: a block-structured field, the box a kernel loops over, and what it leaves untouched.
use selah
implicit none
integer(I4P), parameter :: NV=2, NI=4, NJ=3, NK=2, NGC=2, NB=4   ! components, cells per block, ghosts, pool of blocks
real(R_P),    parameter :: SENTINEL=-9._R_P
real(R_P), allocatable  :: x(:,:,:,:,:), y(:,:,:,:,:)
type(box)               :: bx, part

! V_FIRST layout (v,i,j,k,b): the components innermost (the CPU layout); ghosts 1-NGC..N+NGC on i, j, k
allocate(x(1:NV, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NB))
allocate(y, mold=x)
x = 1._R_P
bx = box_interior(LAYOUT_V_FIRST, NV, NI, NJ, NK, NGC, 3_I4P)   ! the interior of the first 3 blocks of the pool
print '(A,5I4)', 'lb (allocation lower bounds):', bx%lb
print '(A,5I4)', 'lo (loop start):             ', bx%lo
print '(A,5I4)', 'hi (loop end):               ', bx%hi
y = SENTINEL
call selah_axpy(bx, 2._R_P, x, y)                                 ! y = y + 2 x, in the box only
print '(A,I0,A,I0)', 'cells updated: ', count(y /= SENTINEL), ' of ', size(y)
print '(A,I0)', 'expected (2 components x 4x3x2 cells x 3 blocks): ', NV*NI*NJ*NK*3
print '(A,L1)', 'ghosts and block 4 untouched: ', all(y(:,:,:,:,4) == SENTINEL) .and. all(y(:,0,:,:,:) == SENTINEL)
part = box_components(bx, 2_I4P, 2_I4P)                           ! component 2 only: a partition
y = SENTINEL
call selah_copy(part, x, y)
print '(A,I0,A,I0)', 'component range 2..2: cells updated ', count(y /= SENTINEL), ', component 1 untouched: ', &
      merge(1, 0, all(y(1,:,:,:,:) == SENTINEL))
! V_LAST layout (b,i,j,k,v): the blocks innermost (ADAM's device layout); the same box, built for that layout
deallocate(x, y)
allocate(x(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NV))
allocate(y, mold=x)
x = 1._R_P ; y = SENTINEL
bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, 3_I4P)
call selah_axpy(bx, 2._R_P, x, y)
print '(A,I0,A,I0)', 'V_LAST: component dimension ', component_dim(bx), ', cells updated ', count(y /= SENTINEL)
endprogram box_tour
