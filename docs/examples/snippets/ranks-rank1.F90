! a scalar 1-D field q(i), one ghost cell each side: the box of its interior, viewed as rank 5 (1, i, 1, 1, 1)
allocate(q(1-NGC:N+NGC))
bx = box_interior(lb=[1-NGC], ub=[N+NGC], ngc=NGC, components=.false.)
print '(A,5I3)', 'rank 1, lo as rank 5:', bx%lo
print '(A,5I3)', 'rank 1, hi as rank 5:', bx%hi
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1-NGC], ub=[N+NGC], components=.false.)
