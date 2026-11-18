! two components on a 2-D grid, w(v, i, j): viewed as rank 5 (v, i, j, 1, 1)
allocate(w(1:NV, 1-NGC:N+NGC, 1-NGC:N+NGC))
bx3 = box_interior(lb=[1, 1-NGC, 1-NGC], ub=[NV, N+NGC, N+NGC], ngc=NGC, components=.true.)
print '(A,5I3)', 'rank 3, hi as rank 5:', bx3%hi
call st3%initialize(SCHEME_RK4, FORM_BUTCHER, lb=lbound(w), ub=ubound(w), components=.true.)
