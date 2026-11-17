! V_LAST layout (b,i,j,k,v): the blocks innermost (ADAM's device layout); the same box, built for that layout
deallocate(x, y)
allocate(x(1:NB, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NV))
allocate(y, mold=x)
x = 1._R_P ; y = SENTINEL
bx = box_interior(LAYOUT_V_LAST, NV, NI, NJ, NK, NGC, 3_I4P)
call selah_axpy(bx, 2._R_P, x, y)
print '(A,I0,A,I0)', 'V_LAST: component dimension ', component_dim(bx), ', cells updated ', count(y /= SENTINEL)
