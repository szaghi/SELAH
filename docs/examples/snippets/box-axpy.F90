y = SENTINEL
call selah_axpy(bx, 2._R_P, x, y)                                 ! y = y + 2 x, in the box only
print '(A,I0,A,I0)', 'cells updated: ', count(y /= SENTINEL), ' of ', size(y)
print '(A,I0)', 'expected (2 components x 4x3x2 cells x 3 blocks): ', NV*NI*NJ*NK*3
print '(A,L1)', 'ghosts and block 4 untouched: ', all(y(:,:,:,:,4) == SENTINEL) .and. all(y(:,0,:,:,:) == SENTINEL)
