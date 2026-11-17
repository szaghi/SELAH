bx = box_interior(LAYOUT_V_FIRST, NV, NI, NJ, NK, NGC, 3_I4P)   ! the interior of the first 3 blocks of the pool
print '(A,5I4)', 'lb (allocation lower bounds):', bx%lb
print '(A,5I4)', 'lo (loop start):             ', bx%lo
print '(A,5I4)', 'hi (loop end):               ', bx%hi
