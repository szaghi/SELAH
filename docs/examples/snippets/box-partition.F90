part = box_components(bx, 2_I4P, 2_I4P)                           ! component 2 only: a partition
y = SENTINEL
call selah_copy(part, x, y)
print '(A,I0,A,I0)', 'component range 2..2: cells updated ', count(y /= SENTINEL), ', component 1 untouched: ', &
      merge(1, 0, all(y(1,:,:,:,:) == SENTINEL))
