allocate(y(1:2, 0:2, 0:2, 0:2, 1:1))                          ! one cell, two components: q (1) and p (2)
y = 0._R_P ; y(2,1,1,1,1) = 1._R_P
bx = box_interior(LAYOUT_V_FIRST, 2_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
! partition 1 = components 1..1 (q), partition 2 = components 2..2 (p): column p holds (v1, v2)
call st%initialize(name, FORM_SPLITTING, lb=lbound(y), ub=ubound(y), partitions=reshape([1, 1, 2, 2], [2, 2]))
