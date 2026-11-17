allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))              ! (v,i,j,k,b): 1 component, NI x 1 x 1 cells, 1 ghost
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
call st%initialize(scheme, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
