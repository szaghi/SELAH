allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
q = 0._R_P
do b=1, NB
   do i=1, NI
      q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
   enddo
enddo
call dev_alloc(fptr_dev=q_d, ubounds=ubound(q), lbounds=lbound(q), ierr=ierr)   ! the field on the device
call dev_memcpy_to_device(dst=q_d, src=q)
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))      ! registers through FUNDAL: device
