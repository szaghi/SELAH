! the same 1-D field declared rank 5 (1, i, 1, 1, 1): the same values, bit for bit
allocate(q5(1:1, 1-NGC:N+NGC, 1:1, 1:1, 1:1))
q5(1,:,1,1,1) = [(1._R_P + 0.1_R_P * i, i=1-NGC, N+NGC)]
bx5 = box_range(lb=[1, 1-NGC, 1, 1, 1], lo=[1, 1, 1, 1, 1], hi=[1, N, 1, 1, 1], components=.true.)
call st5%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=lbound(q5), ub=ubound(q5))
do step=1, nint(T_END/DT)
   call st5%begin_step(bx5, q5)
   do k=1, st5%stages()
      call st5%begin_stage(k, bx5, DT, q5)
      u5 => q5 ; r5 => st5%rhs_register()
      r5(1,1:N,1,1,1) = -u5(1,1:N,1,1,1)**2
      call st5%end_stage(k, bx5, DT, q5)
   enddo
   call st5%finish_step(bx5, DT, q5)
enddo
print '(A,L1)', 'rank 1 bitwise rank 5: ', all(q == q5(1,:,1,1,1))
