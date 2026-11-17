call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, 0.1_R_P, q)
   u => st%stage_state() ; r => st%rhs_register()
   ! the stage state is q itself or a register, never both; never the residual register
   if (st%stage_state_is_q(k) .eqv. associated(u)) error stop 'stage state ambiguous'
   if (associated(u, r)) error stop 'residual aliases the stage state'
   print '(A,I0,A,L1,A,L1)', 'stage ', k, ': state is q ', st%stage_state_is_q(k), ', stage_state() associated ', &
         associated(u)
   r = 0._R_P
   call st%end_stage(k, bx, 0.1_R_P, q)
enddo
call st%finish_step(bx, 0.1_R_P, q)
