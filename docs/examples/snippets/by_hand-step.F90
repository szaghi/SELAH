do s=1, bt%stages
   a_row = bt%a(s,:)                                           ! row s of A, contiguous
   call selah_stage_state(s-1, dt, a_row, y, regs, u)         ! u = y + dt sum_{j<s} a(s,j) F_j
   regs(:,s) = -LAMBDA * u                                     ! the caller's residual F_s = f(u)
enddo
call selah_step_update(bt%stages, dt, bt%b, regs, y)           ! y = y + dt sum_s b(s) F_s
