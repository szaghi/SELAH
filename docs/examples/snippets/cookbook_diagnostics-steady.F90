do step=1, 10000
   call st%begin_step(bx, u)
   do k=1, st%stages()
      call st%begin_stage(k, bx, 0.5_R_P, u)
      s_ => st%stage_state() ; r => st%rhs_register()
      r(1,1:N,1,1,1) = -(s_(1,1:N,1,1,1) - g(1,1:N,1,1,1))
      call st%end_stage(k, bx, 0.5_R_P, u)
   enddo
   call st%finish_step(bx, 0.5_R_P, u, dq=dq)
   if (maxval(abs(dq(1,1:N,1,1,1))) < TOL) exit             ! converged: the step changes nothing any more
enddo
