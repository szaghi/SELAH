call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   if (st%stage_state_is_q(k)) then                    ! stage 1: its state is q, no copy was made
      u => q
   else
      u => st%stage_state()
   endif
   r => st%rhs_register()
   call fill_ghosts(u)
   call residual(u, r, dx)
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q, dq=dq)                  ! the update writes dq = sum_k b_k F_k in the same pass
