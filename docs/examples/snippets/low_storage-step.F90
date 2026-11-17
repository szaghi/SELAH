call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   if (st%stage_state_is_q(k)) then                ! the low-storage forms keep the stage state in q itself
      u => q
   else                                            ! the Butcher form keeps it in a register
      u => st%stage_state()
   endif
   r => st%rhs_register()
   call fill_ghosts(u)
   call residual(u, r, dx)
   call st%end_stage(k, bx, dt, q)                 ! the low-storage forms advance q here
enddo
call st%finish_step(bx, dt, q)                     ! the Butcher form advances q here
