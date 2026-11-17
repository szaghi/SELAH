call st%begin_step(bx, q, mask=mask)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q, mask=mask)
   r => st%rhs_register()
   call fill_ghosts(q)
   call residual(q, r, dx)                                     ! computed everywhere: the stepper masks the update
   call st%end_stage(k, bx, dt, q, mask=mask)
enddo
call st%finish_step(bx, dt, q, mask=mask)
