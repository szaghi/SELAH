call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)               ! the stepper builds the state of stage k
   u => st%stage_state()                           ! Butcher form: its own register
   r => st%rhs_register()                          ! where the residual of stage k goes
   call fill_ghosts(u)                             ! the caller: ghost cells of the stage state
   call residual(u, r, dx)                         ! the caller: r = -du/dx
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)                     ! q^{n+1}
