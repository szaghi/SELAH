do step=1, nint(T_END/DT)
   call st%begin_step(bx, q)
   do k=1, st%stages()
      call st%begin_stage(k, bx, DT, q)
      if (st%stage_state_is_q(k)) then ; u => q ; else ; call st%stage_state_view(u) ; endif
      call st%rhs_register_view(r)                     ! rank-1 pointers, with the field's bounds
      r(1:N) = -u(1:N)**2                              ! yours: y' = -y^2, cell by cell
      call st%end_stage(k, bx, DT, q)
   enddo
   call st%finish_step(bx, DT, q)
enddo
