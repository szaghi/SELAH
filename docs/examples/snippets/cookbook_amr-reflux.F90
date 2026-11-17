reg = 0._R_P                                         ! the flux register of the shared face, zeroed per step
call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   u => st%stage_state() ; r => st%rhs_register()
   call fill_ghosts(u)
   call residual(u, r, dx)                           ! block 1: central flux; block 2: upwind flux at the face
   reg = reg + st%stage_weight(k) * (face_flux_2(u) - face_flux_1(u))   ! fine minus coarse, weighted
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)
drift_before = max(drift_before, abs(mass(q) - m0))
q(1,NI,1,1,1) = q(1,NI,1,1,1) - dt / dx * reg       ! reflux: block 1's last cell takes block 2's flux
drift_after = max(drift_after, abs(mass(q) - m0))
