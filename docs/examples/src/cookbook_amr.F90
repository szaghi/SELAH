!run amr-reflux cookbook_amr reflux
!run amr-ghosts cookbook_amr ghosts
program cookbook_amr
!< Cookbook, AMR: refluxing with stage_weight; where the ghost cells of a stage come from.
use selah
use gideon, only : gideon_fd_d1_c2, gideon_fv_face_c2, gideon_fv_face_ul1   ! the space operators: GIDEON's
implicit none
integer(I4P), parameter :: NI=16, NB=2
real(R_P),    parameter :: PI=acos(-1._R_P)
character(16)           :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('reflux') ; call reflux
case('ghosts') ; call ghosts
case default   ; error stop 'usage: cookbook_amr reflux|ghosts'
endselect

contains
   subroutine reflux
   !< Two blocks, periodic. At their shared face the two sides compute different fluxes (as a coarse and a fine level
   !< would): the mass drifts. A flux register accumulates the difference with stage_weight; one correction of the
   !< coarse side's cell after the step restores conservation to round-off.
   real(R_P), allocatable, target :: q(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                      :: dx, dt, m0, reg, drift_before, drift_after
   integer(I4P)                   :: step, k, i

   dx = 1._R_P / (NI * NB) ; dt = 0.4_R_P * dx
   allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
   do i=1, NI
      q(1,i,1,1,1) = 1._R_P + 0.5_R_P * sin(2._R_P*PI*(i - 0.5_R_P)*dx)
      q(1,i,1,1,2) = 1._R_P + 0.5_R_P * sin(2._R_P*PI*(NI + i - 0.5_R_P)*dx)
   enddo
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
   call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
   m0 = mass(q) ; drift_before = 0._R_P ; drift_after = 0._R_P
   do step=1, 8
      !region reflux
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
      !endregion reflux
   enddo
   print '(A,ES9.2)', 'mass drift without the correction: ', drift_before
   print '(A,ES9.2)', 'mass drift with it:                ', drift_after
   call st%destroy
   endsubroutine reflux

   subroutine ghosts
   !< The ghost cells of a stage must come from the neighbours' *stage state*, not from their q^n: the order drops.
   real(R_P) :: err(2,2)
   integer(I4P) :: k
   do k=1, 2
      err(k,1) = advect(8 * k, .true.)
      err(k,2) = advect(8 * k, .false.)
   enddo
   print '(A,F5.2)', 'ghosts from the stage state: observed order ', log(err(1,1) / err(2,1)) / log(2._R_P)
   print '(A,F5.2)', 'ghosts from q^n (wrong):     observed order ', log(err(1,2) / err(2,2)) / log(2._R_P)
   endsubroutine ghosts

   function advect(nsteps, right) result(err)
   integer(I4P), intent(in)       :: nsteps
   logical,      intent(in)       :: right
   real(R_P)                      :: err
   real(R_P), allocatable, target :: q(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                      :: dx, dt, c
   integer(I4P)                   :: step, k, i, b

   dx = 1._R_P / (NI * NB) ; dt = 0.25_R_P / nsteps
   allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
   do b=1, NB
      do i=1, NI
         q(1,i,1,1,b) = sin(2._R_P*PI*((b-1)*NI + i - 0.5_R_P)*dx)
      enddo
   enddo
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
   call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
   do step=1, nsteps
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q)
         u => st%stage_state() ; r => st%rhs_register()
         !region ghosts
         if (right) then
            call fill_ghosts(u)                            ! from the neighbours' stage state: right
         else
            call fill_ghosts_from(u, q)                    ! from the neighbours' q^n: a stale exchange
         endif
         !endregion ghosts
         call residual_central(u, r, dx)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
   enddo
   c = sin(2._R_P*PI*dx) / (2._R_P*PI*dx)
   err = 0._R_P
   do b=1, NB
      do i=1, NI
         err = max(err, abs(q(1,i,1,1,b) - sin(2._R_P*PI*(((b-1)*NI + i - 0.5_R_P)*dx - c*0.25_R_P))))
      enddo
   enddo
   call st%destroy
   endfunction advect

   subroutine fill_ghosts(u)
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   call fill_ghosts_from(u, u)
   endsubroutine fill_ghosts

   subroutine fill_ghosts_from(u, src)
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: src(1:,0:,0:,0:,1:)
   integer(I4P)             :: b
   do b=1, NB
      u(1,0,1,1,b)    = src(1,NI,1,1,1 + modulo(b-2, NB))
      u(1,NI+1,1,1,b) = src(1,1,1,1,1 + modulo(b, NB))
   enddo
   endsubroutine fill_ghosts_from

   subroutine residual_central(u, r, dx)
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   integer(I4P)             :: b, i
   real(R_P)                :: d
   do b=1, NB
      do i=1, NI
         call gideon_fd_d1_c2([u(1,i-1,1,1,b), u(1,i,1,1,b), u(1,i+1,1,1,b)], dx, d) ; r(1,i,1,1,b) = -d
      enddo
   enddo
   endsubroutine residual_central

   subroutine residual(u, r, dx)
   !< Face fluxes f = (u_L + u_R)/2, except at block 2's left face, where block 2 uses the upwind f = u_L.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   real(R_P)                :: fl, fr
   integer(I4P)             :: b, i
   do b=1, NB
      do i=1, NI
         call gideon_fv_face_c2([u(1,i-1,1,1,b), u(1,i,1,1,b)], fl)
         call gideon_fv_face_c2([u(1,i,1,1,b), u(1,i+1,1,1,b)], fr)
         if (b == 2 .and. i == 1) fl = face_flux_2(u)
         r(1,i,1,1,b) = -(fr - fl) / dx
      enddo
   enddo
   endsubroutine residual

   pure function face_flux_1(u) result(f)
   !< The shared face as block 1 sees it: central (GIDEON's gideon_fv_face_c2).
   real(R_P), intent(in) :: u(1:,0:,0:,0:,1:)
   real(R_P)             :: f
   call gideon_fv_face_c2([u(1,NI,1,1,1), u(1,NI+1,1,1,1)], f)
   endfunction face_flux_1

   pure function face_flux_2(u) result(f)
   !< The shared face as block 2 sees it: upwind, its left ghost (GIDEON's gideon_fv_face_ul1).
   real(R_P), intent(in) :: u(1:,0:,0:,0:,1:)
   real(R_P)             :: f
   call gideon_fv_face_ul1([u(1,0,1,1,2)], f)
   endfunction face_flux_2

   function mass(q) result(m)
   real(R_P), intent(in) :: q(1:,0:,0:,0:,1:)
   real(R_P)             :: m
   m = sum(q(1,1:NI,1,1,:)) / (NI * NB)
   endfunction mass
endprogram cookbook_amr
