!run flux_register flux_register
program flux_register
!< Tutorial, chapter 8: the boundary fluxes of each stage, weighted by stage_weight, telescope to the step's change of
!< the mass. This is the identity an AMR flux register (reflux) relies on.
use selah
use gideon, only : gideon_fv_face_c2   ! the space operator: GIDEON's centred face value
implicit none
integer(I4P), parameter        :: NI=16, NB=4, NSTEPS=4
real(R_P), allocatable, target :: q(:,:,:,:,:)
type(stepper)                  :: st
type(box)                      :: bx
real(R_P), pointer             :: u(:,:,:,:,:), r(:,:,:,:,:)
real(R_P)                      :: dt, dx, m0, m1, flux_w, flux_u, defect_w, defect_u
integer(I4P)                   :: step, k, b, i

dx = 1._R_P / (NI * NB) ; dt = 0.5_R_P * dx
allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
do b=1, NB
   do i=1, NI
      q(1,i,1,1,b) = 1._R_P + 0.5_R_P * sin(6.283185307179586_R_P * ((b-1)*NI + i - 0.5_R_P) * dx)
   enddo
enddo
bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, 1_I4P, NB)
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
defect_w = 0._R_P ; defect_u = 0._R_P
do step=1, NSTEPS
   m0 = sum(q(1,1:NI,1,1,:)) * dx
   flux_w = 0._R_P ; flux_u = 0._R_P
   call st%begin_step(bx, q)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, q)
      u => st%stage_state() ; r => st%rhs_register()
      call fill_ghosts(u)                                          ! inflow 1 on the left, outflow on the right
      call residual(u, r, dx)
      !region register
      ! the flux register accumulates F_right - F_left of stage k, with the stage's effective weight
      flux_w = flux_w + st%stage_weight(k) * boundary_flux(u)
      flux_u = flux_u + boundary_flux(u) / st%stages()             ! the wrong weights: uniform, for comparison
      !endregion register
      call st%end_stage(k, bx, dt, q)
   enddo
   call st%finish_step(bx, dt, q)
   m1 = sum(q(1,1:NI,1,1,:)) * dx
   defect_w = max(defect_w, abs((m1 - m0) + dt * flux_w))
   defect_u = max(defect_u, abs((m1 - m0) + dt * flux_u))
enddo
print '(A,ES9.2)', 'mass change + dt sum_k stage_weight(k) flux_k: ', defect_w
print '(A,ES9.2)', 'the same with uniform weights 1/s:            ', defect_u
call st%destroy

contains
   subroutine fill_ghosts(u)
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   integer(I4P)             :: b
   do b=1, NB
      u(1,0,1,1,b)    = merge(1._R_P, u(1,NI,1,1,max(b-1,1)), b == 1)
      u(1,NI+1,1,1,b) = merge(u(1,NI,1,1,b), u(1,1,1,1,min(b+1,NB)), b == NB)
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r, dx)
   !< Flux form with central face fluxes f_{i+1/2} = (u_i + u_{i+1})/2.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   real(R_P), intent(in)    :: dx
   integer(I4P)             :: b, i
   real(R_P)                :: fl, fr
   do b=1, NB
      do i=1, NI
         call gideon_fv_face_c2([u(1,i-1,1,1,b), u(1,i,1,1,b)], fl)
         call gideon_fv_face_c2([u(1,i,1,1,b), u(1,i+1,1,1,b)], fr)
         r(1,i,1,1,b) = -(fr - fl) / dx
      enddo
   enddo
   endsubroutine residual

   function boundary_flux(u) result(f)
   !< F_right - F_left: the faces of the domain.
   real(R_P), intent(in) :: u(1:,0:,0:,0:,1:)
   real(R_P)             :: f
   real(R_P)             :: fl, fr
   call gideon_fv_face_c2([u(1,NI,1,1,NB), u(1,NI+1,1,1,NB)], fr)
   call gideon_fv_face_c2([u(1,0,1,1,1), u(1,1,1,1,1)], fl)
   f = fr - fl
   endfunction boundary_flux
endprogram flux_register
