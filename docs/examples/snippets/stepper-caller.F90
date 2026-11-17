subroutine fill_ghosts(u)
!< Periodic ghost cells along i: each block's ghosts from its neighbours.
real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
integer(I4P)             :: b
do b=1, NB
   u(1,0,1,1,b)    = u(1,NI,1,1,1 + modulo(b-2, NB))
   u(1,NI+1,1,1,b) = u(1,1,1,1,1 + modulo(b, NB))
enddo
endsubroutine fill_ghosts

subroutine residual(u, r, dx)
!< Central flux form: r_i = -((u_i + u_{i+1})/2 - (u_{i-1} + u_i)/2)/dx.
real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
real(R_P), intent(in)    :: dx
integer(I4P)             :: b, i
do b=1, NB
   do i=1, NI
      r(1,i,1,1,b) = -(u(1,i+1,1,1,b) - u(1,i-1,1,1,b)) / (2._R_P*dx)
   enddo
enddo
endsubroutine residual
