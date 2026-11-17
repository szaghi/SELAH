subroutine fill_ghosts(u)
!< Periodic ghost cells, one device thread per block.
real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
integer(I4P)             :: b

!$acc parallel loop independent DEVICEVAR(u)
!$omp OMPLOOP DEVICEPTR(u)
do b=1, NB
   u(1,0,1,1,b)    = u(1,NI,1,1,1 + modulo(b-2, NB))
   u(1,NI+1,1,1,b) = u(1,1,1,1,1 + modulo(b, NB))
enddo
endsubroutine fill_ghosts

subroutine residual(u, r, dx)
!< The central residual, one device thread per cell.
real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
real(R_P), intent(in)    :: dx
integer(I4P)             :: b, i

!$acc parallel loop independent collapse(2) DEVICEVAR(u, r)
!$omp OMPLOOP collapse(2) DEVICEPTR(u, r)
do b=1, NB
   do i=1, NI
      r(1,i,1,1,b) = -(u(1,i+1,1,1,b) - u(1,i-1,1,1,b)) / (2._R_P*dx)
   enddo
enddo
endsubroutine residual
