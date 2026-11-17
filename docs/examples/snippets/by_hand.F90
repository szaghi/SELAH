program by_hand
!< Tutorial, chapter 1: one step of SSP(5,4) by hand, with the flat kernels, on y' = -lambda y.
use selah
implicit none
integer(I4P), parameter :: N=3                                      ! three independent equations
real(R_P),    parameter :: LAMBDA(N) = [0.5_R_P, 1._R_P, 2._R_P]    ! their rates
type(butcher_tableau)   :: bt                                       ! the scheme's coefficients
real(R_P), allocatable  :: regs(:,:)                                ! one residual register per stage
real(R_P)               :: y(N), u(N), err(4), dt, a_row(5)
integer(I4P)            :: nsteps, step, s, k

bt = butcher(SCHEME_SSP_54)
allocate(regs(N, bt%stages))
print '(A,I0,A,I0)', 'ssp-54: stages ', bt%stages, ', order ', bt%order
do k=1, 4
   nsteps = 5 * 2**(k-1)
   dt = 1._R_P / nsteps
   y = 1._R_P
   do step=1, nsteps
      do s=1, bt%stages
         a_row = bt%a(s,:)                                           ! row s of A, contiguous
         call selah_stage_state(s-1, dt, a_row, y, regs, u)         ! u = y + dt sum_{j<s} a(s,j) F_j
         regs(:,s) = -LAMBDA * u                                     ! the caller's residual F_s = f(u)
      enddo
      call selah_step_update(bt%stages, dt, bt%b, regs, y)           ! y = y + dt sum_s b(s) F_s
   enddo
   err(k) = maxval(abs(y - exp(-LAMBDA)))
   print '(A,I3,A,ES10.3)', 'steps ', nsteps, ': error at t = 1 ', err(k)
enddo
print '(A,3F6.2)', 'observed orders:', log(err(1:3) / err(2:4)) / log(2._R_P)
endprogram by_hand
