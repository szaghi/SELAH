!run ranks ranks
program ranks
!< Tutorial, chapter 4: the stepper on fields of rank 1 and 3, and the same field as rank 5.
use selah
implicit none
integer(I4P), parameter :: N=6, NGC=1, NV=2
real(R_P),    parameter :: DT=0.05_R_P, T_END=1._R_P
real(R_P), allocatable, target :: q(:), w(:,:,:), q5(:,:,:,:,:)
real(R_P), pointer      :: u(:), r(:), u3(:,:,:), r3(:,:,:), u5(:,:,:,:,:), r5(:,:,:,:,:)
type(stepper)           :: st, st3, st5
type(box)               :: bx, bx3, bx5
integer(I4P)            :: i, k, step

!region rank1
! a scalar 1-D field q(i), one ghost cell each side: the box of its interior, viewed as rank 5 (1, i, 1, 1, 1)
allocate(q(1-NGC:N+NGC))
bx = box_interior(lb=[1-NGC], ub=[N+NGC], ngc=NGC, components=.false.)
print '(A,5I3)', 'rank 1, lo as rank 5:', bx%lo
print '(A,5I3)', 'rank 1, hi as rank 5:', bx%hi
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1-NGC], ub=[N+NGC], components=.false.)
!endregion rank1
!region rank3
! two components on a 2-D grid, w(v, i, j): viewed as rank 5 (v, i, j, 1, 1)
allocate(w(1:NV, 1-NGC:N+NGC, 1-NGC:N+NGC))
bx3 = box_interior(lb=[1, 1-NGC, 1-NGC], ub=[NV, N+NGC, N+NGC], ngc=NGC, components=.true.)
print '(A,5I3)', 'rank 3, hi as rank 5:', bx3%hi
call st3%initialize(SCHEME_RK4, FORM_BUTCHER, lb=lbound(w), ub=ubound(w), components=.true.)
!endregion rank3
q = [(1._R_P + 0.1_R_P * i, i=1-NGC, N+NGC)]
w = 1._R_P
!region step
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
!endregion step
print '(A,F12.9)', 'rank 1, q(3) at t = 1:        ', q(3)
print '(A,F12.9)', 'exact 1/(1/q0 + t):           ', 1._R_P / (1._R_P / 1.3_R_P + T_END)
do step=1, nint(T_END/DT)
   call st3%begin_step(bx3, w)
   do k=1, st3%stages()
      call st3%begin_stage(k, bx3, DT, w)
      call st3%stage_state_view(u3) ; call st3%rhs_register_view(r3)
      r3(:,1:N,1:N) = -u3(:,1:N,1:N)                   ! w' = -w
      call st3%end_stage(k, bx3, DT, w)
   enddo
   call st3%finish_step(bx3, DT, w)
enddo
print '(A,F12.9,A,F12.9)', 'rank 3, w(2,1,1) at t = 1:    ', w(2,1,1), '   exp(-1): ', exp(-1._R_P)
!region rank5
! the same 1-D field declared rank 5 (1, i, 1, 1, 1): the same values, bit for bit
allocate(q5(1:1, 1-NGC:N+NGC, 1:1, 1:1, 1:1))
q5(1,:,1,1,1) = [(1._R_P + 0.1_R_P * i, i=1-NGC, N+NGC)]
bx5 = box_range(lb=[1, 1-NGC, 1, 1, 1], lo=[1, 1, 1, 1, 1], hi=[1, N, 1, 1, 1], components=.true.)
call st5%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=lbound(q5), ub=ubound(q5))
do step=1, nint(T_END/DT)
   call st5%begin_step(bx5, q5)
   do k=1, st5%stages()
      call st5%begin_stage(k, bx5, DT, q5)
      u5 => q5 ; r5 => st5%rhs_register()
      r5(1,1:N,1,1,1) = -u5(1,1:N,1,1,1)**2
      call st5%end_stage(k, bx5, DT, q5)
   enddo
   call st5%finish_step(bx5, DT, q5)
enddo
print '(A,L1)', 'rank 1 bitwise rank 5: ', all(q == q5(1,:,1,1,1))
!endregion rank5
endprogram ranks
