program splittings
!< Tutorial, chapter 6: a splitting scheme through the stepper, on the oscillator q' = p, p' = -q.
use selah
implicit none
character(15)           :: names(3) = [character(15) :: SCHEME_STRANG, SCHEME_YOSHIDA_4, SCHEME_BLANES_MOAN_S10]
real(R_P)               :: err(3)
integer(I4P)            :: s, k, n

do s=1, 3
   n = merge(4, 16, s == 3)
   do k=1, 3
      err(k) = oscillate(trim(names(s)), n * 2**(k-1))
   enddo
   print '(A15,A,3ES10.2,A,2F5.2)', names(s), ' errors', err, '  orders', log(err(1:2) / err(2:3)) / log(2._R_P)
enddo

contains
   function oscillate(name, nsteps) result(err)
   !< nsteps steps to t = 1 from (q, p) = (0, 1); the error against (sin t, cos t).
   character(*), intent(in)       :: name
   integer(I4P), intent(in)       :: nsteps
   real(R_P)                      :: err
   real(R_P), allocatable, target :: y(:,:,:,:,:)
   type(stepper)                  :: st
   type(box)                      :: bx
   real(R_P), pointer             :: r(:,:,:,:,:)
   real(R_P)                      :: dt
   integer(I4P)                   :: step, k

   allocate(y(1:2, 0:2, 0:2, 0:2, 1:1))                          ! one cell, two components: q (1) and p (2)
   y = 0._R_P ; y(2,1,1,1,1) = 1._R_P
   bx = box_interior(LAYOUT_V_FIRST, 2_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
   ! partition 1 = components 1..1 (q), partition 2 = components 2..2 (p): column p holds (v1, v2)
   call st%initialize(name, FORM_SPLITTING, lb=lbound(y), ub=ubound(y), partitions=reshape([1, 1, 2, 2], [2, 2]))
   dt = 1._R_P / nsteps
   do step=1, nsteps
      call st%begin_step(bx, y)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, y)
         r => st%rhs_register()
         select case(st%stage_partition(k))                      ! the flow of stage k advances one partition:
         case(1) ; r(1,1,1,1,1) = y(2,1,1,1,1)                   !   q' = p  (its residual depends on p only)
         case(2) ; r(2,1,1,1,1) = -y(1,1,1,1,1)                  !   p' = -q (on q only)
         endselect
         call st%end_stage(k, bx, dt, y)                         ! the partition's components += c_k dt r
      enddo
      call st%finish_step(bx, dt, y)
   enddo
   err = max(abs(y(1,1,1,1,1) - sin(1._R_P)), abs(y(2,1,1,1,1) - cos(1._R_P)))
   call st%destroy
   endfunction oscillate
endprogram splittings
