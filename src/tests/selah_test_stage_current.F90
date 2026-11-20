!< SELAH test: the look-ahead state and the charge-conserving stage currents, on the host or the device.

#include "fundal.H"
program selah_test_stage_current
!< SELAH test: the look-ahead state and the charge-conserving stage currents, on the host or the device.
!<
!< A particle-in-cell code keeps Gauss's law at every stage when the weighted sum of its stage currents is the
!< charge-conserving current of the displacement q^n -> (the next stage's state) (notes/companions-pushers-pml.md §4).
!< A stencil-free model of it: a particle x' = v, v' = -x (a companion), a "field" D' = -J (one cell), and the "deposit"
!< H(x^n -> X) = (sin X - sin x^n)/dt, a nonlinear stand-in for Esirkepov's current: the conserved charge is
!< D + sin x.
!<
!<+ **Gauss's law at every stage and step.** For every Runge–Kutta scheme (Butcher form), 10 steps of dt = 0.3: at
!<  stage k the particle's residual, its `lookahead(k)` X, H = (sin X - sin x^n)/dt, the field's `stage_current(k)`
!<  J_k, the field's residual -J_k. D + sin x is constant to round-off in every stage state and after every step.
!<+ **The look-ahead is the next stage's state.** lookahead(k) equals the state of stage k+1 bit for bit, and at the
!<  last stage the step's result to round-off.
!<+ **Admissibility.** Every catalogue scheme is admissible; a tableau with a(3,2) = 0 is not; a splitting is not.
use penf, only : I4P
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
implicit none
integer(I4P), parameter :: L5(5) = 1_I4P                 !< Bounds of the one-cell field (rank 5).
integer(I4P)            :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
print '(A)', 'Gauss''s law (D + sin x) at every stage and step:'
do s=1, size(RK_SCHEMES)
   call gauss(trim(RK_SCHEMES(s)))
enddo
call admissibility
if (failures > 0) error stop 'selah_test_stage_current: failed'
print '(A)', 'selah_test_stage_current: ok'

contains
   subroutine gauss(name)
   !< 10 steps of the particle-field model: the invariant, and the look-ahead against the next stage.
   character(*), intent(in) :: name
   real(R_P),    parameter  :: DT = 0.3_R_P
   type(stepper)            :: st, sp
   type(box)                :: bd, bp
   real(R_P), pointer       :: d(:,:,:,:,:)=>null(), p(:,:)=>null(), pt(:,:)=>null(), h(:,:,:,:,:)=>null()
   real(R_P), pointer       :: jall(:,:,:,:,:,:)=>null(), jk(:,:,:,:,:)=>null()
   real(R_P), pointer       :: ud(:,:,:,:,:), up(:,:), rd(:,:,:,:,:), rp(:,:)
   real(R_P)                :: hd(1,1,1,1,1), hp(2,1), hpt(2,1), c0, worst, xn, prev(2,1)
   integer(I4P)             :: step, k, ns, ierr
   logical                  :: next_same, last_close

   call dev_alloc(fptr_dev=d, ubounds=L5, lbounds=L5, ierr=ierr)
   call dev_alloc(fptr_dev=h, ubounds=L5, lbounds=L5, ierr=ierr)
   call dev_alloc(fptr_dev=jk, ubounds=L5, lbounds=L5, ierr=ierr)
   call dev_alloc(fptr_dev=p, ubounds=[2_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   call dev_alloc(fptr_dev=pt, ubounds=[2_I4P, 1_I4P], lbounds=[1_I4P, 1_I4P], ierr=ierr)
   bd = box_interior(LAYOUT_V_FIRST, 1_I4P, 1_I4P, 1_I4P, 1_I4P, 0_I4P, 1_I4P)
   bp = box_interior(lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], ngc=0_I4P, components=.true.)
   call st%initialize(name, FORM_BUTCHER, lb=L5, ub=L5)
   call sp%initialize_companion(st, lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], components=.true.)
   ns = st%stages()
   call dev_alloc(fptr_dev=jall, ubounds=[L5, ns], lbounds=[L5, 1_I4P], ierr=ierr)
   hd = 0.2_R_P ; hp(:,1) = [0.7_R_P, 0.4_R_P]
   call dev_memcpy_to_device(dst=d, src=hd) ; call dev_memcpy_to_device(dst=p, src=hp)
   c0 = hd(1,1,1,1,1) + sin(hp(1,1))
   worst = 0._R_P ; next_same = .true. ; last_close = .true.
   do step=1, 10
      xn = hp(1,1)
      call st%begin_step(bd, d) ; call sp%begin_step(bp, p)
      do k=1, ns
         call st%begin_stage(k, bd, DT, d) ; call sp%begin_stage(k, bp, DT, p)
         if (st%stage_state_is_q(k)) then
            ud => d
         else
            ud => st%stage_state()
         endif
         if (sp%stage_state_is_q(k)) then
            up => p
         else
            call sp%stage_state_view(up)
         endif
         call dev_memcpy_from_device(dst=hd, src=ud) ; call dev_memcpy_from_device(dst=hp, src=up)
         worst = max(worst, abs(hd(1,1,1,1,1) + sin(hp(1,1)) - c0))
         if (k > 1) next_same = next_same .and. all(hp == prev)
         ! the particle's residual, then its look-ahead, the deposit H, the stage current, the field's residual
         call sp%rhs_register_view(rp)
         call oscillator(up, rp)
         call sp%lookahead(k, bp, DT, p, pt)
         call dev_memcpy_from_device(dst=hpt, src=pt)
         prev = hpt
         hd = (sin(hpt(1,1)) - sin(xn)) / DT
         call dev_memcpy_to_device(dst=h, src=hd)
         call st%stage_current(k, bd, h, jall, jk)
         call selah_copy(bd, jk, jall(:,:,:,:,:,k))
         rd => st%rhs_register()
         call minus(jk, rd)
         call st%end_stage(k, bd, DT, d) ; call sp%end_stage(k, bp, DT, p)
      enddo
      call st%finish_step(bd, DT, d) ; call sp%finish_step(bp, DT, p)
      call dev_memcpy_from_device(dst=hd, src=d) ; call dev_memcpy_from_device(dst=hp, src=p)
      worst = max(worst, abs(hd(1,1,1,1,1) + sin(hp(1,1)) - c0))
      last_close = last_close .and. all(abs(hp - prev) <= 8*epsilon(1._R_P))
   enddo
   print '(A,ES10.2)', '        '//name//': max |D + sin x - (D + sin x)^0| over every stage and step', worst
   call check(worst <= 64*epsilon(1._R_P), name//': Gauss''s law at every stage and step')
   call check(next_same, name//': lookahead(k) is the state of stage k+1, bitwise')
   call check(last_close, name//': lookahead(s) is the step''s result to round-off')
   call st%destroy ; call sp%destroy
   call dev_free(d) ; call dev_free(h) ; call dev_free(jk) ; call dev_free(p) ; call dev_free(pt) ; call dev_free(jall)
   endsubroutine gauss

   subroutine oscillator(u, r)
   !< r = (v, -x).
   real(R_P), intent(in)    :: u(1:,1:)
   real(R_P), intent(inout) :: r(1:,1:)

   !$acc serial DEVICEVAR(u, r)
   r(1,1) = u(2,1)
   r(2,1) = -u(1,1)
   !$acc end serial
   endsubroutine oscillator

   subroutine minus(j, r)
   !< r = -j (one cell).
   real(R_P), intent(in)    :: j(1:,1:,1:,1:,1:)
   real(R_P), intent(inout) :: r(1:,1:,1:,1:,1:)

   !$acc serial DEVICEVAR(j, r)
   r(1,1,1,1,1) = -j(1,1,1,1,1)
   !$acc end serial
   endsubroutine minus

   subroutine admissibility
   !< Every catalogue scheme admissible; a(3,2) = 0 not; a splitting not.
   type(stepper) :: st
   integer(I4P)  :: k
   logical       :: all_ok

   all_ok = .true.
   do k=1, size(RK_SCHEMES)
      call st%initialize(trim(RK_SCHEMES(k)), FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
      all_ok = all_ok .and. st%stage_current_admissible()
      call st%destroy
   enddo
   call check(all_ok, 'every catalogue scheme admits stage currents')
   call st%initialize(SCHEME_SSP_33, FORM_BUTCHER, lb=[1_I4P], ub=[1_I4P], components=.false.)
   st%bt%a(3,2) = 0._R_P
   call check(.not. st%stage_current_admissible(), 'a(3,2) = 0: no stage currents')
   call st%destroy
   call st%initialize(SCHEME_STRANG, FORM_SPLITTING, lb=[1_I4P, 1_I4P], ub=[2_I4P, 1_I4P], components=.true., &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   call check(.not. st%stage_current_admissible(), 'a splitting: no stage currents')
   call st%destroy
   endsubroutine admissibility
endprogram selah_test_stage_current
