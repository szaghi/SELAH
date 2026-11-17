!< SELAH test: the stepper, driven by a caller as a solver would drive it, on the host or the device.

#include "fundal.H"

program selah_test_stepper
!< SELAH test: the stepper, driven by a caller as a solver would drive it, on the host or the device.
!<
!< The caller here is a 1-D finite-volume advection u_t + u_x = 0 on NB blocks of NI cells with one ghost cell, its
!< ghost fill and residual written as device loops (FUNDAL arrays: device memory in an OpenACC build).
!<
!<+ **Order.** On the periodic domain the central residual of a sine has an exact semi-discrete solution,
!<  sin(2 pi (x - c t)) with c = sin(2 pi dx)/(2 pi dx): the error is the time integrator's alone. Every scheme in
!<  every form it has reaches its order; the splittings reach theirs on a two-partition oscillator (q' = p, p' = -q).
!<+ **Contract.** Ghost cells hold sentinels through a whole step wherever the caller does not fill them; the residual
!<  register is distinct from the stage state; state_aliases_q is false for the Butcher form only.
!<+ **Stages on q.** The Butcher form with stages_on_q=.true. ('butcher-q'): stage_state_is_q(k) is true exactly for
!<  the zero rows of A, stage_state() is null during those stages, and the solution is bitwise the one of the plain
!<  Butcher form (the stage-1 copy only moves the same values).
!<+ **Telescoping.** With inflow on the left and outflow on the right, the change of the mass over a step equals
!<  -dt sum_k stage_weight(k) (F_right - F_left)_k, the stage fluxes at the boundary faces: the identity a flux
!<  register relies on, for every Runge–Kutta form (round-off).
!<
!< The orders are checked in R8P builds only (R4P coefficients cannot show them).
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use selah
use selah_test_tools, only : check, failures
use selah_vv_cases
implicit none
integer                 :: u_ord, u_tel !< Report files (SELAH_VV_DIR): stepper-orders.dat, telescoping.dat.
integer(I4P), parameter :: NI=16, NB=4, NGC=1                 !< Cells per block, blocks, ghosts.
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P),    parameter :: DX = 1._R_P / real(NI*NB, R_P)     !< Cell size.
real(R_P),    parameter :: SENTINEL = -777._R_P               !< Ghost value nobody must touch.
character(9), parameter :: RK_FORMS(5) = [character(9) :: 'butcher', 'butcher-q', 'two-reg', '2n', 'ket']
integer(I4P)            :: s

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
print '(A)', 'device build'
#endif
if (vv_report()) then
   call vv_open('stepper-orders.dat', u_ord, 'scheme form steps error (periodic advection, 4 halvings)')
   call vv_open('telescoping.dat', u_tel, 'scheme form step mass_change+dt*sum_k(w_k*flux_k) uniform_weights_defect')
endif
do s=1, size(RK_SCHEMES)
   call rk_scheme(trim(RK_SCHEMES(s)))
enddo
if (vv_report()) then
   close(u_ord) ; close(u_tel)
endif
do s=1, size(SPLITTING_SCHEMES)
   call split_scheme(trim(SPLITTING_SCHEMES(s)))
enddo
if (failures > 0) error stop 'selah_test_stepper: failed'
print '(A)', 'selah_test_stepper: ok'

contains
   function form_of(tag) result(form)
   character(*), intent(in) :: tag
   character(:), allocatable :: form
   select case(tag)
   case('butcher', 'butcher-q') ; form = FORM_BUTCHER
   case('two-reg') ; form = FORM_SHU_OSHER_2R
   case('2n')      ; form = FORM_WILLIAMSON_2N
   case default    ; form = FORM_KETCHESON_2R
   endselect
   endfunction form_of

   subroutine rk_scheme(name)
   !< Every form of an RK scheme: order (periodic), contract and telescoping (inflow/outflow).
   character(*), intent(in) :: name
   real(R_P)                :: errors(4), orders(3)
   integer(I4P)             :: f, n, k
   type(butcher_tableau)    :: bt
   real(R_P), allocatable   :: q_plain(:,:,:,:,:), q_on_q(:,:,:,:,:)
   logical                  :: on_q

   bt = butcher(name)
   do f=1, size(RK_FORMS)
      if (registers_needed(name, form_of(trim(RK_FORMS(f)))) < 0) cycle
      associate(tag=>name//' '//trim(RK_FORMS(f)))
      on_q = trim(RK_FORMS(f)) == 'butcher-q'
#if !defined _R_P_IS_R4P
      n = merge(8_I4P, 16_I4P, bt%order >= 4)
      do k=1, 4
         errors(k) = advect(name, form_of(trim(RK_FORMS(f))), n * 2_I4P**(k-1), .true., on_q)
      enddo
      orders = log(errors(1:3) / errors(2:4)) / log(2._R_P)
      if (vv_report()) write(u_ord, '(A,1X,A,4(I6,ES14.6))') name, trim(RK_FORMS(f)), (n * 2**(k-1), errors(k), k=1, 4)
      print '(A,3F7.3)', '        '//tag//' observed orders:', orders
      call check(orders(3) >= real(bt%order, R_P) - 0.2_R_P, tag//': observed order through the stepper')
#endif
      errors(1) = advect(name, form_of(trim(RK_FORMS(f))), 4_I4P, .false., on_q) ! contract and telescoping, inside
      endassociate
   enddo
   ! the Butcher form with and without the stages on q: the same solution in the cells, bit for bit (the ghost cells
   ! differ by design: with the stages on q the caller fills those of q at stage 1)
   errors(1) = advect(name, FORM_BUTCHER, 8_I4P, .true., .false., q_plain)
   errors(2) = advect(name, FORM_BUTCHER, 8_I4P, .true., .true., q_on_q)
   call check(all(q_plain(:,1:NI,1,1,:) == q_on_q(:,1:NI,1,1,:)), &
              name//' butcher-q: the solution bitwise that of the plain Butcher form')
   endsubroutine rk_scheme

   function advect(name, form, nsteps, periodic, on_q, q_final) result(err)
   !< Advect to t = 0.25 in nsteps steps. Periodic: return the error against the semi-discrete solution. Otherwise
   !< (inflow 1 on the left, outflow on the right, a contract check): check the ghosts, the registers and the
   !< telescoping of the mass at every step.
   character(*), intent(in) :: name, form
   integer(I4P), intent(in) :: nsteps
   logical,      intent(in) :: periodic
   logical,      intent(in) :: on_q                    !< Butcher form: initialise with stages_on_q.
   real(R_P), allocatable, intent(out), optional :: q_final(:,:,:,:,:) !< The solution at the end.
   real(R_P)                :: err
   type(butcher_tableau)    :: bt
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q_d(:,:,:,:,:)=>null(), u_d(:,:,:,:,:)=>null(), r_d(:,:,:,:,:)=>null()
   real(R_P), allocatable   :: q(:,:,:,:,:), h(:,:,:,:,:)
   real(R_P)                :: dt, t, c, m0, m1, flux, flux_u, fb, x
   integer(I4P)             :: step, k, i, b, ierr
   logical                  :: ghosts_ok, distinct_ok, telescoping_ok, on_q_ok

   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, NI, 1_I4P, 1_I4P, NGC, NB)
   allocate(q(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
   q = SENTINEL
   do b=1, NB
      do i=1, NI
         x = (real((b-1)*NI + i, R_P) - 0.5_R_P) * DX
         q(1,i,1,1,b) = merge(sin(2._R_P*PI*x), 1._R_P + 0.5_R_P*sin(2._R_P*PI*x), periodic)
      enddo
   enddo
   call dev_alloc(fptr_dev=q_d, ubounds=[1_I4P, NI+1_I4P, 2_I4P, 2_I4P, NB], lbounds=[1_I4P, 0_I4P, 0_I4P, 0_I4P, 1_I4P], &
                  ierr=ierr)
   call dev_memcpy_to_device(dst=q_d, src=q)
   call st%initialize(name, form, lb=[1_I4P, 0_I4P, 0_I4P, 0_I4P, 1_I4P], ub=[1_I4P, NI+1_I4P, 2_I4P, 2_I4P, NB], &
                      stages_on_q=on_q)
   dt = 0.25_R_P / real(nsteps, R_P)
   ghosts_ok = .true. ; distinct_ok = .true. ; telescoping_ok = .true. ; on_q_ok = .true.
   call check(st%state_aliases_q() .eqv. (form /= FORM_BUTCHER), name//' '//form//': state_aliases_q')
   bt = butcher(name)
   do k=1, st%stages() ! q is the stage state: every stage of the aliasing forms, the zero rows of A with on_q
      if (st%stage_state_is_q(k) .neqv. (form /= FORM_BUTCHER .or. (on_q .and. all(bt%a(k,:) == 0._R_P)))) &
         on_q_ok = .false.
   enddo
   if (.not. periodic) call sentinel_registers(st)
   do step=1, nsteps
      if (.not. periodic) m0 = mass(q_d)
      flux = 0._R_P ; flux_u = 0._R_P
      call st%begin_step(bx, q_d)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q_d)
         if (st%stage_state_is_q(k)) then
            if (associated(st%stage_state())) on_q_ok = .false. ! null on a stage on q, never a stale register
            u_d => q_d
         else
            u_d => st%stage_state()
         endif
         r_d => st%rhs_register()
         if (associated(r_d, u_d)) distinct_ok = .false.
         call fill_ghosts(u_d, periodic)                   ! the caller's ghost fill (clause 3)
         call residual(u_d, r_d)                           ! the caller's residual (clause 4)
         if (.not. periodic) then
            fb = boundary_flux(u_d)
            flux = flux + st%stage_weight(k) * fb
            flux_u = flux_u + fb / st%stages()   ! uniform weights: the control of the report
         endif
         call st%end_stage(k, bx, dt, q_d)
      enddo
      call st%finish_step(bx, dt, q_d)
      if (.not. periodic) then
         m1 = mass(q_d)
         if (abs((m1 - m0) + dt * flux) > 256._R_P * epsilon(1._R_P) * max(abs(m1), 1._R_P)) telescoping_ok = .false.
         if (vv_report() .and. .not. on_q) write(u_tel, '(A,1X,A,I4,2ES14.5)') name, form, step, &
                                               abs((m1 - m0) + dt * flux), abs((m1 - m0) + dt * flux_u)
         if (.not. registers_untouched(st)) ghosts_ok = .false.
      endif
   enddo
   err = 0._R_P
   if (periodic) then
      allocate(h, mold=q) ; call dev_memcpy_from_device(dst=h, src=q_d)
      t = 0.25_R_P
      c = sin(2._R_P*PI*DX) / (2._R_P*PI*DX)
      do b=1, NB
         do i=1, NI
            x = (real((b-1)*NI + i, R_P) - 0.5_R_P) * DX
            err = max(err, abs(h(1,i,1,1,b) - sin(2._R_P*PI*(x - c*t))))
         enddo
      enddo
   else
      allocate(h, mold=q) ; call dev_memcpy_from_device(dst=h, src=q_d)
      ! the j/k ghost layers of q are never filled by the caller: still the sentinel (clause 2)
      if (any(h(:,:,0,:,:) /= SENTINEL) .or. any(h(:,:,2,:,:) /= SENTINEL)) ghosts_ok = .false.
      call check(ghosts_ok, name//' '//form//': no ghost cell written by the stepper (q and registers)')
      call check(distinct_ok, name//' '//form//': residual register distinct from the stage state')
      call check(telescoping_ok, name//' '//form//': mass change = -dt sum_k w_k (boundary fluxes)_k')
      call check(on_q_ok, name//' '//form//merge(' on q', '     ', on_q)//': stage_state_is_q, stage_state() null on q')
   endif
   if (present(q_final)) then
      allocate(q_final, mold=q) ; call dev_memcpy_from_device(dst=q_final, src=q_d)
   endif
   call st%destroy
   call dev_free(q_d)
   endfunction advect

   subroutine fill_ghosts(u, periodic)
   !< The caller's ghost fill along i: periodic, or inflow 1 (left) and outflow (right); j and k ghosts untouched.
   real(R_P), intent(inout) :: u(1:,0:,0:,0:,1:)
   logical,   intent(in)    :: periodic
   integer(I4P)             :: b, bl, br

   !$acc parallel loop independent DEVICEVAR(u) private(bl, br)
   !$omp OMPLOOP DEVICEPTR(u) private(bl, br)
   do b=1, NB
      bl = b - 1 ; if (bl < 1) bl = NB
      br = b + 1 ; if (br > NB) br = 1
      if (periodic .or. b > 1) then
         u(1,0,1,1,b) = u(1,NI,1,1,bl)
      else
         u(1,0,1,1,b) = 1._R_P
      endif
      if (periodic .or. b < NB) then
         u(1,NI+1,1,1,b) = u(1,1,1,1,br)
      else
         u(1,NI+1,1,1,b) = u(1,NI,1,1,b)
      endif
   enddo
   endsubroutine fill_ghosts

   subroutine residual(u, r)
   !< The caller's residual: flux form, central face fluxes (u_L + u_R)/2.
   real(R_P), intent(in)    :: u(1:,0:,0:,0:,1:)
   real(R_P), intent(inout) :: r(1:,0:,0:,0:,1:)
   integer(I4P)             :: b, i

   !$acc parallel loop independent collapse(2) DEVICEVAR(u, r)
   !$omp OMPLOOP collapse(2) DEVICEPTR(u, r)
   do b=1, NB
      do i=1, NI
         r(1,i,1,1,b) = -(0.5_R_P*(u(1,i,1,1,b) + u(1,i+1,1,1,b)) - 0.5_R_P*(u(1,i-1,1,1,b) + u(1,i,1,1,b))) / DX
      enddo
   enddo
   endsubroutine residual

   function boundary_flux(u_d) result(f)
   !< F_right - F_left of the stage state: (u_N + ghost)/2 - (1 + u_1)/2.
   real(R_P), intent(in) :: u_d(1:,0:,0:,0:,1:)
   real(R_P)             :: f
   real(R_P)             :: h(1:1, 0:NI+1, 0:2, 0:2, 1:NB)

   call dev_memcpy_from_device(dst=h, src=u_d)
   f = 0.5_R_P*(h(1,NI,1,1,NB) + h(1,NI+1,1,1,NB)) - 0.5_R_P*(h(1,0,1,1,1) + h(1,1,1,1,1))
   endfunction boundary_flux

   function mass(q_d) result(m)
   !< Sum of the interior values times dx.
   real(R_P), intent(in) :: q_d(1:,0:,0:,0:,1:)
   real(R_P)             :: m
   real(R_P)             :: h(1:1, 0:NI+1, 0:2, 0:2, 1:NB)

   call dev_memcpy_from_device(dst=h, src=q_d)
   m = sum(h(1,1:NI,1,1,:)) * DX
   endfunction mass

   subroutine sentinel_registers(st)
   !< Fill every register with the sentinel (the j/k ghost layers must keep it: nobody writes them).
   type(stepper), intent(inout) :: st
   real(R_P), allocatable       :: h(:,:,:,:,:,:)

   allocate(h(lbound(st%regs,1):ubound(st%regs,1), lbound(st%regs,2):ubound(st%regs,2), lbound(st%regs,3): &
              ubound(st%regs,3), lbound(st%regs,4):ubound(st%regs,4), lbound(st%regs,5):ubound(st%regs,5), &
              lbound(st%regs,6):ubound(st%regs,6)))
   h = SENTINEL
   call dev_memcpy_to_device(dst=st%regs, src=h)
   endsubroutine sentinel_registers

   logical function registers_untouched(st) result(ok)
   !< The j/k ghost layers of every register still hold the sentinel.
   type(stepper), intent(in) :: st
   real(R_P), allocatable    :: h(:,:,:,:,:,:)

   allocate(h(lbound(st%regs,1):ubound(st%regs,1), lbound(st%regs,2):ubound(st%regs,2), lbound(st%regs,3): &
              ubound(st%regs,3), lbound(st%regs,4):ubound(st%regs,4), lbound(st%regs,5):ubound(st%regs,5), &
              lbound(st%regs,6):ubound(st%regs,6)))
   call dev_memcpy_from_device(dst=h, src=st%regs)
   ok = all(h(:,:,0,:,:,:) == SENTINEL) .and. all(h(:,:,2,:,:,:) == SENTINEL)
   endfunction registers_untouched

   subroutine split_scheme(name)
   !< A splitting through the stepper: q' = p (partition 1, component 1), p' = -q (partition 2, component 2).
   character(*), intent(in)    :: name
   type(splitting_composition) :: sc
   real(R_P)                   :: errors(4), orders(3)
   integer(I4P)                :: k, n

   sc = splitting(name)
#if !defined _R_P_IS_R4P
   n = merge(2_I4P, 8_I4P, sc%order >= 6)
   do k=1, 4
      errors(k) = oscillate(name, n * 2_I4P**(k-1))
   enddo
   orders = log(errors(1:3) / errors(2:4)) / log(2._R_P)
   print '(A,3F7.3)', '        '//name//' (stepper) observed orders:', orders
   call check(orders(3) >= real(sc%order, R_P) - 0.2_R_P, name//': observed order through the stepper')
#else
   errors(1) = oscillate(name, 4_I4P)
   call check(errors(1) < 1.e-2_R_P, name//': runs through the stepper (R4P: no order check)')
#endif
   endsubroutine split_scheme

   function oscillate(name, nsteps) result(err)
   !< nsteps steps to t = 1 of the oscillator in every cell (q0 = 0.1 b, p0 = 0.5); max error.
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: nsteps
   real(R_P)                :: err
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: q_d(:,:,:,:,:)=>null(), r_d(:,:,:,:,:)=>null()
   real(R_P)                :: h(1:2, 0:2, 0:2, 0:2, 1:NB), dt, q0, p0
   integer(I4P)             :: step, k, b, ierr

   bx = box_interior(LAYOUT_V_FIRST, 2_I4P, 1_I4P, 1_I4P, 1_I4P, NGC, NB)
   h = 0._R_P
   do b=1, NB
      h(1,1,1,1,b) = 0.1_R_P * b ; h(2,1,1,1,b) = 0.5_R_P
   enddo
   call dev_alloc(fptr_dev=q_d, ubounds=[2_I4P, 2_I4P, 2_I4P, 2_I4P, NB], lbounds=[1_I4P, 0_I4P, 0_I4P, 0_I4P, 1_I4P], &
                  ierr=ierr)
   call dev_memcpy_to_device(dst=q_d, src=h)
   call st%initialize(name, FORM_SPLITTING, lb=[1_I4P, 0_I4P, 0_I4P, 0_I4P, 1_I4P], ub=[2_I4P, 2_I4P, 2_I4P, 2_I4P, NB], &
                      partitions=reshape([1_I4P, 1_I4P, 2_I4P, 2_I4P], [2, 2]))
   dt = 1._R_P / real(nsteps, R_P)
   do step=1, nsteps
      call st%begin_step(bx, q_d)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, q_d)
         r_d => st%rhs_register()
         call osc_residual(q_d, r_d, st%stage_partition(k))
         call st%end_stage(k, bx, dt, q_d)
      enddo
      call st%finish_step(bx, dt, q_d)
   enddo
   call dev_memcpy_from_device(dst=h, src=q_d)
   err = 0._R_P
   do b=1, NB
      q0 = 0.1_R_P * b ; p0 = 0.5_R_P
      err = max(err, abs(h(1,1,1,1,b) - (q0*cos(1._R_P) + p0*sin(1._R_P))), &
                     abs(h(2,1,1,1,b) - (-q0*sin(1._R_P) + p0*cos(1._R_P))))
   enddo
   call st%destroy
   call dev_free(q_d)
   endfunction oscillate

   subroutine osc_residual(q, r, part)
   !< The residual of one partition: R_1 = p (component 1), R_2 = -q (component 2).
   real(R_P),    intent(in)    :: q(1:,0:,0:,0:,1:)
   real(R_P),    intent(inout) :: r(1:,0:,0:,0:,1:)
   integer(I4P), intent(in)    :: part
   integer(I4P)                :: b

   !$acc parallel loop independent DEVICEVAR(q, r)
   !$omp OMPLOOP DEVICEPTR(q, r)
   do b=1, NB
      if (part == 1) then
         r(1,1,1,1,b) = q(2,1,1,1,b)
      else
         r(2,1,1,1,b) = -q(1,1,1,1,b)
      endif
   enddo
   endsubroutine osc_residual
endprogram selah_test_stepper
