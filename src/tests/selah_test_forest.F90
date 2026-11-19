!< SELAH test: the stepper driven as ADAM's forest drives its integrator (the M7 precondition, PLAN §7).
program selah_test_forest
!< SELAH test: the stepper driven as ADAM's forest drives its integrator (the M7 precondition, PLAN §7).
!<
!< Two realms, each with its own stepper over a rank-5 field (v, i, j, k, b) of two blocks of N cells, make a periodic
!< line of four blocks: realm 1 blocks 1-2, realm 2 blocks 3-4. u_t + u_x = 0, the residual GIDEON's fourth-order
!< finite-volume derivative, the face fluxes GIDEON's fourth-order face values. Each step follows ADAM's forest
!< (adam_realm_object, adam_flume_cpu_object): open the step (the alpha fill writes the ghost cells of every q from its
!< neighbours' committed state), begin_stage(k) on every realm, then the seam fill of every realm's stage state, then
!< the residuals, then end_stage(k) on every realm; finish_step. The seams between realms are either filled at every
!< stage from the neighbour's stage state, or lagged (ADAM's end_of_step cadence, issue #113): every stage reads the
!< neighbour's state at the start of the step, which `seam_ghosts_to_stages` copies from the ghosts of q into the stage
!< buffers. Here the stage state is q itself in every form but Butcher, whose ghosts keep the alpha fill (the stepper
!< never writes a ghost cell); in the Butcher form the caller copies them into the stage state at every stage not on q.
!< Every Runge–Kutta scheme in the Butcher form (with and without stages on q) and in its other forms; asserted:
!<
!<+ **Per-stage seams.** The four blocks give the single-block solution of the same line (one stepper, one block of 4N
!<  cells), bitwise.
!<+ **Lagged seams.** The solution of the lagged scheme defined directly from the Butcher tableau (every stage stored,
!<  the lagged ghosts the neighbour's u^n, as ADAM's q_rk), to 1000 eps; and not the single-block one (the lag acts),
!<  but for a one-stage scheme, whose only stage state is u^n: there the lag is no lag (bitwise the single block).
!<+ **Reflux register.** sum_k stage_weight(k) F_k at a block's two faces telescopes its mass change, every step, every
!<  block, to round-off. With lagged seams the two sides of a seam register different fluxes and the mass drifts (more
!<  than one stage); the
!<  end-of-step reflux (ADAM's Berger–Colella correction, at equal resolution: one side takes the other's register)
!<  conserves it to round-off; so does the per-stage reflux (the caller corrects rhs_register() before end_stage).
!<+ **Registers as scratch (the D_flux diagnostic).** With the caller's registers, writing them between steps, as ADAM's
!<  D_flux diagnostic does (slot 1 a copy of q, slot 2 its residual), or filling them with NaN, leaves the solution
!<  bitwise unchanged.
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
use selah
use selah_test_tools, only : check, failures
use gideon, only : gideon_fv_d1_c4, gideon_fv_face_c4
implicit none
integer(I4P), parameter :: N = 16, G = 2, NB = 2, NR = 2, NBT = NB*NR, M = N*NBT, NSTEPS = 10
real(R_P),    parameter :: DX = 1._R_P / M, DT = 0.3_R_P * DX
real(R_P),    parameter :: TOL = 1000._R_P * epsilon(1._R_P)
character(7),  parameter :: OTHER_SCHEMES(6) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, &
                                                SCHEME_LSRK_33, SCHEME_LSRK_54, SCHEME_SSP_104]
character(14), parameter :: OTHER_FORMS(6) = [character(14) :: FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, &
                                              FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_WILLIAMSON_2N, &
                                              FORM_KETCHESON_2R]
integer(I4P), parameter :: REFLUX_NONE = 0, REFLUX_STEP = 1, REFLUX_STAGE = 2
integer(I4P), parameter :: SCRATCH_NONE = 0, SCRATCH_DFLUX = 1, SCRATCH_NAN = 2

type :: realm_t
   !< A realm: its field, its stepper, the registers when the caller owns them, the face fluxes of its blocks.
   real(R_P), allocatable :: q(:,:,:,:,:)            !< Field (1, 1-G:N+G, 1-G:1+G, 1-G:1+G, NB).
   real(R_P), allocatable :: qn(:,:,:,:,:)           !< The field at the start of the step.
   type(stepper)          :: st                      !< Stepper.
   real(R_P), pointer     :: regs(:,:,:,:,:,:) => null() !< Caller registers (or null).
   real(R_P), pointer     :: u(:,:,:,:,:) => null()  !< Stage state of the active stage.
   logical                :: u_is_q = .true.         !< The stage state is q itself (a flag, not associated(u, q):
                                                     !< nvfortran 26.9's front end crashes on that, see the night log).
   real(R_P), pointer     :: r(:,:,:,:,:) => null()  !< Residual register of the active stage.
   real(R_P)              :: fst(2,NB) = 0._R_P      !< Face fluxes (left, right) of the active stage.
   real(R_P)              :: freg(2,NB) = 0._R_P     !< Register: sum_k stage_weight(k) fst over the step.
endtype realm_t

real(R_P)    :: u0(M)
integer(I4P) :: s, c0

do c0=1, M
   u0(c0) = sin(6.283185307179586_R_P * (c0 - 0.5_R_P) * DX) + 0.5_R_P * exp(-100._R_P * ((c0 - 0.5_R_P) * DX - 0.5_R_P)**2)
enddo
do s=1, size(RK_SCHEMES)
   call configuration(trim(RK_SCHEMES(s)), FORM_BUTCHER, .false.)
   call configuration(trim(RK_SCHEMES(s)), FORM_BUTCHER, .true.)
enddo
do s=1, size(OTHER_SCHEMES)
   call configuration(trim(OTHER_SCHEMES(s)), trim(OTHER_FORMS(s)), .false.)
enddo

if (failures > 0) error stop 'selah_test_forest: failed'
print '(A)', 'selah_test_forest: ok'

contains
   subroutine configuration(name, form, onq)
   !< Every assertion for one scheme, form and stages-on-q choice.
   character(*), intent(in) :: name, form
   logical,      intent(in) :: onq
   character(:), allocatable :: tag
   real(R_P)    :: qs(M), qf(M), qref(M), qa(M), qb(M), tele, drift, mass0
   real(R_P)    :: tele_l, drift_none, drift_step, drift_stage

   tag = name//' '//form//merge(' stages_on_q', '            ', onq)
   mass0 = sum(u0) * DX
   call single_block(name, form, onq, qs)
   call forest(name, form, onq, .false., REFLUX_NONE, SCRATCH_NONE, qf, tele, drift)
   call check(all(qf == qs), tag//': per-stage seams, the single-block solution bitwise')
   call reference(name, qref)
   call forest(name, form, onq, .true., REFLUX_NONE, SCRATCH_NONE, qa, tele_l, drift_none)
   print '(A,ES10.2,A,ES10.2,A,ES10.2)', '        '//tag//': lagged vs reference', maxval(abs(qa - qref)), &
         ', vs single block', maxval(abs(qa - qs)), ', telescoping', max(tele, tele_l)
   call check(maxval(abs(qa - qref)) <= TOL * maxval(abs(qref)), tag//': lagged seams, the reference solution')
   if (rl_stages(name, form) > 1) then
      ! the lag's effect is physical (1.6e-2 of max|q| at the least), far above the round-off of R4P builds
      call check(maxval(abs(qa - qs)) > 1.e-3_R_P * maxval(abs(qs)), tag//': lagged seams act (not the single-block solution)')
   else ! one stage: its state is u^n, so a lagged seam reads what a per-stage seam reads
      call check(all(qa == qs), tag//': one stage, lagged seams give the single-block solution bitwise')
   endif
   call check(max(tele, tele_l) <= TOL * mass0, tag//': the register telescopes every block''s mass change')
   call forest(name, form, onq, .true., REFLUX_STEP, SCRATCH_NONE, qb, tele, drift_step)
   call forest(name, form, onq, .true., REFLUX_STAGE, SCRATCH_NONE, qb, tele, drift_stage)
   print '(A,3ES10.2)', '        '//tag//': mass drift, no reflux / end of step / per stage', drift_none, drift_step, &
         drift_stage
   if (rl_stages(name, form) > 1) then
      ! the drift is physical (5e-4 of the mass at the least), far above the round-off of R4P builds
      call check(drift_none > 1.e-4_R_P * mass0, tag//': lagged seams without reflux do not conserve')
   else
      call check(drift_none <= TOL * mass0, tag//': one stage, lagged seams conserve without reflux')
   endif
   call check(drift_step <= TOL * mass0, tag//': end-of-step reflux conserves the mass')
   call check(drift_stage <= TOL * mass0, tag//': per-stage reflux conserves the mass')
   call forest(name, form, onq, .true., REFLUX_STEP, SCRATCH_DFLUX, qb, tele, drift)
   call forest(name, form, onq, .true., REFLUX_STEP, SCRATCH_NONE, qa, tele, drift, own=.false.)
   call check(all(qb == qa), tag//': caller registers written between steps (D_flux), the solution bitwise')
   call forest(name, form, onq, .true., REFLUX_STEP, SCRATCH_NAN, qb, tele, drift)
   call check(all(qb == qa), tag//': caller registers NaN-filled between steps, the solution bitwise')
   call forest(name, form, onq, .true., REFLUX_STEP, SCRATCH_NONE, qb, tele, drift)
   call check(all(qb == qa), tag//': caller registers, the solution of the own registers bitwise')
   endsubroutine configuration

   subroutine forest(name, form, onq, lagged, reflux, scratch, qout, tele, drift, own)
   !< NSTEPS steps of the two realms, as ADAM's forest drives them; qout the line of 4N cells at the end, tele the
   !< largest telescoping mismatch, drift the largest |mass - mass(0)|.
   character(*),      intent(in)  :: name, form
   logical,           intent(in)  :: onq     !< Butcher form: stages on q.
   logical,           intent(in)  :: lagged  !< The seams between realms are lagged (end_of_step cadence).
   integer(I4P),      intent(in)  :: reflux  !< REFLUX_NONE, REFLUX_STEP, REFLUX_STAGE.
   integer(I4P),      intent(in)  :: scratch !< SCRATCH_NONE, SCRATCH_DFLUX, SCRATCH_NAN (caller registers).
   real(R_P),         intent(out) :: qout(M), tele, drift
   logical, optional, intent(in)  :: own     !< .false.: caller registers even without scratch (default: own
                                             !< registers unless scratch is asked).
   type(realm_t), target :: rl(NR)
   type(box)             :: bx
   real(R_P)             :: mass0, mass
   integer(I4P)          :: r, b, k, step
   logical               :: caller

   caller = scratch /= SCRATCH_NONE
   if (present(own)) caller = caller .or. .not. own
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, N, 1_I4P, 1_I4P, G, NB)
   do r=1, NR
      allocate(rl(r)%q(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NB), rl(r)%qn(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NB))
      rl(r)%q = 0._R_P
      do b=1, NB
         rl(r)%q(1,1:N,1,1,b) = u0(((r-1)*NB + b - 1)*N + 1 : ((r-1)*NB + b)*N)
      enddo
      if (caller) then
         allocate(rl(r)%regs(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NB, 1:registers_needed(name, form)))
         rl(r)%regs = ieee_value(1._R_P, ieee_quiet_nan)
         call rl(r)%st%initialize(name, form, lb=lbound(rl(r)%q), ub=ubound(rl(r)%q), registers=rl(r)%regs, &
                                  stages_on_q=onq)
      else
         call rl(r)%st%initialize(name, form, lb=lbound(rl(r)%q), ub=ubound(rl(r)%q), stages_on_q=onq)
      endif
   enddo
   mass0 = line_mass(rl) ; tele = 0._R_P ; drift = 0._R_P
   do step=1, NSTEPS
      ! open the step: the alpha fill of every q from the neighbours' committed state, every seam
      do r=1, NR
         rl(r)%qn = rl(r)%q
         rl(r)%freg = 0._R_P
         rl(r)%u => rl(r)%q
         rl(r)%u_is_q = .true.
      enddo
      call seam_fill(rl, .false., .true.)
      do r=1, NR
         call rl(r)%st%begin_step(bx, rl(r)%q)
      enddo
      do k=1, rl(1)%st%stages()
         do r=1, NR ! begin_stage on every realm before any seam fill
            call rl(r)%st%begin_stage(k, bx, DT, rl(r)%q)
            rl(r)%u_is_q = rl(r)%st%stage_state_is_q(k)
            if (rl(r)%u_is_q) then
               rl(r)%u => rl(r)%q
            else
               rl(r)%u => rl(r)%st%stage_state()
            endif
         enddo
         call seam_fill(rl, lagged, .false.)
         do r=1, NR
            rl(r)%r => rl(r)%st%rhs_register()
            call residual(N, NB, rl(r)%u, rl(r)%r, rl(r)%fst)
         enddo
         if (reflux == REFLUX_STAGE .and. lagged) call reflux_stage_correction(rl)
         do r=1, NR
            rl(r)%freg = rl(r)%freg + rl(r)%st%stage_weight(k) * rl(r)%fst
            call rl(r)%st%end_stage(k, bx, DT, rl(r)%q)
         enddo
      enddo
      do r=1, NR
         call rl(r)%st%finish_step(bx, DT, rl(r)%q)
         do b=1, NB ! the register telescopes the block's mass change (before any reflux)
            tele = max(tele, abs(sum(rl(r)%q(1,1:N,1,1,b) - rl(r)%qn(1,1:N,1,1,b)) * DX + &
                                 DT * (rl(r)%freg(2,b) - rl(r)%freg(1,b))))
         enddo
      enddo
      if (reflux == REFLUX_STEP .and. lagged) call reflux_step_correction(rl)
      mass = line_mass(rl)
      drift = max(drift, abs(mass - mass0))
      select case(scratch) ! between steps: the registers are the caller's scratch
      case(SCRATCH_DFLUX) ! ADAM's D_flux diagnostic: slot 1 a copy of q, slot 2 its residual
         do r=1, NR
            rl(r)%regs(:,:,:,:,:,1) = rl(r)%q
            call residual(N, NB, rl(r)%regs(:,:,:,:,:,1), rl(r)%regs(:,:,:,:,:,2), rl(r)%fst)
         enddo
      case(SCRATCH_NAN)
         do r=1, NR
            rl(r)%regs = ieee_value(1._R_P, ieee_quiet_nan)
         enddo
      endselect
   enddo
   do r=1, NR
      do b=1, NB
         qout(((r-1)*NB + b - 1)*N + 1 : ((r-1)*NB + b)*N) = rl(r)%q(1,1:N,1,1,b)
      enddo
      call rl(r)%st%destroy
      if (associated(rl(r)%regs)) deallocate(rl(r)%regs)
   enddo
   endsubroutine forest

   subroutine seam_fill(rl, lagged, alpha)
   !< The ghost cells along i of every block. alpha: the step's opening fill of q from the neighbours' q. Otherwise
   !< the stage fill of every stage state: from the neighbour's stage state, but for a lagged seam (between realms)
   !< from the ghosts of the realm's own q (seam_ghosts_to_stages), nothing to do when the stage state is q.
   type(realm_t), intent(inout), target :: rl(NR)
   logical,       intent(in)            :: lagged, alpha
   integer(I4P) :: r, b, gb, gl, gr, rn, bn

   do r=1, NR
      do b=1, NB
         gb = (r - 1)*NB + b
         gl = modulo(gb - 2, NBT) + 1 ; gr = modulo(gb, NBT) + 1
         ! left ghosts from the left neighbour's last G cells
         rn = (gl - 1)/NB + 1 ; bn = gl - (rn - 1)*NB
         if (alpha) then
            rl(r)%q(1,1-G:0,1,1,b) = rl(rn)%q(1,N-G+1:N,1,1,bn)
         elseif (lagged .and. rn /= r) then
            if (.not. rl(r)%u_is_q) rl(r)%u(1,1-G:0,1,1,b) = rl(r)%q(1,1-G:0,1,1,b)
         else
            rl(r)%u(1,1-G:0,1,1,b) = rl(rn)%u(1,N-G+1:N,1,1,bn)
         endif
         ! right ghosts from the right neighbour's first G cells
         rn = (gr - 1)/NB + 1 ; bn = gr - (rn - 1)*NB
         if (alpha) then
            rl(r)%q(1,N+1:N+G,1,1,b) = rl(rn)%q(1,1:G,1,1,bn)
         elseif (lagged .and. rn /= r) then
            if (.not. rl(r)%u_is_q) rl(r)%u(1,N+1:N+G,1,1,b) = rl(r)%q(1,N+1:N+G,1,1,b)
         else
            rl(r)%u(1,N+1:N+G,1,1,b) = rl(rn)%u(1,1:G,1,1,bn)
         endif
      enddo
   enddo
   endsubroutine seam_fill

   subroutine reflux_step_correction(rl)
   !< End-of-step reflux at every seam between realms: the block on the right of the seam takes the register of the
   !< block on its left (equal resolution: the Berger–Colella correction dt/dx (F_other - F_own) of its first cell).
   type(realm_t), intent(inout) :: rl(NR)
   integer(I4P) :: r, b, gb, gl, rn, bn

   do r=1, NR
      do b=1, NB
         gb = (r - 1)*NB + b ; gl = modulo(gb - 2, NBT) + 1
         rn = (gl - 1)/NB + 1 ; bn = gl - (rn - 1)*NB
         if (rn == r) cycle
         rl(r)%q(1,1,1,1,b) = rl(r)%q(1,1,1,1,b) + DT / DX * (rl(rn)%freg(2,bn) - rl(r)%freg(1,b))
      enddo
   enddo
   endsubroutine reflux_step_correction

   subroutine reflux_stage_correction(rl)
   !< Per-stage reflux at every seam between realms: the residual of the first cell on the right of the seam is
   !< corrected to the left block's face flux of this stage, before end_stage; its register follows.
   type(realm_t), intent(inout) :: rl(NR)
   integer(I4P) :: r, b, gb, gl, rn, bn

   do r=1, NR
      do b=1, NB
         gb = (r - 1)*NB + b ; gl = modulo(gb - 2, NBT) + 1
         rn = (gl - 1)/NB + 1 ; bn = gl - (rn - 1)*NB
         if (rn == r) cycle
         rl(r)%r(1,1,1,1,b) = rl(r)%r(1,1,1,1,b) + (rl(rn)%fst(2,bn) - rl(r)%fst(1,b)) / DX
         rl(r)%fst(1,b) = rl(rn)%fst(2,bn)
      enddo
   enddo
   endsubroutine reflux_stage_correction

   function rl_stages(name, form) result(n_)
   !< Stages of a step of the scheme in the form.
   character(*), intent(in) :: name, form
   integer(I4P)             :: n_
   type(stepper)            :: st

   call st%initialize(name, form, lb=[1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P], ub=[1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P])
   n_ = st%stages()
   call st%destroy
   endfunction rl_stages

   function line_mass(rl) result(mass)
   !< sum u dx over the four blocks.
   type(realm_t), intent(in) :: rl(NR)
   real(R_P)                 :: mass
   integer(I4P)              :: r

   mass = 0._R_P
   do r=1, NR
      mass = mass + sum(rl(r)%q(1,1:N,1,1,:)) * DX
   enddo
   endfunction line_mass

   subroutine residual(n_, nb_, u, r, f)
   !< r = -du/dx (GIDEON's fourth-order finite-volume derivative) in the cells; f the faces' fluxes (left, right) of
   !< every block (GIDEON's fourth-order face values, the wave speed 1).
   integer(I4P), intent(in)    :: n_, nb_
   real(R_P),    intent(in)    :: u(1:,1-G:,1-G:,1-G:,1:)
   real(R_P),    intent(inout) :: r(1:,1-G:,1-G:,1-G:,1:)
   real(R_P),    intent(out)   :: f(1:,1:)
   real(R_P)                   :: d
   integer(I4P)                :: b, c

   do b=1, nb_
      do c=1, n_
         call gideon_fv_d1_c4(u(1,c-2:c+2,1,1,b), DX, d)
         r(1,c,1,1,b) = -d
      enddo
      call gideon_fv_face_c4(u(1,-1:2,1,1,b), f(1,b))
      call gideon_fv_face_c4(u(1,n_-1:n_+2,1,1,b), f(2,b))
   enddo
   endsubroutine residual

   subroutine single_block(name, form, onq, qout)
   !< The same line as one block of 4N cells, periodic, one stepper.
   character(*), intent(in)  :: name, form
   logical,      intent(in)  :: onq
   real(R_P),    intent(out) :: qout(M)
   real(R_P), allocatable, target :: q(:,:,:,:,:)
   real(R_P), pointer        :: u(:,:,:,:,:), r(:,:,:,:,:)
   real(R_P)                 :: f(2,1)
   type(stepper)             :: st
   type(box)                 :: bx
   integer(I4P)              :: step, k

   allocate(q(1:1, 1-G:M+G, 1-G:1+G, 1-G:1+G, 1:1)) ; q = 0._R_P
   q(1,1:M,1,1,1) = u0
   bx = box_interior(LAYOUT_V_FIRST, 1_I4P, M, 1_I4P, 1_I4P, G, 1_I4P)
   call st%initialize(name, form, lb=lbound(q), ub=ubound(q), stages_on_q=onq)
   do step=1, NSTEPS
      call st%begin_step(bx, q)
      do k=1, st%stages()
         call st%begin_stage(k, bx, DT, q)
         if (st%stage_state_is_q(k)) then
            u => q
         else
            u => st%stage_state()
         endif
         u(1,1-G:0,1,1,1) = u(1,M-G+1:M,1,1,1) ; u(1,M+1:M+G,1,1,1) = u(1,1:G,1,1,1)
         r => st%rhs_register()
         call residual(M, 1_I4P, u, r, f)
         call st%end_stage(k, bx, DT, q)
      enddo
      call st%finish_step(bx, DT, q)
   enddo
   qout = q(1,1:M,1,1,1)
   call st%destroy
   endsubroutine single_block

   subroutine reference(name, qout)
   !< The lagged scheme from its definition: every stage state stored (as ADAM's q_rk), its seam ghosts between realms
   !< the neighbour's u^n, the others the neighbour's stage state; u^{n+1} = u^n + dt sum_s b_s F_s.
   character(*), intent(in)  :: name
   real(R_P),    intent(out) :: qout(M)
   type(butcher_tableau)     :: bt
   real(R_P), allocatable    :: y(:,:,:,:,:,:), fs(:,:,:,:,:,:), qn(:,:,:,:,:)
   real(R_P)                 :: f(2,NBT)
   integer(I4P)              :: step, s, sm, b, bl, br

   bt = butcher(name)
   allocate(y(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NBT, bt%stages), fs(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NBT, bt%stages))
   allocate(qn(1:1, 1-G:N+G, 1-G:1+G, 1-G:1+G, 1:NBT))
   y = 0._R_P ; fs = 0._R_P ; qn = 0._R_P
   do b=1, NBT
      qn(1,1:N,1,1,b) = u0((b - 1)*N + 1 : b*N)
   enddo
   do step=1, NSTEPS
      do b=1, NBT ! the committed ghosts: the neighbours' u^n
         bl = modulo(b - 2, NBT) + 1 ; br = modulo(b, NBT) + 1
         qn(1,1-G:0,1,1,b) = qn(1,N-G+1:N,1,1,bl) ; qn(1,N+1:N+G,1,1,b) = qn(1,1:G,1,1,br)
      enddo
      do s=1, bt%stages
         y(:,:,:,:,:,s) = qn
         do sm=1, s-1
            y(1,1:N,1,1,:,s) = y(1,1:N,1,1,:,s) + DT * bt%a(s,sm) * fs(1,1:N,1,1,:,sm)
         enddo
         do b=1, NBT
            bl = modulo(b - 2, NBT) + 1 ; br = modulo(b, NBT) + 1
            if ((bl - 1)/NB == (b - 1)/NB) y(1,1-G:0,1,1,b,s) = y(1,N-G+1:N,1,1,bl,s) ! same realm: per stage
            if ((br - 1)/NB == (b - 1)/NB) y(1,N+1:N+G,1,1,b,s) = y(1,1:G,1,1,br,s)
         enddo
         call residual(N, NBT, y(:,:,:,:,:,s), fs(:,:,:,:,:,s), f)
      enddo
      do s=1, bt%stages
         qn(1,1:N,1,1,:) = qn(1,1:N,1,1,:) + DT * bt%b(s) * fs(1,1:N,1,1,:,s)
      enddo
   enddo
   do b=1, NBT
      qout((b - 1)*N + 1 : b*N) = qn(1,1:N,1,1,b)
   enddo
   endsubroutine reference
endprogram selah_test_forest
