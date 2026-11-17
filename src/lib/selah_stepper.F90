!< SELAH, L2: the stepper, a time step the caller drives one stage at a time.

#include "fundal.H"

module selah_stepper
!< SELAH, L2: the stepper, a time step the caller drives one stage at a time.
!<
!< The stepper owns the registers and the coefficients; the caller owns the solution `q`, the time loop, the ghost
!< cells and the residual. One step of any scheme is the same sequence of calls:
!<
!<```fortran
!< call st%begin_step(bx, q)
!< do k=1, st%stages()
!<    call st%begin_stage(k, bx, dt, q)                  ! builds the state of stage k
!<    ! fill the ghost cells of the stage state (q itself when st%stage_state_is_q(k), else st%stage_state())
!<    ! write its residual into st%rhs_register() (the components of st%stage_partition(k) for a splitting)
!<    call st%end_stage(k, bx, dt, q)
!< enddo
!< call st%finish_step(bx, dt, q)
!<```
!<
!< **The stage contract.**
!<
!< 1. Between `begin_stage(k)` and `end_stage(k)` the stage state holds the state of stage k in the cells of the box:
!<    `q` itself when `stage_state_is_q(k)` (every stage of every form but the Butcher one, where `state_aliases_q()`
!<    is true; in the Butcher form the stages whose row of A is zero, when the stepper was initialised with
!<    `stages_on_q=.true.`), else `stage_state()`.
!< 2. The stepper never writes a ghost cell, of any register or of `q`.
!< 3. The caller may write the ghost cells of the stage state during that window (ghost exchange, seam fill, physical
!<    boundaries).
!< 4. `rhs_register()` is written by the caller only, in the cells of the box, and is distinct from the stage state.
!< 5. After `finish_step` the registers are scratch.
!<
!< `stage_weight(k)` is the effective Butcher weight b_k of stage k in every Runge–Kutta form (the weight a flux
!< register needs to accumulate the step's boundary fluxes); for a splitting it is the fraction of the step of flow k.
!<
!< **Stages on q.** The state of a stage whose row of A is zero (the first stage of every explicit scheme) is q itself.
!< With `stages_on_q=.true.` the Butcher form does not copy q into its stage-state register for such a stage: the
!< caller evaluates the residual on q, and `stage_state()` is null during that stage, so that a caller that does not
!< ask `stage_state_is_q(k)` fails at once instead of reading a stale register. Default `.false.`: every stage of the
!< Butcher form has its register (v0.1.0 behaviour).
!<
!< **Storage.** `registers_needed(scheme, form)` field-sized buffers, allocated through FUNDAL (device memory in a
!< device build) or passed by the caller (`registers`, any rank-6 array with the field's bounds and enough slots).
use fundal, only : dev_alloc, dev_free, dev_assign_to_device
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_tableaux
use selah_schemes
use selah_splittings
use selah_storage
use selah_box
use selah_kernels
implicit none
private
public :: stepper

type :: stepper
   !< A time step driven stage by stage.
   character(:), allocatable   :: scheme                  !< Scheme name.
   character(:), allocatable   :: form                    !< Storage form.
   integer(I4P)                :: nstages = 0_I4P         !< Stages (flows for a splitting).
   integer(I4P)                :: k_active = 0_I4P        !< Stage between begin_stage and end_stage, else 0.
   integer(I4P)                :: lb(5) = 1_I4P           !< Lower bounds of the field.
   integer(I4P)                :: part(2,2) = 0_I4P       !< Component ranges of the two partitions (splitting).
   type(butcher_tableau)       :: bt                      !< Butcher tableau (every RK form).
   type(shu_osher_2r_form)     :: r2                      !< Two-register form.
   type(williamson_2n_form)    :: w2n                     !< 2N form.
   type(splitting_composition) :: sc                      !< Splitting.
   real(R_P), pointer          :: regs(:,:,:,:,:,:) => null() !< Registers.
   real(R_P), pointer          :: at_d(:,:) => null()     !< Transposed A (column k = row k), device resident.
   real(R_P), pointer          :: b_d(:) => null()        !< Weights b, device resident.
   logical                     :: owned = .false.         !< The registers are the stepper's.
   logical, allocatable        :: on_q(:)                 !< Stage k's state is q itself.
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: destroy
      procedure, pass(self) :: begin_step
      procedure, pass(self) :: begin_stage
      procedure, pass(self) :: end_stage
      procedure, pass(self) :: finish_step
      procedure, pass(self) :: stage_state
      procedure, pass(self) :: rhs_register
      procedure, pass(self) :: stages
      procedure, pass(self) :: stage_time
      procedure, pass(self) :: stage_weight
      procedure, pass(self) :: stage_partition
      procedure, pass(self) :: state_aliases_q
      procedure, pass(self) :: stage_state_is_q
endtype stepper

contains
   subroutine initialize(self, scheme, form, lb, ub, partitions, registers, stages_on_q)
   !< Set the scheme and its storage form; allocate the registers (or take the caller's).
   class(stepper),     intent(inout)                  :: self          !< Stepper.
   character(*),       intent(in)                     :: scheme        !< Scheme name (selah_schemes, selah_splittings).
   character(*),       intent(in)                     :: form          !< Storage form (selah_storage).
   integer(I4P),       intent(in)                     :: lb(5), ub(5)  !< Bounds of the field, ghosts included.
   integer(I4P),       intent(in), optional           :: partitions(2,2) !< Splitting: (v1, v2) of partitions 1, 2.
   real(R_P),          intent(in), optional, pointer  :: registers(:,:,:,:,:,:) !< The caller's registers.
   logical,            intent(in), optional           :: stages_on_q   !< Butcher form: the stages whose row of A is
                                                                       !< zero are evaluated on q (default .false.).
   integer(I4P)                                       :: nreg, ierr  !< Registers needed, error status.
   integer(I4P)                                       :: k           !< Stage.

   call self%destroy
   self%scheme = trim(adjustl(scheme))
   self%form = trim(adjustl(form))
   self%lb = lb
   nreg = registers_needed(self%scheme, self%form)
   if (nreg < 0) error stop 'selah stepper: scheme "'//self%scheme//'" has no form "'//self%form//'"'
   self%bt = butcher(self%scheme)
   select case(self%form)
   case(FORM_BUTCHER)
      self%nstages = self%bt%stages
      call dev_assign_to_device(dst=self%at_d, src=transpose(self%bt%a))
      call dev_assign_to_device(dst=self%b_d, src=self%bt%b)
   case(FORM_SHU_OSHER_2R)
      self%r2 = shu_osher_2r(self%scheme)
      self%nstages = self%r2%stages
   case(FORM_WILLIAMSON_2N)
      self%w2n = williamson_2n(self%scheme)
      self%nstages = self%w2n%stages
   case(FORM_KETCHESON_2R)
      self%nstages = self%bt%stages
   case(FORM_SPLITTING)
      if (.not. present(partitions)) error stop 'selah stepper: a splitting needs the components of its partitions'
      self%sc = splitting(self%scheme)
      self%nstages = self%sc%flows
      self%part = partitions
   endselect
   allocate(self%on_q(self%nstages))
   if (self%form == FORM_BUTCHER) then
      self%on_q = .false.
      if (present(stages_on_q)) then
         if (stages_on_q) then
            do k=1, self%nstages
               self%on_q(k) = all(self%bt%a(k,:) == 0._R_P)
            enddo
         endif
      endif
   else
      self%on_q = .true.
   endif
   if (present(registers)) then
      if (.not. associated(registers)) error stop 'selah stepper: the registers passed are not associated'
      if (size(registers, 6) < nreg) error stop 'selah stepper: too few registers passed'
      if (lbound(registers, 6) /= 1) error stop 'selah stepper: the registers passed must be numbered from 1 (dim 6)'
      if (any(lbound(registers) /= [lb, 1_I4P]) .or. any(ubound(registers, kind=I4P) /= [ub, ubound(registers, 6)])) &
         error stop 'selah stepper: the registers passed do not have the bounds of the field (lb, ub)'
      self%regs => registers
      self%owned = .false.
   else
      ! zeroed, so that every register starts finite (the 2N form's first stage no longer reads Delta Q: A_1 = 0)
      call dev_alloc(fptr_dev=self%regs, ubounds=[ub, nreg], lbounds=[lb, 1_I4P], init_value=0._R_P, ierr=ierr)
      if (ierr /= 0) error stop 'selah stepper: the registers could not be allocated'
      self%owned = .true.
   endif
   endsubroutine initialize

   subroutine destroy(self)
   !< Free what the stepper owns.
   class(stepper), intent(inout) :: self !< Stepper.

   if (self%owned .and. associated(self%regs)) call dev_free(self%regs)
   if (associated(self%at_d)) call dev_free(self%at_d)
   if (associated(self%b_d)) call dev_free(self%b_d)
   nullify(self%regs, self%at_d, self%b_d)
   if (allocated(self%on_q)) deallocate(self%on_q)
   self%owned = .false. ; self%nstages = 0_I4P ; self%k_active = 0_I4P
   endsubroutine destroy

   subroutine begin_step(self, bx, q, mask, mask_comp)
   !< Open a step: keep what the form needs of the solution at the start of the step.
   class(stepper), intent(inout)        :: self                 !< Stepper.
   type(box),      intent(in)           :: bx                   !< Box.
   real(R_P),      intent(in)           :: q(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Solution.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp            !< Mask component.

   select case(self%form)
   case(FORM_SHU_OSHER_2R, FORM_KETCHESON_2R) ! u^n (two-register) or q2 (Ketcheson) start as the solution
      call selah_copy(bx, q, self%regs(:,:,:,:,:,1), mask, mask_comp)
   endselect
   self%k_active = 0_I4P
   endsubroutine begin_step

   subroutine begin_stage(self, k, bx, dt, q, mask, mask_comp)
   !< Open stage k: build its state (Butcher form, unless the stage's state is q; every other form keeps it in q).
   class(stepper), intent(inout)        :: self                 !< Stepper.
   integer(I4P),   intent(in)           :: k                    !< Stage.
   type(box),      intent(in)           :: bx                   !< Box.
   real(R_P),      intent(in)           :: dt                   !< Time step.
   real(R_P),      intent(in)           :: q(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Solution.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp            !< Mask component.

   if (self%form == FORM_BUTCHER .and. .not. self%on_q(k)) &
      call selah_stage_state(bx, k-1_I4P, dt, self%at_d(:,k), q, self%regs(:,:,:,:,:,2:), self%regs(:,:,:,:,:,1), &
                             mask, mask_comp)
   self%k_active = k
   endsubroutine begin_stage

   subroutine end_stage(self, k, bx, dt, q, mask, mask_comp)
   !< Close stage k, its residual written into rhs_register(): the forms that advance q in place do it here.
   class(stepper), intent(inout)        :: self                 !< Stepper.
   integer(I4P),   intent(in)           :: k                    !< Stage.
   type(box),      intent(in)           :: bx                   !< Box.
   real(R_P),      intent(in)           :: dt                   !< Time step.
   real(R_P),      intent(inout)        :: q(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Solution.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp            !< Mask component.
   integer(I4P)                         :: p                    !< Partition.

   associate(r1=>self%regs(:,:,:,:,:,1))
   select case(self%form)
   case(FORM_SHU_OSHER_2R) ! q = a_k u^n + b_k q + dt c_k R
      call selah_shu_osher_2r_stage(bx, dt, self%r2%a(k), self%r2%b(k), self%r2%c(k), r1, self%regs(:,:,:,:,:,2), q, &
                                    mask, mask_comp)
   case(FORM_WILLIAMSON_2N) ! dQ = A_k dQ + dt R; q = q + B_k dQ
      call selah_williamson_2n_stage(bx, dt, self%w2n%a(k), self%w2n%b(k), self%regs(:,:,:,:,:,2), r1, q, mask, mask_comp)
   case(FORM_KETCHESON_2R) ! Ketcheson (2008): Euler steps of dt/6, the mid-step combination, the last stage
      if (k < 10_I4P) call selah_axpy(bx, dt/6._R_P, self%regs(:,:,:,:,:,2), q, mask, mask_comp)
      if (k == 5_I4P) then ! q2 = (u^n + 9 q1)/25; q1 = 15 q2 - 5 q1
         call selah_lincomb(bx, 9._R_P/25._R_P, q, 1._R_P/25._R_P, r1, mask, mask_comp)
         call selah_lincomb(bx, 15._R_P, r1, -5._R_P, q, mask, mask_comp)
      endif
      if (k == 10_I4P) & ! u^{n+1} = q2 + 3/5 q1 + dt/10 F(q1)
         call selah_shu_osher_2r_stage(bx, dt, 1._R_P, 3._R_P/5._R_P, 1._R_P/10._R_P, r1, self%regs(:,:,:,:,:,2), q, &
                                       mask, mask_comp)
   case(FORM_SPLITTING) ! the components of the flow's partition advance by c_k dt R
      p = self%sc%partition(k)
      call selah_axpy(box_components(bx, self%part(1,p), self%part(2,p)), self%sc%coefficient(k) * dt, r1, q, &
                      mask, mask_comp)
   endselect
   endassociate
   self%k_active = 0_I4P
   endsubroutine end_stage

   subroutine finish_step(self, bx, dt, q, dq, mask, mask_comp)
   !< Close the step: the Butcher form advances q now (with the step residual dq if asked); the others already have.
   class(stepper), intent(inout)           :: self               !< Stepper.
   type(box),      intent(in)              :: bx                 !< Box.
   real(R_P),      intent(in)              :: dt                 !< Time step.
   real(R_P),      intent(inout)           :: q(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):)  !< Solution.
   real(R_P),      intent(inout), optional :: dq(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Step residual.
   real(R_P),      intent(in),    optional :: mask(bx%lb(1):,bx%lb(2):,bx%lb(3):,bx%lb(4):,bx%lb(5):) !< Mask.
   integer(I4P),   intent(in),    optional :: mask_comp          !< Mask component.

   if (self%form == FORM_BUTCHER) then
      call selah_step_update(bx, self%nstages, dt, self%b_d, self%regs(:,:,:,:,:,2:), q, dq, mask, mask_comp)
   elseif (present(dq)) then
      error stop 'selah stepper: the step residual dq is computed by the Butcher form only'
   endif
   endsubroutine finish_step

   function stage_state(self) result(u)
   !< The state of the active stage, with the field's bounds: the Butcher form's register; null when the stage's state
   !< is q itself (stage_state_is_q), so in every other form.
   class(stepper), intent(in) :: self         !< Stepper.
   real(R_P),      pointer    :: u(:,:,:,:,:) !< Stage state.

   u => null()
   if (self%k_active > 0_I4P) then
      if (self%on_q(self%k_active)) return
   endif
   if (self%form == FORM_BUTCHER) &
      u(self%lb(1):,self%lb(2):,self%lb(3):,self%lb(4):,self%lb(5):) => self%regs(:,:,:,:,:,1)
   endfunction stage_state

   function rhs_register(self) result(r)
   !< Where the caller writes the residual of the active stage, with the field's bounds.
   class(stepper), intent(in) :: self         !< Stepper.
   real(R_P),      pointer    :: r(:,:,:,:,:) !< Residual register.
   integer(I4P)               :: slot         !< Register slot.

   select case(self%form)
   case(FORM_BUTCHER)
      slot = self%k_active + 1_I4P
   case(FORM_SPLITTING)
      slot = 1_I4P
   case default
      slot = 2_I4P
   endselect
   r(self%lb(1):,self%lb(2):,self%lb(3):,self%lb(4):,self%lb(5):) => self%regs(:,:,:,:,:,slot)
   endfunction rhs_register

   pure function stages(self) result(n)
   !< Stages of a step (flows of a splitting).
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P)               :: n    !< Stages.

   n = self%nstages
   endfunction stages

   pure function stage_time(self, k) result(c)
   !< Abscissa c_k of stage k (fraction of the step); 0 for a splitting.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   real(R_P)                  :: c    !< Abscissa.

   c = 0._R_P
   if (self%form /= FORM_SPLITTING) c = self%bt%c(k)
   endfunction stage_time

   pure function stage_weight(self, k) result(w)
   !< Effective Butcher weight b_k of stage k; for a splitting, the fraction of the step of flow k.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   real(R_P)                  :: w    !< Weight.

   if (self%form == FORM_SPLITTING) then
      w = self%sc%coefficient(k)
   else
      w = self%bt%b(k)
   endif
   endfunction stage_weight

   pure function stage_partition(self, k) result(p)
   !< Partition advanced by flow k of a splitting (1 or 2); 0 for a Runge–Kutta form.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   integer(I4P)               :: p    !< Partition.

   p = 0_I4P
   if (self%form == FORM_SPLITTING) p = self%sc%partition(k)
   endfunction stage_partition

   pure function state_aliases_q(self) result(yes)
   !< The stage state is q itself (every form but the Butcher one).
   class(stepper), intent(in) :: self !< Stepper.
   logical                    :: yes  !< Aliasing.

   yes = self%form /= FORM_BUTCHER
   endfunction state_aliases_q

   pure function stage_state_is_q(self, k) result(yes)
   !< The state of stage k is q itself: the caller evaluates the residual of stage k on q. True for every stage of
   !< every form but the Butcher one; in the Butcher form, for the stages whose row of A is zero when the stepper was
   !< initialised with stages_on_q=.true.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   logical                    :: yes  !< The stage's state is q.

   yes = self%on_q(k)
   endfunction stage_state_is_q
endmodule selah_stepper
