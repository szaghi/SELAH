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
!<    Both `rhs_register()` and `stage_state()` are null outside that window (before `begin_stage(k)`, after
!<    `end_stage(k)`): a caller that keeps them across stages fails at once instead of writing a register that is no
!<    longer the one of the active stage.
!< 5. After `finish_step` the registers are scratch.
!<
!< `stage_weight(k)` is the effective Butcher weight b_k of stage k in every Runge–Kutta form (the weight a flux
!< register needs to accumulate the step's boundary fluxes); for a splitting it is the fraction of the step of flow k.
!<
!< **Stage time.** `stage_time(k)` is the fraction of the step at which the residual of stage k is evaluated: the
!< abscissa c_k in every Runge–Kutta form. For a splitting, time is a variable with t' = 1, and a sub-flow is exact only
!< when its partition's residual does not depend on its own variables; so time advances with the flows of the partition
!< whose residual does not depend on t, and stays frozen through the other's. `stage_time(k)` is then the sum of the
!< fractions of the flows of the other partition before flow k: the time flow k needs when its own residual is the one
!< that depends on t, whichever partition that is (the composition keeps its order). When both residuals depend on t no
!< sub-flow is exact, and the composition gives no order guarantee for the time dependence.
!<
!< **Forward-Euler steps.** `ssp_coefficient()` is the SSP coefficient C of the scheme, from its Shu–Osher form (0 for a
!< scheme that is not SSP, and for a splitting). `stage_euler_fraction(k)` is the largest forward-Euler step, as a
!< fraction of dt, that the residual of stage k enters in the Shu–Osher form: f_k = max_i beta_ik / alpha_ik
!< (selah_tableaux `euler_step_fractions`; max_k f_k = 1/C). A limiter that makes the state of stage k plus
!< f_k dt times its residual admissible (positive, bounded) keeps every later stage and the step admissible. For a
!< splitting, flow k is a forward-Euler step of fraction c_k on its partition (0 for a backward flow).
!<
!< **Mask policy.** With `mask` present, the kernels update the active cells only. With `masked_q=.true.` at
!< initialisation, the masked cells of every stage state hold q^n: the Butcher form writes q^n into them as it builds
!< the stage state (in the same pass); the other forms keep the stage state in q, whose masked cells no update touches.
!< A caller whose residual reads masked cells (an immersed boundary extrapolating from the solid side) then sees q^n
!< there at every stage. Default `.false.`: the Butcher register's masked cells are left untouched (v0.1.0 behaviour).
!<
!< **Companion states.** A state that advances with the field in the same step, its residual built inside the field's
!< (particles pushed by the fields of each stage, a perfectly matched layer's auxiliary fields), is the same scheme
!< applied to another state: a companion is a second stepper, `initialize_companion(field, lb, ub, ...)`, with the
!< field stepper's scheme, form and stages on q, driven through the same calls in the same stage window (begin_stage
!< of both, the residuals of both, end_stage of both). Its bounds are a capacity: the box of each call is the active
!< part (particles 1..n of a pool, the components that evolve), and the kernels never touch the rest. Re-initialise it
!< when its shape changes (a regrid). Stepped together, field and companion are bitwise the two stepped as one state.
!<
!< **Stages on q.** The state of a stage whose row of A is zero (the first stage of every explicit scheme) is q itself.
!< With `stages_on_q=.true.` the Butcher form does not copy q into its stage-state register for such a stage: the
!< caller evaluates the residual on q, and `stage_state()` is null during that stage, so that a caller that does not
!< ask `stage_state_is_q(k)` fails at once instead of reading a stale register. Default `.false.`: every stage of the
!< Butcher form has its register (v0.1.0 behaviour).
!<
!< **Storage.** `registers_needed(scheme, form)` field-sized buffers, allocated through FUNDAL (device memory in a
!< device build) or passed by the caller (`registers`, a contiguous rank-6 array with the field's bounds and enough
!< slots; rank-5 fields only).
!<
!< **Fields of rank 1 to 4.** `initialize` takes the bounds of a field of rank k = size(lb) and, below rank 5, the
!< `components` flag (selah_box: the field viewed as rank 5 with singleton dimensions). `begin_step`, `begin_stage`,
!< `end_stage` and `finish_step` are generic over the rank of q (dq and mask have q's rank); the box must have the
!< field's rank (`box_interior(lb, ub, ngc, components)`). The field reaches the one rank-5 implementation through an
!< explicit-shape dummy, with no copy when it is contiguous; it must be contiguous in a device build. The stage state
!< and the residual register of a rank-k field are `stage_state_view(u)` and `rhs_register_view(r)`: rank-k pointers
!< with the field's bounds, null when `stage_state()` and `rhs_register()` are (the functions keep their rank-5
!< result).
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
   integer(I4P)                :: lb(5) = 1_I4P           !< Lower bounds of the field, as rank 5.
   integer(I4P)                :: ub(5) = 1_I4P           !< Upper bounds of the field, as rank 5.
   integer(I4P)                :: rank = 5_I4P            !< Rank of the field.
   logical                     :: components = .true.     !< The field's first dimension holds its components.
   integer(I4P)                :: lbk(5) = 1_I4P          !< Lower bounds of the field, rank k (the first rank values).
   integer(I4P)                :: ubk(5) = 1_I4P          !< Upper bounds of the field, rank k (the first rank values).
   integer(I4P)                :: part(2,2) = 0_I4P       !< Component ranges of the two partitions (splitting).
   type(butcher_tableau)       :: bt                      !< Butcher tableau (every RK form).
   type(shu_osher_2r_form)     :: r2                      !< Two-register form.
   type(williamson_2n_form)    :: w2n                     !< 2N form.
   type(splitting_composition) :: sc                      !< Splitting.
   real(R_P), pointer, contiguous :: regs(:,:,:,:,:,:) => null() !< Registers (contiguous: the rank-k views remap them).
   real(R_P), pointer          :: at_d(:,:) => null()     !< Transposed A (column k = row k), device resident.
   real(R_P), pointer          :: b_d(:) => null()        !< Weights b, device resident.
   logical                     :: owned = .false.         !< The registers are the stepper's.
   logical, allocatable        :: on_q(:)                 !< Stage k's state is q itself.
   real(R_P), allocatable      :: ef(:)                   !< Forward-Euler fraction of each stage (Shu–Osher form).
   real(R_P)                   :: cssp = 0._R_P           !< SSP coefficient of the Shu–Osher form, 0 if none.
   logical                     :: masked_q = .false.      !< The masked cells of every stage state are q^n.
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: initialize_companion
      procedure, pass(self) :: destroy
      procedure, pass(self), private :: begin_step_r1, begin_step_r2, begin_step_r3, begin_step_r4, begin_step_r5
      procedure, pass(self), private :: begin_stage_r1, begin_stage_r2, begin_stage_r3, begin_stage_r4, begin_stage_r5
      procedure, pass(self), private :: end_stage_r1, end_stage_r2, end_stage_r3, end_stage_r4, end_stage_r5
      procedure, pass(self), private :: finish_step_r1, finish_step_r2, finish_step_r3, finish_step_r4, finish_step_r5
      procedure, pass(self), private :: stage_state_view_r1, stage_state_view_r2, stage_state_view_r3, &
                                        stage_state_view_r4, stage_state_view_r5
      procedure, pass(self), private :: rhs_register_view_r1, rhs_register_view_r2, rhs_register_view_r3, &
                                        rhs_register_view_r4, rhs_register_view_r5
      generic :: begin_step => begin_step_r1, begin_step_r2, begin_step_r3, begin_step_r4, begin_step_r5
      generic :: begin_stage => begin_stage_r1, begin_stage_r2, begin_stage_r3, begin_stage_r4, begin_stage_r5
      generic :: end_stage => end_stage_r1, end_stage_r2, end_stage_r3, end_stage_r4, end_stage_r5
      generic :: finish_step => finish_step_r1, finish_step_r2, finish_step_r3, finish_step_r4, finish_step_r5
      generic :: stage_state_view => stage_state_view_r1, stage_state_view_r2, stage_state_view_r3, &
                                     stage_state_view_r4, stage_state_view_r5
      generic :: rhs_register_view => rhs_register_view_r1, rhs_register_view_r2, rhs_register_view_r3, &
                                      rhs_register_view_r4, rhs_register_view_r5
      procedure, pass(self) :: stage_state
      procedure, pass(self) :: rhs_register
      procedure, pass(self) :: stages
      procedure, pass(self) :: stage_time
      procedure, pass(self) :: stage_weight
      procedure, pass(self) :: stage_partition
      procedure, pass(self) :: stage_euler_fraction
      procedure, pass(self) :: ssp_coefficient => stepper_ssp_coefficient
      procedure, pass(self) :: imag_axis_bound => stepper_imag_axis_bound
      procedure, pass(self) :: real_axis_bound => stepper_real_axis_bound
      procedure, pass(self) :: state_aliases_q
      procedure, pass(self) :: stage_state_is_q
endtype stepper

contains
   subroutine initialize(self, scheme, form, lb, ub, partitions, registers, stages_on_q, components, masked_q)
   !< Set the scheme and its storage form for a field of rank size(lb); allocate the registers (or take the caller's).
   class(stepper),     intent(inout)                  :: self          !< Stepper.
   character(*),       intent(in)                     :: scheme        !< Scheme name (selah_schemes, selah_splittings).
   character(*),       intent(in)                     :: form          !< Storage form (selah_storage).
   integer(I4P),       intent(in)                     :: lb(:), ub(:)  !< Bounds of the field, ghosts included (rank k).
   integer(I4P),       intent(in), optional           :: partitions(2,2) !< Splitting: (v1, v2) of partitions 1, 2.
   real(R_P),          intent(in), optional, pointer  :: registers(:,:,:,:,:,:) !< The caller's registers.
   logical,            intent(in), optional           :: stages_on_q   !< Butcher form: the stages whose row of A is
                                                                       !< zero are evaluated on q (default .false.).
   logical,            intent(in), optional           :: components    !< Rank below 5: the field's first dimension
                                                                       !< holds its components (selah_box).
   logical,            intent(in), optional           :: masked_q      !< The masked cells of every stage state are
                                                                       !< q^n (default .false.: left untouched).
   integer(I4P)                                       :: nreg, ierr  !< Registers needed, error status.
   integer(I4P)                                       :: k           !< Stage.
   type(box)                                          :: bx          !< The field's dimensions as rank 5.
   real(R_P), allocatable                             :: at(:,:)     !< Transposed A, a contiguous array.

   call self%destroy
   self%scheme = trim(adjustl(scheme))
   self%form = trim(adjustl(form))
   self%rank = size(lb, kind=I4P)
   if (size(ub) /= self%rank) error stop 'selah stepper: lb and ub must have the rank of the field'
   if (self%rank < 5_I4P .and. .not. present(components)) &
      error stop 'selah stepper: a field of rank below 5 needs components= (its first dimension holds components?)'
   self%components = .true.
   if (present(components)) self%components = components
   self%masked_q = .false.
   if (present(masked_q)) self%masked_q = masked_q
   if (self%rank == 5_I4P .and. .not. self%components) error stop 'selah stepper: a rank-5 field has components'
   bx = box_range(lb, lb, ub, self%components)
   self%lb = bx%lb
   self%ub = bx%hi
   self%lbk = 1_I4P ; self%ubk = 1_I4P
   self%lbk(1:self%rank) = lb
   self%ubk(1:self%rank) = ub
   nreg = registers_needed(self%scheme, self%form)
   if (nreg < 0) error stop 'selah stepper: scheme "'//self%scheme//'" has no form "'//self%form//'"'
   self%bt = butcher(self%scheme)
   select case(self%form)
   case(FORM_BUTCHER)
      self%nstages = self%bt%stages
      ! A^T through a contiguous array, not the expression transpose(A): gfortran 14.3 passes that expression on to
      ! FUNDAL's contiguous dummy without the copy it needs, so at_d held A itself (every Butcher stage state but the
      ! first wrong; seen in ADAM's build, -cpp -O0/-O2; SELAH's -std=f2018 build hides it)
      at = transpose(self%bt%a)
      call dev_assign_to_device(dst=self%at_d, src=at)
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
   allocate(self%on_q(self%nstages), self%ef(self%nstages))
   call euler_steps(self)
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
      if (self%rank /= 5_I4P) error stop 'selah stepper: the caller passes registers for rank-5 fields only'
      if (.not. associated(registers)) error stop 'selah stepper: the registers passed are not associated'
      if (.not. is_contiguous(registers)) error stop 'selah stepper: the registers passed must be contiguous'
      if (size(registers, 6) < nreg) error stop 'selah stepper: too few registers passed'
      if (lbound(registers, 6) /= 1) error stop 'selah stepper: the registers passed must be numbered from 1 (dim 6)'
      if (any(lbound(registers) /= [self%lb, 1_I4P]) .or. &
          any(ubound(registers, kind=I4P) /= [self%ub, ubound(registers, 6, kind=I4P)])) &
         error stop 'selah stepper: the registers passed do not have the bounds of the field (lb, ub)'
      self%regs => registers
      self%owned = .false.
   else
      ! zeroed, so that every register starts finite (the 2N form's first stage no longer reads Delta Q: A_1 = 0)
      call dev_alloc(fptr_dev=self%regs, ubounds=[self%ub, nreg], lbounds=[self%lb, 1_I4P], init_value=0._R_P, ierr=ierr)
      if (ierr /= 0) error stop 'selah stepper: the registers could not be allocated'
      self%owned = .true.
   endif
   endsubroutine initialize

   subroutine initialize_companion(self, field, lb, ub, partitions, registers, components)
   !< A companion of the stepper `field`: a stepper of the same scheme, form and stages on q, for another state that
   !< advances with the field in the same step (see "Companion states" above).
   class(stepper),     intent(inout)                  :: self          !< Companion.
   class(stepper),     intent(in)                     :: field         !< The field's stepper, initialised.
   integer(I4P),       intent(in)                     :: lb(:), ub(:)  !< Bounds of the companion's state (rank k).
   integer(I4P),       intent(in), optional           :: partitions(2,2) !< Splitting: the companion's partitions.
   real(R_P),          intent(in), optional, pointer  :: registers(:,:,:,:,:,:) !< The caller's registers (rank 5).
   logical,            intent(in), optional           :: components    !< Rank below 5: the first dimension holds
                                                                       !< components.

   if (field%nstages == 0_I4P) error stop 'selah stepper: a companion needs an initialised field stepper'
   if (field%form == FORM_SPLITTING .and. .not. present(partitions)) &
      error stop 'selah stepper: a companion of a splitting needs the components of its own partitions'
   call self%initialize(field%scheme, field%form, lb, ub, partitions=partitions, registers=registers, &
                        stages_on_q=any(field%on_q) .and. field%form == FORM_BUTCHER, components=components)
   endsubroutine initialize_companion

   subroutine euler_steps(self)
   !< The forward-Euler fraction of each stage and the SSP coefficient: from the scheme's Shu–Osher form, whatever the
   !< storage form (the stages are the same, in the same order); a splitting's flow k is an Euler step of fraction c_k
   !< on its partition (0 for a backward flow, c_k < 0).
   class(stepper), intent(inout) :: self !< Stepper.
   type(shu_osher_form)          :: so   !< Shu–Osher form of the scheme.

   self%ef = 0._R_P
   self%cssp = 0._R_P
   if (self%form == FORM_SPLITTING) then
      self%ef = max(self%sc%coefficient, 0._R_P)
   else
      so = shu_osher(self%scheme)
      if (so%stages == self%nstages) then
         self%ef = euler_step_fractions(so)
         self%cssp = ssp_coefficient(so)
      endif
   endif
   endsubroutine euler_steps

   subroutine destroy(self)
   !< Free what the stepper owns.
   class(stepper), intent(inout) :: self !< Stepper.

   if (self%owned .and. associated(self%regs)) call dev_free(self%regs)
   if (associated(self%at_d)) call dev_free(self%at_d)
   if (associated(self%b_d)) call dev_free(self%b_d)
   nullify(self%regs, self%at_d, self%b_d)
   if (allocated(self%on_q)) deallocate(self%on_q)
   if (allocated(self%ef)) deallocate(self%ef)
   self%cssp = 0._R_P
   self%owned = .false. ; self%nstages = 0_I4P ; self%k_active = 0_I4P
   self%rank = 5_I4P ; self%components = .true. ; self%masked_q = .false.
   endsubroutine destroy

   subroutine begin_step_r5(self, bx, q, mask, mask_comp)
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
   endsubroutine begin_step_r5

   subroutine begin_stage_r5(self, k, bx, dt, q, mask, mask_comp)
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
                             mask, mask_comp, masked_q=self%masked_q)
   self%k_active = k
   endsubroutine begin_stage_r5

   subroutine end_stage_r5(self, k, bx, dt, q, mask, mask_comp)
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
   endsubroutine end_stage_r5

   subroutine finish_step_r5(self, bx, dt, q, dq, mask, mask_comp)
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
   endsubroutine finish_step_r5

   function stage_state(self) result(u)
   !< The state of the active stage, with the field's bounds: the Butcher form's register; null when the stage's state
   !< is q itself (stage_state_is_q), so in every other form, and outside a stage window.
   class(stepper), intent(in) :: self         !< Stepper.
   real(R_P),      pointer    :: u(:,:,:,:,:) !< Stage state.

   integer(I4P)               :: slot         !< Register slot.

   u => null()
   slot = state_slot(self)
   if (slot > 0_I4P) u(self%lb(1):,self%lb(2):,self%lb(3):,self%lb(4):,self%lb(5):) => self%regs(:,:,:,:,:,slot)
   endfunction stage_state

   function rhs_register(self) result(r)
   !< Where the caller writes the residual of the active stage, with the field's bounds; null outside a stage window
   !< (in the Butcher form, slot k_active + 1 would otherwise be the stage-state register).
   class(stepper), intent(in) :: self         !< Stepper.
   real(R_P),      pointer    :: r(:,:,:,:,:) !< Residual register.
   integer(I4P)               :: slot         !< Register slot.

   r => null()
   slot = rhs_slot(self)
   if (slot > 0_I4P) r(self%lb(1):,self%lb(2):,self%lb(3):,self%lb(4):,self%lb(5):) => self%regs(:,:,:,:,:,slot)
   endfunction rhs_register

   pure function state_slot(self) result(slot)
   !< The register slot of the stage state, 0 when it is null (no active stage, a stage on q, a form other than Butcher).
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P)               :: slot !< Register slot.

   slot = 0_I4P
   if (self%k_active == 0_I4P) return
   if (self%on_q(self%k_active)) return
   if (self%form == FORM_BUTCHER) slot = 1_I4P
   endfunction state_slot

   pure function rhs_slot(self) result(slot)
   !< The register slot of the active stage's residual, 0 outside a stage window.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P)               :: slot !< Register slot.

   slot = 0_I4P
   if (self%k_active == 0_I4P) return
   select case(self%form)
   case(FORM_BUTCHER)
      slot = self%k_active + 1_I4P
   case(FORM_SPLITTING)
      slot = 1_I4P
   case default
      slot = 2_I4P
   endselect
   endfunction rhs_slot

   subroutine check_front(self, bx, rank, shp, ub)
   !< A rank-k front: the box and the stepper have the field's rank, q has the field's shape; ub, q's upper bounds as
   !< rank 5.
   class(stepper), intent(in)  :: self     !< Stepper.
   type(box),      intent(in)  :: bx       !< Box.
   integer(I4P),   intent(in)  :: rank     !< Rank of the front.
   integer(I4P),   intent(in)  :: shp(:)   !< Shape of q.
   integer(I4P),   intent(out) :: ub(5)    !< Upper bounds of q as rank 5.

   if (self%rank /= rank) error stop 'selah stepper: q does not have the rank the stepper was initialised with'
   if (bx%rank /= rank .or. (bx%components .neqv. self%components)) &
      error stop 'selah stepper: the box does not have the rank (or the components) of q'
   ub = front_ub(bx, shp)
   if (any(bx%lb /= self%lb) .or. any(ub /= self%ub)) &
      error stop 'selah stepper: q does not have the bounds the stepper was initialised with'
   endsubroutine check_front

   pure function front_ub(bx, shp) result(ub)
   !< Upper bounds, as rank 5, of an array of shape shp with the box's lower bounds (a mask: its own components).
   type(box),    intent(in) :: bx     !< Box.
   integer(I4P), intent(in) :: shp(:) !< Shape.
   integer(I4P)             :: ub(5)  !< Upper bounds.
   integer(I4P)             :: e(5), d !< Extents as rank 5, dimension.

   e = 1_I4P
   do d=1_I4P, size(shp, kind=I4P)
      e(box_dim(bx, d)) = shp(d)
   enddo
   ub = bx%lb + e - 1_I4P
   endfunction front_ub

   ! ---- the explicit-shape bridges: a rank-k field reaches the rank-5 procedures as the rank-5 array it is in memory
   subroutine begin_step_x(self, bx, ubq, q, ubm, mask, mask_comp)
   !< begin_step through an explicit-shape rank-5 view.
   class(stepper), intent(inout)        :: self    !< Stepper.
   type(box),      intent(in)           :: bx      !< Box.
   integer(I4P),   intent(in)           :: ubq(5)  !< Upper bounds of q.
   real(R_P),      intent(in)           :: q(bx%lb(1):ubq(1),bx%lb(2):ubq(2),bx%lb(3):ubq(3),bx%lb(4):ubq(4), &
                                             bx%lb(5):ubq(5)) !< Solution.
   integer(I4P),   intent(in)           :: ubm(5)  !< Upper bounds of the mask.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):ubm(1),bx%lb(2):ubm(2),bx%lb(3):ubm(3),bx%lb(4):ubm(4), &
                                                bx%lb(5):ubm(5)) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp !< Mask component.

   call self%begin_step_r5(bx, q, mask, mask_comp)
   endsubroutine begin_step_x

   subroutine begin_stage_x(self, k, bx, dt, ubq, q, ubm, mask, mask_comp)
   !< begin_stage through an explicit-shape rank-5 view.
   class(stepper), intent(inout)        :: self    !< Stepper.
   integer(I4P),   intent(in)           :: k       !< Stage.
   type(box),      intent(in)           :: bx      !< Box.
   real(R_P),      intent(in)           :: dt      !< Time step.
   integer(I4P),   intent(in)           :: ubq(5)  !< Upper bounds of q.
   real(R_P),      intent(in)           :: q(bx%lb(1):ubq(1),bx%lb(2):ubq(2),bx%lb(3):ubq(3),bx%lb(4):ubq(4), &
                                             bx%lb(5):ubq(5)) !< Solution.
   integer(I4P),   intent(in)           :: ubm(5)  !< Upper bounds of the mask.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):ubm(1),bx%lb(2):ubm(2),bx%lb(3):ubm(3),bx%lb(4):ubm(4), &
                                                bx%lb(5):ubm(5)) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp !< Mask component.

   call self%begin_stage_r5(k, bx, dt, q, mask, mask_comp)
   endsubroutine begin_stage_x

   subroutine end_stage_x(self, k, bx, dt, ubq, q, ubm, mask, mask_comp)
   !< end_stage through an explicit-shape rank-5 view.
   class(stepper), intent(inout)        :: self    !< Stepper.
   integer(I4P),   intent(in)           :: k       !< Stage.
   type(box),      intent(in)           :: bx      !< Box.
   real(R_P),      intent(in)           :: dt      !< Time step.
   integer(I4P),   intent(in)           :: ubq(5)  !< Upper bounds of q.
   real(R_P),      intent(inout)        :: q(bx%lb(1):ubq(1),bx%lb(2):ubq(2),bx%lb(3):ubq(3),bx%lb(4):ubq(4), &
                                             bx%lb(5):ubq(5)) !< Solution.
   integer(I4P),   intent(in)           :: ubm(5)  !< Upper bounds of the mask.
   real(R_P),      intent(in), optional :: mask(bx%lb(1):ubm(1),bx%lb(2):ubm(2),bx%lb(3):ubm(3),bx%lb(4):ubm(4), &
                                                bx%lb(5):ubm(5)) !< Mask.
   integer(I4P),   intent(in), optional :: mask_comp !< Mask component.

   call self%end_stage_r5(k, bx, dt, q, mask, mask_comp)
   endsubroutine end_stage_x

   subroutine finish_step_x(self, bx, dt, ubq, q, dq, ubm, mask, mask_comp)
   !< finish_step through an explicit-shape rank-5 view.
   class(stepper), intent(inout)           :: self   !< Stepper.
   type(box),      intent(in)              :: bx     !< Box.
   real(R_P),      intent(in)              :: dt     !< Time step.
   integer(I4P),   intent(in)              :: ubq(5) !< Upper bounds of q (and dq).
   real(R_P),      intent(inout)           :: q(bx%lb(1):ubq(1),bx%lb(2):ubq(2),bx%lb(3):ubq(3),bx%lb(4):ubq(4), &
                                                bx%lb(5):ubq(5)) !< Solution.
   real(R_P),      intent(inout), optional :: dq(bx%lb(1):ubq(1),bx%lb(2):ubq(2),bx%lb(3):ubq(3),bx%lb(4):ubq(4), &
                                                 bx%lb(5):ubq(5)) !< Step residual.
   integer(I4P),   intent(in)              :: ubm(5) !< Upper bounds of the mask.
   real(R_P),      intent(in),    optional :: mask(bx%lb(1):ubm(1),bx%lb(2):ubm(2),bx%lb(3):ubm(3),bx%lb(4):ubm(4), &
                                                   bx%lb(5):ubm(5)) !< Mask.
   integer(I4P),   intent(in),    optional :: mask_comp !< Mask component.

   call self%finish_step_r5(bx, dt, q, dq, mask, mask_comp)
   endsubroutine finish_step_x

   ! ---- the rank-k fronts (selah_stepper_front.INC, once per rank)
#define SELAH_FRONT_RANK 1_I4P
#define SELAH_FRONT_DIMS :
#define SELAH_BEGIN_STEP begin_step_r1
#define SELAH_BEGIN_STAGE begin_stage_r1
#define SELAH_END_STAGE end_stage_r1
#define SELAH_FINISH_STEP finish_step_r1
#define SELAH_STATE_VIEW stage_state_view_r1
#define SELAH_RHS_VIEW rhs_register_view_r1
#define SELAH_VIEW_BOUNDS l(1):h(1)
#include "selah_stepper_front.INC"
#undef SELAH_FRONT_RANK
#undef SELAH_FRONT_DIMS
#undef SELAH_BEGIN_STEP
#undef SELAH_BEGIN_STAGE
#undef SELAH_END_STAGE
#undef SELAH_FINISH_STEP
#undef SELAH_STATE_VIEW
#undef SELAH_RHS_VIEW
#undef SELAH_VIEW_BOUNDS
#define SELAH_FRONT_RANK 2_I4P
#define SELAH_FRONT_DIMS :,:
#define SELAH_BEGIN_STEP begin_step_r2
#define SELAH_BEGIN_STAGE begin_stage_r2
#define SELAH_END_STAGE end_stage_r2
#define SELAH_FINISH_STEP finish_step_r2
#define SELAH_STATE_VIEW stage_state_view_r2
#define SELAH_RHS_VIEW rhs_register_view_r2
#define SELAH_VIEW_BOUNDS l(1):h(1),l(2):h(2)
#include "selah_stepper_front.INC"
#undef SELAH_FRONT_RANK
#undef SELAH_FRONT_DIMS
#undef SELAH_BEGIN_STEP
#undef SELAH_BEGIN_STAGE
#undef SELAH_END_STAGE
#undef SELAH_FINISH_STEP
#undef SELAH_STATE_VIEW
#undef SELAH_RHS_VIEW
#undef SELAH_VIEW_BOUNDS
#define SELAH_FRONT_RANK 3_I4P
#define SELAH_FRONT_DIMS :,:,:
#define SELAH_BEGIN_STEP begin_step_r3
#define SELAH_BEGIN_STAGE begin_stage_r3
#define SELAH_END_STAGE end_stage_r3
#define SELAH_FINISH_STEP finish_step_r3
#define SELAH_STATE_VIEW stage_state_view_r3
#define SELAH_RHS_VIEW rhs_register_view_r3
#define SELAH_VIEW_BOUNDS l(1):h(1),l(2):h(2),l(3):h(3)
#include "selah_stepper_front.INC"
#undef SELAH_FRONT_RANK
#undef SELAH_FRONT_DIMS
#undef SELAH_BEGIN_STEP
#undef SELAH_BEGIN_STAGE
#undef SELAH_END_STAGE
#undef SELAH_FINISH_STEP
#undef SELAH_STATE_VIEW
#undef SELAH_RHS_VIEW
#undef SELAH_VIEW_BOUNDS
#define SELAH_FRONT_RANK 4_I4P
#define SELAH_FRONT_DIMS :,:,:,:
#define SELAH_BEGIN_STEP begin_step_r4
#define SELAH_BEGIN_STAGE begin_stage_r4
#define SELAH_END_STAGE end_stage_r4
#define SELAH_FINISH_STEP finish_step_r4
#define SELAH_STATE_VIEW stage_state_view_r4
#define SELAH_RHS_VIEW rhs_register_view_r4
#define SELAH_VIEW_BOUNDS l(1):h(1),l(2):h(2),l(3):h(3),l(4):h(4)
#include "selah_stepper_front.INC"
#undef SELAH_FRONT_RANK
#undef SELAH_FRONT_DIMS
#undef SELAH_BEGIN_STEP
#undef SELAH_BEGIN_STAGE
#undef SELAH_END_STAGE
#undef SELAH_FINISH_STEP
#undef SELAH_STATE_VIEW
#undef SELAH_RHS_VIEW
#undef SELAH_VIEW_BOUNDS

   subroutine stage_state_view_r5(self, u)
   !< The stage state of a rank-5 field (as stage_state()).
   class(stepper), intent(in)          :: self         !< Stepper.
   real(R_P), pointer, intent(out)     :: u(:,:,:,:,:) !< Stage state, null when stage_state() is.

   u => self%stage_state()
   endsubroutine stage_state_view_r5

   subroutine rhs_register_view_r5(self, r)
   !< The residual register of a rank-5 field (as rhs_register()).
   class(stepper), intent(in)          :: self         !< Stepper.
   real(R_P), pointer, intent(out)     :: r(:,:,:,:,:) !< Residual register, null when rhs_register() is.

   r => self%rhs_register()
   endsubroutine rhs_register_view_r5

   pure function stages(self) result(n)
   !< Stages of a step (flows of a splitting).
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P)               :: n    !< Stages.

   n = self%nstages
   endfunction stages

   pure function stage_time(self, k) result(c)
   !< Fraction of the step at which the residual of stage k is evaluated: the abscissa c_k; for a splitting, the
   !< fractions of the other partition's flows before flow k (see "Stage time" above).
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   real(R_P)                  :: c    !< Abscissa.

   if (self%form == FORM_SPLITTING) then
      c = sum(self%sc%coefficient(1:k-1), mask=self%sc%partition(1:k-1) /= self%sc%partition(k))
   else
      c = self%bt%c(k)
   endif
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

   pure function stage_euler_fraction(self, k) result(f)
   !< The forward-Euler step stage k's residual takes, as a fraction of the step (see "Forward-Euler steps" above): the
   !< step a limiter makes admissible on the state of stage k; 0 when the scheme has no such decomposition.
   class(stepper), intent(in) :: self !< Stepper.
   integer(I4P),   intent(in) :: k    !< Stage.
   real(R_P)                  :: f    !< Fraction of the step.

   f = self%ef(k)
   endfunction stage_euler_fraction

   pure function stepper_ssp_coefficient(self) result(c)
   !< The SSP coefficient C of the scheme (its Shu–Osher form): the scheme keeps every convex property of forward Euler
   !< for dt <= C dt_FE; 0 when the scheme is not SSP, and for a splitting.
   class(stepper), intent(in) :: self !< Stepper.
   real(R_P)                  :: c    !< SSP coefficient.

   c = self%cssp
   endfunction stepper_ssp_coefficient

   pure function stepper_imag_axis_bound(self) result(beta)
   !< The scheme's imaginary-axis stability bound (selah_tableaux `imag_axis_bound`): an oscillatory mode of frequency
   !< omega (a gyration, a wave) is stable for omega dt <= beta; 0 for euler, ssp-22 and a splitting.
   class(stepper), intent(in) :: self !< Stepper.
   real(R_P)                  :: beta !< Bound.

   beta = 0._R_P
   if (self%form /= FORM_SPLITTING) beta = imag_axis_bound(self%bt)
   endfunction stepper_imag_axis_bound

   pure function stepper_real_axis_bound(self) result(beta)
   !< The scheme's real-axis stability bound (selah_tableaux `real_axis_bound`): a decay of rate lambda is stable for
   !< lambda dt <= beta; 0 for a splitting.
   class(stepper), intent(in) :: self !< Stepper.
   real(R_P)                  :: beta !< Bound.

   beta = 0._R_P
   if (self%form /= FORM_SPLITTING) beta = real_axis_bound(self%bt)
   endfunction stepper_real_axis_bound

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
