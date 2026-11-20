!< SELAH, L2: additive splitting of a field step and a caller's flow.
module selah_additive_splitting
!< SELAH, L2: additive splitting of a field step and a caller's flow.
!<
!< For \(u' = F(u) + G(u)\), where F is the residual the stepper integrates (the field step, a Runge–Kutta step) and G
!< is a term the caller integrates itself (an exact or approximate sub-flow: radiative cooling with Townsend's exact
!< flow, a stiff relaxation, a particle step), one step is a composition of the two operators' flows over fractions of
!< the step. SELAH owns the composition (which operator, what fraction, at what time); the caller owns both flows.
!<
!<```fortran
!< call sp%initialize('strang')
!< do m=1, sp%substeps()
!<    h = sp%substep_fraction(m) * dt ; tm = t + sp%substep_time(m) * dt
!<    if (sp%substep_operator(m) == OPERATOR_FIELD) then
!<       ! one step of the stepper of length h, from time tm
!<    else
!<       ! the caller's flow over h, from time tm
!<    endif
!< enddo
!<```
!<
!<| Scheme | Substeps | Order | Composition (flow 1 first) |
!<|--------|----------|-------|----------------------------|
!<| `lie`    | 2 | 1 | field h, flow h (`first=OPERATOR_FLOW`: flow h, field h) |
!<| `strang` | 3 | 2 | flow h/2, field h, flow h/2 (`first=OPERATOR_FIELD`: field h/2, flow h, field h/2) |
!<| `yoshida-4`, `blanes-moan-s6`, `blanes-moan-s10` | 7, 13, 21 | 4, 4, 6 | the palindromes of selah_splittings |
!<
!< **Backward fractions.** Every real splitting of order 3 or more has a negative fraction in each operator (Sheng
!< 1989, Suzuki 1991, Goldman & Kaper 1996, Blanes & Casas 2005): a sub-flow run backward in time. That is fine for
!< reversible flows (a particle Hamiltonian) and ill-posed for an irreversible one: an upwind or limited field step,
!< cooling, relaxation, diffusion, accretion. `initialize` refuses such a scheme unless the caller declares both
!< operators reversible (`reversible=.true.`); `backward()` tells.
!<
!< **Order.** The order is the composition's when both sub-flows are exact, or each at least as accurate; the field
!< step of a scheme of order p contributes its own \(O(\Delta t^p)\). With a stiff flow the composition may lose order
!< (order reduction): see the verification of `selah_test_additive_splitting`.
!<
!< **Time.** Each operator integrates its own time dependence over its substep (the field step through its stages'
!< `stage_time`, the caller's flow from its start time), so each keeps its own clock: substep m starts at
!< `substep_time(m)`, the sum of the fractions of the **same** operator's substeps before m. This is the composition
!< of the autonomous system (u, t_1, t_2)' = (F(u, t_1), 1, 0) + (G(u, t_2), 0, 1), whose exact solution keeps
!< t_1 = t_2 = t: the composition keeps its order when either term or both depend on t. (The stepper's own splitting
!< form differs: its sub-flows are Euler steps of one partition, which cannot carry time; its `stage_time` takes the
!< other partition's fractions.)
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_splittings
implicit none
private
public :: additive_splitting
public :: OPERATOR_FIELD, OPERATOR_FLOW
public :: SCHEME_LIE
public :: ADDITIVE_SCHEMES

integer(I4P), parameter :: OPERATOR_FIELD = 1_I4P !< The field step: the stepper.
integer(I4P), parameter :: OPERATOR_FLOW  = 2_I4P !< The caller's flow.
character(*), parameter :: SCHEME_LIE = 'lie'     !< Lie (Lie–Trotter) splitting, order 1.
character(15), parameter :: ADDITIVE_SCHEMES(5) = [character(15) :: SCHEME_LIE, SCHEME_STRANG, SCHEME_YOSHIDA_4, &
                                                   SCHEME_BLANES_MOAN_S6, SCHEME_BLANES_MOAN_S10] !< Every scheme.

type :: additive_splitting
   !< A composition of the field step and the caller's flow.
   character(:), allocatable   :: scheme !< Scheme name.
   type(splitting_composition) :: sc     !< Composition: partition 1 = OPERATOR_FIELD, 2 = OPERATOR_FLOW.
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: substeps
      procedure, pass(self) :: substep_operator
      procedure, pass(self) :: substep_fraction
      procedure, pass(self) :: substep_time
      procedure, pass(self) :: order
      procedure, pass(self) :: backward
endtype additive_splitting

contains
   subroutine initialize(self, scheme, first, reversible)
   !< Set the composition; `first` is the operator of the first substep (default: the field for `lie`, ADAM's order,
   !< and the flow for the palindromes, so that the field step, the costly one, runs once in `strang`). A scheme with
   !< backward substeps stops the program unless `reversible` (see "Backward fractions" above).
   class(additive_splitting), intent(inout)        :: self       !< Splitting.
   character(*),              intent(in)           :: scheme     !< Scheme name (ADDITIVE_SCHEMES).
   integer(I4P),              intent(in), optional :: first      !< OPERATOR_FIELD or OPERATOR_FLOW.
   logical,                   intent(in), optional :: reversible !< Both operators may run backward (default .false.).
   integer(I4P)                                    :: f          !< First operator.
   logical                                         :: rev        !< Both operators reversible.

   self%scheme = trim(adjustl(scheme))
   if (self%scheme == SCHEME_LIE) then
      f = OPERATOR_FIELD ; if (present(first)) f = first
      self%sc%flows = 2_I4P
      self%sc%order = 1_I4P
      self%sc%partition = [f, 3_I4P - f]
      self%sc%coefficient = [1._R_P, 1._R_P]
   else
      f = OPERATOR_FLOW ; if (present(first)) f = first
      self%sc = splitting(self%scheme, first=f)
   endif
   if (f /= OPERATOR_FIELD .and. f /= OPERATOR_FLOW) error stop 'selah additive splitting: first must be 1 or 2'
   if (self%sc%flows == 0_I4P) error stop 'selah additive splitting: unknown scheme "'//self%scheme//'"'
   rev = .false. ; if (present(reversible)) rev = reversible
   if (self%backward() .and. .not. rev) &
      error stop 'selah additive splitting: "'//self%scheme//'" runs substeps backward in time: lie or strang for an '// &
                 'irreversible field step or flow (reversible=.true. if both operators are reversible)'
   endsubroutine initialize

   pure function substeps(self) result(n)
   !< Substeps of a step.
   class(additive_splitting), intent(in) :: self !< Splitting.
   integer(I4P)                          :: n    !< Substeps.

   n = self%sc%flows
   endfunction substeps

   pure function substep_operator(self, m) result(op)
   !< The operator of substep m: OPERATOR_FIELD (a step of the stepper) or OPERATOR_FLOW (the caller's flow).
   class(additive_splitting), intent(in) :: self !< Splitting.
   integer(I4P),              intent(in) :: m    !< Substep.
   integer(I4P)                          :: op   !< Operator.

   op = self%sc%partition(m)
   endfunction substep_operator

   pure function substep_fraction(self, m) result(c)
   !< The length of substep m as a fraction of the step (negative: a backward substep, see `backward`).
   class(additive_splitting), intent(in) :: self !< Splitting.
   integer(I4P),              intent(in) :: m    !< Substep.
   real(R_P)                             :: c    !< Fraction.

   c = self%sc%coefficient(m)
   endfunction substep_fraction

   pure function substep_time(self, m) result(c)
   !< The time at which substep m starts, as a fraction of the step: the fractions of the same operator's substeps
   !< before m (see "Time" above).
   class(additive_splitting), intent(in) :: self !< Splitting.
   integer(I4P),              intent(in) :: m    !< Substep.
   real(R_P)                             :: c    !< Fraction.

   c = sum(self%sc%coefficient(1:m-1), mask=self%sc%partition(1:m-1) == self%sc%partition(m))
   endfunction substep_time

   pure function order(self) result(p)
   !< The order of the composition with exact sub-flows.
   class(additive_splitting), intent(in) :: self !< Splitting.
   integer(I4P)                          :: p    !< Order.

   p = self%sc%order
   endfunction order

   pure function backward(self) result(yes)
   !< Some substep runs backward in time: admissible for reversible flows only.
   class(additive_splitting), intent(in) :: self !< Splitting.
   logical                               :: yes  !< A negative fraction.

   yes = any(self%sc%coefficient < 0._R_P)
   endfunction backward
endmodule selah_additive_splitting
