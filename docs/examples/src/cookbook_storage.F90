!run storage-pool cookbook_storage pool
!run storage-precision cookbook_storage precision
program cookbook_storage
!< Cookbook, storage: registers from the caller's pool; the working precision.
use selah
implicit none
character(16) :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('pool')      ; call pool
case('precision') ; call precision_of_the_build
case default      ; error stop 'usage: cookbook_storage pool|precision'
endselect

contains
   subroutine pool
   !< One register pool, sized once for the largest need, serving steppers of different schemes in turn.
   !region pool
   integer(I4P), parameter :: LB(5) = [1, -2, -2, -2, 1], UB(5) = [5, 9, 9, 9, 8]   ! the field's bounds
   character(7), parameter :: SCHEMES(3) = [character(7) :: SCHEME_SSP_54, SCHEME_RK4, SCHEME_SSP_33]
   real(R_P), pointer      :: regs(:,:,:,:,:,:) => null()
   type(stepper)           :: st
   integer(I4P)            :: s, slots

   slots = 0
   do s=1, size(SCHEMES)
      slots = max(slots, registers_needed(trim(SCHEMES(s)), FORM_BUTCHER))
   enddo
   allocate(regs(LB(1):UB(1), LB(2):UB(2), LB(3):UB(3), LB(4):UB(4), LB(5):UB(5), 1:slots))  ! slots from 1
   print '(A,I0,A,F6.1,A)', 'pool: ', slots, ' registers, ', real(size(regs), R_P) * storage_size(regs) / 8 / 2**20, ' MiB'
   do s=1, size(SCHEMES)
      call st%initialize(trim(SCHEMES(s)), FORM_BUTCHER, lb=LB, ub=UB, registers=regs)   ! the stepper uses, never frees
      print '(A,A,I0,A)', trim(SCHEMES(s)), ': uses ', registers_needed(trim(SCHEMES(s)), FORM_BUTCHER), ' of the pool'
      call st%destroy
   enddo
   deallocate(regs)
   !endregion pool
   endsubroutine pool

   subroutine precision_of_the_build
   !< The kind SELAH computes in: R8P by default, R4P in the *-r4 build modes (-D_R_P_IS_R4P for SELAH and PENF).
   !region precision
   print '(A,I0,A,I0)', 'R_P: kind ', R_P, ', decimal digits ', precision(1._R_P)
   print '(A,L1)', 'double precision build: ', R_P == R8P
   !endregion precision
   endsubroutine precision_of_the_build
endprogram cookbook_storage
