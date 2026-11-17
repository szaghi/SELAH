!< SELAH expected failure: caller registers without the field's bounds stop the program.

#include "fundal.H"

program selah_test_xfail_registers
!< SELAH expected failure: caller registers without the field's bounds stop the program.
!<
!< `initialize(registers=)` takes a rank-6 array with the bounds of the field and the slots numbered from 1. Here the
!< slots are numbered from 0, which the stepper would otherwise index one slot off. The runner
!< (scripts/run_tests.sh) passes this program when it exits with a non-zero status.
use fundal, only : dev_alloc
use selah
implicit none
type(stepper)      :: st
real(R_P), pointer :: regs_d(:,:,:,:,:,:)=>null()
integer(I4P)       :: ierr

call dev_alloc(fptr_dev=regs_d, ubounds=[1, 5, 2, 2, 1, 5], lbounds=[1, 0, 0, 0, 1, 0], init_value=0._R_P, ierr=ierr)
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=[1, 0, 0, 0, 1], ub=[1, 5, 2, 2, 1], registers=regs_d) ! stops
print '(A)', 'selah_test_xfail_registers: the stepper accepted registers numbered from 0'
endprogram selah_test_xfail_registers
