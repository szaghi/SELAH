!< SELAH expected failure: a scheme asked for in a storage form it does not have stops the program.

program selah_test_xfail_form
!< SELAH expected failure: a scheme asked for in a storage form it does not have stops the program.
!<
!< rk4 has a Butcher form only (registers_needed(rk4, two-register) = -1). The runner (scripts/run_tests.sh) passes
!< this program when it exits with a non-zero status.
use selah
implicit none
type(stepper) :: st

call st%initialize(SCHEME_RK4, FORM_SHU_OSHER_2R, lb=[1, 0, 0, 0, 1], ub=[1, 5, 2, 2, 1]) ! must stop here
print '(A)', 'selah_test_xfail_form: the stepper accepted rk4 in the two-register form'
endprogram selah_test_xfail_form
