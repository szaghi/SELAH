!< SELAH expected failure: the step residual dq asked of a form other than the Butcher one stops the program.

#include "fundal.H"

program selah_test_xfail_dq_form
!< SELAH expected failure: the step residual dq asked of a form other than the Butcher one stops the program.
!<
!< `finish_step(dq=)` is the Butcher form's only (the in-place forms never hold every stage residual). The runner
!< (scripts/run_tests.sh) passes this program when it exits with a non-zero status.
use fundal, only : dev_alloc
use selah
implicit none
type(stepper)      :: st
type(box)          :: bx
real(R_P), pointer :: q_d(:,:,:,:,:)=>null(), dq_d(:,:,:,:,:)=>null()
integer(I4P)       :: ierr

bx = box_interior(LAYOUT_V_FIRST, 1_I4P, 4_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P)
call dev_alloc(fptr_dev=q_d, ubounds=[1, 5, 2, 2, 1], lbounds=[1, 0, 0, 0, 1], init_value=0._R_P, ierr=ierr)
call dev_alloc(fptr_dev=dq_d, ubounds=[1, 5, 2, 2, 1], lbounds=[1, 0, 0, 0, 1], init_value=0._R_P, ierr=ierr)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=lbound(q_d), ub=ubound(q_d))
call st%finish_step(bx, 0.1_R_P, q_d, dq=dq_d) ! must stop here
print '(A)', 'selah_test_xfail_dq_form: the stepper accepted dq for the two-register form'
endprogram selah_test_xfail_dq_form
