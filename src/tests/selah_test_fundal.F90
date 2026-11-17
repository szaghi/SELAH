!< SELAH test: the FUNDAL dependency is wired (header, modules, host mode, no MPI).

#include "fundal.H"

program selah_test_fundal
!< SELAH test: the FUNDAL dependency is wired (header, modules, host mode, no MPI).
!<
!< `fundal.H` is found on the include path, the FUNDAL modules compile with `-DFUNDAL_NO_MPI` (no MPI in SELAH's
!< modes), and the allocator works in host mode (no DEV_* macro): a rank-2 array with a lower bound of -1, initialised.
use fundal, only : dev_alloc, dev_free
use penf, only : I4P, R8P
implicit none
real(R8P), pointer :: a(:,:) => null()
integer(I4P)       :: ierr

call dev_alloc(fptr_dev=a, ubounds=[3_I4P, 4_I4P], lbounds=[1_I4P, -1_I4P], init_value=1._R8P, ierr=ierr)
if (ierr /= 0) error stop 'selah_test_fundal: dev_alloc failed'
if (lbound(a, 2) /= -1 .or. ubound(a, 2) /= 4) error stop 'selah_test_fundal: wrong bounds'
if (any(a /= 1._R8P)) error stop 'selah_test_fundal: wrong initial value'
call dev_free(a)
print '(A)', 'selah_test_fundal: ok'
endprogram selah_test_fundal
