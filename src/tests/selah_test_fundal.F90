!< SELAH test: the FUNDAL dependency is wired (header, modules, no MPI), on the host and on the device.

#include "fundal.H"

program selah_test_fundal
!< SELAH test: the FUNDAL dependency is wired (header, modules, no MPI), on the host and on the device.
!<
!< `fundal.H` is found on the include path, the FUNDAL modules compile with `-DFUNDAL_NO_MPI` (no MPI in SELAH's
!< modes), and the allocator works: a rank-2 array with a lower bound of -1, initialised. In a device build the array
!< is device memory: its bounds are checked on the descriptor, its values after a copy back, and
!< `dev_init(require_device=.true.)` forbids a silent host fallback.
use fundal, only : dev_alloc, dev_free, dev_memcpy_from_device
#if defined DEV_OAC || defined DEV_OMP
use fundal, only : dev_init
#endif
use penf, only : I4P, R8P
implicit none
real(R8P), pointer :: a(:,:) => null()
real(R8P)          :: h(1:3,-1:4)
integer(I4P)       :: ierr

#if defined DEV_OAC || defined DEV_OMP
call dev_init(require_device=.true.)
#endif
call dev_alloc(fptr_dev=a, ubounds=[3_I4P, 4_I4P], lbounds=[1_I4P, -1_I4P], init_value=1._R8P, ierr=ierr)
if (ierr /= 0) error stop 'selah_test_fundal: dev_alloc failed'
if (lbound(a, 2) /= -1 .or. ubound(a, 2) /= 4) error stop 'selah_test_fundal: wrong bounds'
h = 0._R8P
call dev_memcpy_from_device(dst=h, src=a)
if (any(h /= 1._R8P)) error stop 'selah_test_fundal: wrong initial value'
call dev_free(a)
print '(A)', 'selah_test_fundal: ok'
endprogram selah_test_fundal
