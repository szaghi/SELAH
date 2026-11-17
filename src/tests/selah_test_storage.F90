!< SELAH test: the buffers each storage form needs, through `use selah` only.
program selah_test_storage
!< SELAH test: the buffers each storage form needs, through `use selah` only.
!<
!< The counts of every scheme in every form it has, -1 for a form it does not have, and ADAM's memory as the check:
!< its SSP path holds nrk stage buffers plus dq (ssp-54: 6 = s + 1), its low-storage path q_rk(1) plus dq (2).
!< Using only `use selah` also checks that the facade exports the whole of L0.
use selah
use selah_test_tools, only : check, failures
implicit none
type(butcher_tableau) :: bt
integer(I4P)          :: s

do s=1, size(RK_SCHEMES)
   bt = butcher(RK_SCHEMES(s))
   call check(registers_needed(RK_SCHEMES(s), FORM_BUTCHER) == bt%stages + 1, &
              trim(RK_SCHEMES(s))//': Butcher form needs stages + 1')
enddo
call check(registers_needed(SCHEME_SSP_54, FORM_BUTCHER) == 6, 'ssp-54 Butcher: 6, as ADAM (5 stage buffers + dq)')
call check(registers_needed(SCHEME_SSP_33, FORM_SHU_OSHER_2R) == 2, 'ssp-33 two-register: 2, as ADAM runge-kutta-3')
call check(registers_needed(SCHEME_EULER, FORM_SHU_OSHER_2R) == 2, 'euler two-register: 2')
call check(registers_needed(SCHEME_LSRK_54, FORM_WILLIAMSON_2N) == 2, 'lsrk-54 2N: 2')
call check(registers_needed(SCHEME_SSP_104, FORM_KETCHESON_2R) == 2, 'ssp-104 Ketcheson two-register: 2')
call check(registers_needed(SCHEME_SSP_104, FORM_BUTCHER) == 11, 'ssp-104 Butcher: 11 (why its 2-register form matters)')
do s=1, size(SPLITTING_SCHEMES)
   call check(registers_needed(SPLITTING_SCHEMES(s), FORM_SPLITTING) == 1, trim(SPLITTING_SCHEMES(s))//': 1 buffer')
enddo
call check(registers_needed(SCHEME_SSP_54, FORM_SHU_OSHER_2R) == -1, 'ssp-54 has no two-register form: -1')
call check(registers_needed(SCHEME_RK4, FORM_WILLIAMSON_2N) == -1, 'rk4 has no 2N form in the catalogue: -1')
call check(registers_needed(SCHEME_YOSHIDA_4, FORM_BUTCHER) == -1, 'a splitting has no Butcher form: -1')
call check(registers_needed('no-such-scheme', FORM_BUTCHER) == -1, 'unknown scheme: -1')
call check(registers_needed(SCHEME_SSP_54, 'no-such-form') == -1, 'unknown form: -1')

if (failures > 0) error stop 'selah_test_storage: failed'
print '(A)', 'selah_test_storage: ok'
endprogram selah_test_storage
