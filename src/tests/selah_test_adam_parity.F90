!< SELAH test: the SSP tableaux and the two-register coefficients are bitwise those of ADAM.
program selah_test_adam_parity
!< SELAH test: the SSP tableaux and the two-register coefficients are bitwise those of ADAM.
!<
!< The oracle is ADAM's code, verbatim (selah_adam_reference). Bitwise means `==` on every coefficient: SELAH must be a
!< drop-in source of ADAM's tableaux, so that wiring it into ADAM (notes/PLAN.md, M7) changes no golden output. R8P
!< builds only.
use penf, only : I4P, R8P
use selah_kind_parameters, only : R_P
use selah_tableaux
use selah_schemes
use selah_adam_reference
use selah_test_tools, only : check, failures
implicit none
character(18), parameter :: ADAM_SSP(4)  = [character(18) :: RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54]
character(7),  parameter :: SELAH_SSP(4) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54]
character(13), parameter :: ADAM_LS(3)   = [character(13) :: RK_1, RK_2, RK_3]
character(7),  parameter :: SELAH_LS(3)  = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33]
type(butcher_tableau)    :: bt
type(shu_osher_2r_form)  :: r
real(R8P), allocatable   :: alph(:,:), beta(:), gamm(:), ark(:), brk(:), crk(:)
integer(I4P)             :: nrk, s

if (R_P /= R8P) error stop 'selah_test_adam_parity: R8P builds only'
do s=1, size(ADAM_SSP)
   nrk = 0
   call ssp_tableau(scheme=ADAM_SSP(s), nrk=nrk, alph=alph, beta=beta, gamm=gamm)
   bt = butcher(SELAH_SSP(s))
   call check(bt%stages == nrk, trim(SELAH_SSP(s))//' vs '//trim(ADAM_SSP(s))//': stages')
   call check(all(bt%a == alph) .and. all(bt%b == beta) .and. all(bt%c == gamm), &
              trim(SELAH_SSP(s))//' vs '//trim(ADAM_SSP(s))//': A, b, c bitwise equal')
enddo
do s=1, size(ADAM_LS)
   call low_storage(scheme=ADAM_LS(s), nrk=nrk, ark=ark, brk=brk, crk=crk)
   r = shu_osher_2r(SELAH_LS(s))
   call check(r%stages == nrk .and. all(r%a == ark) .and. all(r%b == brk) .and. all(r%c == crk), &
              trim(SELAH_LS(s))//' two-register vs '//trim(ADAM_LS(s))//' ark, brk, crk: bitwise equal')
enddo

if (failures > 0) error stop 'selah_test_adam_parity: failed'
print '(A)', 'selah_test_adam_parity: ok'
endprogram selah_test_adam_parity
