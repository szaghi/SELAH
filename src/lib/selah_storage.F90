!< SELAH, L0: the storage forms of the schemes and the field-sized buffers each one needs.
module selah_storage
!< SELAH, L0: the storage forms of the schemes and the field-sized buffers each one needs.
!<
!< `registers_needed(scheme, form)` counts the field-sized buffers an integration step needs besides the solution
!< itself, **including the buffer the residual is written to**. It is known before anything is allocated, so a solver
!< can size its memory budget first (ADAM sizes its block pool this way). For comparison, ADAM's SSP path holds `nrk`
!< stage buffers plus `dq`, which is the Butcher count s + 1, and its low-storage path holds `q_rk(1)` plus `dq`, which
!< is the two-register count 2.
!<
!<| Form | Buffers | Schemes | What they hold |
!<|------|---------|---------|----------------|
!<| `FORM_BUTCHER`       | s + 1 | every RK scheme | the stage state and the s stage residuals |
!<| `FORM_SHU_OSHER_2R`  | 2     | euler, ssp-22, ssp-33 | u^n and the residual (the current stage is the solution) |
!<| `FORM_WILLIAMSON_2N` | 2     | lsrk-33, lsrk-54 | Delta Q and the residual |
!<| `FORM_KETCHESON_2R`  | 2     | ssp-104 | the second register q2 and the residual |
!<| `FORM_SPLITTING`     | 1     | every splitting scheme | the residual of the partition being advanced |
!<
!< A scheme that has no such form gives -1.
use penf, only : I4P
use selah_tableaux, only : butcher_tableau
use selah_schemes
use selah_splittings
implicit none
private
public :: FORM_BUTCHER, FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, FORM_KETCHESON_2R, FORM_SPLITTING
public :: registers_needed

character(*), parameter :: FORM_BUTCHER       = 'butcher'       !< Stage state plus every stage residual.
character(*), parameter :: FORM_SHU_OSHER_2R  = 'shu-osher-2r'  !< Two-register Shu–Osher (u^n and the residual).
character(*), parameter :: FORM_WILLIAMSON_2N = 'williamson-2n' !< Williamson 2N (Delta Q and the residual).
character(*), parameter :: FORM_KETCHESON_2R  = 'ketcheson-2r'  !< Ketcheson's two-register SSPRK(10,4).
character(*), parameter :: FORM_SPLITTING     = 'splitting'     !< One partition residual.

contains
   pure function registers_needed(scheme, form) result(n)
   !< Field-sized buffers a step of `scheme` in `form` needs besides the solution, the residual buffer included.
   character(*), intent(in) :: scheme !< Scheme name.
   character(*), intent(in) :: form   !< Storage form.
   integer(I4P)             :: n      !< Buffers, -1 if the scheme has no such form.
   type(butcher_tableau)    :: bt     !< Tableau of an RK scheme.
   character(:), allocatable :: name  !< Trimmed scheme name.

   n = -1_I4P
   name = trim(adjustl(scheme))
   select case(trim(adjustl(form)))
   case(FORM_BUTCHER)
      bt = butcher(name)
      if (bt%stages > 0) n = bt%stages + 1_I4P
   case(FORM_SHU_OSHER_2R)
      if (any(name == [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33])) n = 2_I4P
   case(FORM_WILLIAMSON_2N)
      if (any(name == [character(7) :: SCHEME_LSRK_33, SCHEME_LSRK_54])) n = 2_I4P
   case(FORM_KETCHESON_2R)
      if (name == SCHEME_SSP_104) n = 2_I4P
   case(FORM_SPLITTING)
      if (any(name == SPLITTING_SCHEMES)) n = 1_I4P
   endselect
   endfunction registers_needed
endmodule selah_storage
