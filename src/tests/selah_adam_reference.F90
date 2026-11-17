!< SELAH tests: ADAM's Runge–Kutta coefficients, verbatim, the oracle of the bitwise parity test.
module selah_adam_reference
!< SELAH tests: ADAM's Runge–Kutta coefficients, verbatim, the oracle of the bitwise parity test.
!<
!< `ssp_tableau` is copied verbatim from ADAM, `src/lib/common/adam_rk_object.F90` lines 72-159 at commit `ff879596`
!< (2026-11-17); `low_storage` holds the `runge-kutta-1/2/3` coefficients of `rk_object%initialize`, lines 450-464 of
!< the same file. R8P only, as in ADAM: excluded from the R4P builds (fobos, template-r4-tests).
use penf, only : I4P, R8P
implicit none
private
public :: RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
public :: ssp_tableau
public :: low_storage

character(len=13), parameter :: RK_1       ="runge-kutta-1"       !< Parameter of time scheme, Runge-Kutta 1.
character(len=13), parameter :: RK_2       ="runge-kutta-2"       !< Parameter of time scheme, Runge-Kutta 2.
character(len=13), parameter :: RK_3       ="runge-kutta-3"       !< Parameter of time scheme, Runge-Kutta 3.
character(len=18), parameter :: RK_SSP_11  ="runge-kutta-ssp-11"  !< Parameter of time scheme, Runge-Kutta SSP 11 (forward Euler).
character(len=18), parameter :: RK_SSP_22  ="runge-kutta-ssp-22"  !< Parameter of time scheme, Runge-Kutta SSP 22.
character(len=18), parameter :: RK_SSP_33  ="runge-kutta-ssp-33"  !< Parameter of time scheme, Runge-Kutta SSP 33.
character(len=18), parameter :: RK_SSP_54  ="runge-kutta-ssp-54"  !< Parameter of time scheme, Runge-Kutta SSP 54.

contains
   ! ---- verbatim from ADAM (ff879596), src/lib/common/adam_rk_object.F90:72-159 ----
   pure subroutine ssp_tableau(scheme, nrk, alph, beta, gamm)
   !< Return the Butcher tableau of an SSP Runge-Kutta scheme (issue #104): `alph` the matrix A (strictly lower
   !< triangular), `beta` the weights b, `gamm` the abscissas c. The one source of the SSP tableaux: `rk_object` and the
   !< PRISM boundary and PIC integrators call it. A scheme that is not SSP leaves `nrk` unchanged and the arrays
   !< unallocated.
   character(*),           intent(in)    :: scheme    !< Runge-Kutta scheme.
   integer(I4P),           intent(inout) :: nrk       !< Runge-Kutta stages number.
   real(R8P), allocatable, intent(out)   :: alph(:,:) !< Butcher matrix A.
   real(R8P), allocatable, intent(out)   :: beta(:)   !< Butcher weights b.
   real(R8P), allocatable, intent(out)   :: gamm(:)   !< Butcher abscissas c.
   real(R8P)                             :: so_a(5,0:4) !< SSP-54 Shu-Osher alpha.
   real(R8P)                             :: so_b(5,0:4) !< SSP-54 Shu-Osher beta.
   real(R8P)                             :: bk(0:5,5)   !< SSP-54 stage coefficients: u(i) = u^n + dt sum_j bk(i,j) F_j.
   integer(I4P)                          :: i, k        !< Counters.

   select case(scheme)
   case(RK_SSP_11) ! 1 stage, 1st order SSP (forward Euler): the 1-stage debugging
                   ! instrument of the staged forest path (issue #25); alph = Butcher A = 0,
                   ! beta = Butcher b = 1, gamm = Butcher c = 0.
      nrk = 1
      allocate(alph(nrk,nrk), beta(nrk), gamm(nrk))
      alph = 0._R8P
      beta = 0._R8P
      gamm = 0._R8P

      beta(1) = 1._R8P
   case(RK_SSP_22) ! 2 stages, 2nd order SSP
      nrk = 2
      allocate(alph(nrk,nrk), beta(nrk), gamm(nrk))
      alph = 0._R8P
      beta = 0._R8P
      gamm = 0._R8P

      beta(1) = 0.5_R8P
      beta(2) = 0.5_R8P

      alph(2,1) = 1._R8P

      gamm(2) = 1._R8P
   case(RK_SSP_33) ! 3 stages, 3rd order SSP
      nrk = 3 ! was 2 (issue #25 finding): beta(3)/alph(3,:)/gamm(3) below wrote OUT OF
                   ! BOUNDS of nrk=2 allocations -- silent heap corruption for any SSP-33 user
      allocate(alph(nrk,nrk), beta(nrk), gamm(nrk))
      alph = 0._R8P
      beta = 0._R8P
      gamm = 0._R8P

      beta(1) = 1._R8P/6._R8P
      beta(2) = 1._R8P/6._R8P
      beta(3) = 2._R8P/3._R8P

      alph(2,1) = 1._R8P
      alph(3,1) = 0.25_R8P ; alph(3,2) = 0.25_R8P

      gamm(2) = 1._R8P
      gamm(3) = 0.5_R8P
   case(RK_SSP_54) ! 5 stages, 4th order SSP: Spiteri & Ruuth (2002, SIAM J. Numer. Anal. 40, 469), SSP coeff. 1.508
      ! The Shu-Osher form u(i) = sum_k a_ik u(k) + dt b_ik F(u(k)), u(0) = u^n, u^(n+1) = u(5), as printed (15 digits)
      ! by Isherwood, Grant & Gottlieb (2018, arXiv:1810.04800, section 4); every a_ik, b_ik >= 0 (the SSP property) and
      ! sum_k a_ik = 1. The Butcher tableau is derived from it, row i of the stage coefficients being
      ! sum_k a_ik (row k) + b_ik e_(k+1): the order conditions then hold to ~5e-16 (issue #104; the former 14-digit
      ! Butcher literals were accurate to ~1e-10 only, their weights summing to 1 - 8.8e-11).
      nrk = 5
      allocate(alph(nrk,nrk), beta(nrk), gamm(nrk))
      so_a = 0._R8P
      so_b = 0._R8P
      so_a(1,0) = 1._R8P
      so_b(1,0) = 0.391752226571890_R8P
      so_a(2,0) = 0.444370493651235_R8P ; so_a(2,1) = 0.555629506348765_R8P
      so_b(2,1) = 0.368410593050371_R8P
      so_a(3,0) = 0.620101851488403_R8P ; so_a(3,2) = 0.379898148511597_R8P
      so_b(3,2) = 0.251891774271694_R8P
      so_a(4,0) = 0.178079954393132_R8P ; so_a(4,3) = 0.821920045606868_R8P
      so_b(4,3) = 0.544974750228521_R8P
      so_a(5,2) = 0.517231671970585_R8P ; so_a(5,3) = 0.096059710526147_R8P ; so_a(5,4) = 0.386708617503268_R8P
      so_b(5,3) = 0.063692468666290_R8P ; so_b(5,4) = 0.226007483236906_R8P
      bk = 0._R8P
      do i=1, 5
         do k=0, i - 1
            bk(i,:) = bk(i,:) + so_a(i,k) * bk(k,:)
            bk(i,k+1) = bk(i,k+1) + so_b(i,k)
         enddo
      enddo
      alph = bk(0:4,:)
      beta = bk(5,:)
      gamm = sum(alph, dim=2)
   endselect
   endsubroutine ssp_tableau
   ! ---- end of the verbatim copy ----

   pure subroutine low_storage(scheme, nrk, ark, brk, crk)
   !< ADAM's `runge-kutta-1/2/3` coefficients, from `rk_object%initialize` (ff879596, adam_rk_object.F90:450-464).
   character(*),           intent(in)  :: scheme !< Runge-Kutta scheme.
   integer(I4P),           intent(out) :: nrk    !< Stages number.
   real(R8P), allocatable, intent(out) :: ark(:) !< Coefficient of q^n.
   real(R8P), allocatable, intent(out) :: brk(:) !< Coefficient of the current stage.
   real(R8P), allocatable, intent(out) :: crk(:) !< Coefficient of dt times the residual.

   nrk = 0
   select case(scheme)
   case(RK_1) ! 1 stage, 1st order, Euler
      nrk = 1
      allocate(ark(nrk), brk(nrk), crk(nrk))
      ark(1) = 1._R8P ; brk(1) = 0._R8P ; crk(1) = 1._R8P
   case(RK_2) ! 2 stages, low storage, 2nd order TVD
      nrk = 2
      allocate(ark(nrk), brk(nrk), crk(nrk))
      ark(1) = 1._R8P  ; brk(1) = 0._R8P  ; crk(1) = 1._R8P
      ark(2) = 0.5_R8P ; brk(2) = 0.5_R8P ; crk(2) = 0.5_R8P
   case(RK_3) ! 3 stages, low storage, 3rd order TVD
      nrk = 3
      allocate(ark(nrk), brk(nrk), crk(nrk))
      ark(1) = 1._R8P        ; brk(1) = 0._R8P        ; crk(1) = 1._R8P
      ark(2) = 0.75_R8P      ; brk(2) = 0.25_R8P      ; crk(2) = 0.25_R8P
      ark(3) = 1._R8P/3._R8P ; brk(3) = 2._R8P/3._R8P ; crk(3) = 2._R8P/3._R8P
   endselect
   endsubroutine low_storage
endmodule selah_adam_reference
