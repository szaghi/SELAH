!< SELAH, L0: the Runge–Kutta scheme catalogue.
module selah_schemes
!< SELAH, L0: the Runge–Kutta scheme catalogue.
!<
!< Each scheme is stored in the form its source publishes, and its other forms are derived from that one:
!<
!<| Scheme | Stages | Order | SSP coefficient | Primary form | Other forms | Source |
!<|--------|--------|-------|-----------------|--------------|-------------|--------|
!<| `euler`   | 1  | 1 | 1     | Butcher | Shu–Osher, two-register | — |
!<| `ssp-22`  | 2  | 2 | 1     | Butcher | Shu–Osher, two-register | Shu & Osher (1988) |
!<| `ssp-33`  | 3  | 3 | 1     | Butcher | Shu–Osher, two-register | Shu & Osher (1988) |
!<| `ssp-54`  | 5  | 4 | 1.508 | Shu–Osher | Butcher | Spiteri & Ruuth (2002), digits of Isherwood, Grant & Gottlieb (2018) |
!<| `ssp-104` | 10 | 4 | 6     | Shu–Osher | Butcher | Ketcheson (2008) |
!<| `rk4`     | 4  | 4 | 0     | Butcher | — | Kutta (1901) |
!<| `lsrk-33` | 3  | 3 | 0     | Williamson 2N | Butcher | Williamson (1980) |
!<| `lsrk-54` | 5  | 4 | 0     | Williamson 2N | Butcher | Carpenter & Kennedy (1994) |
!<
!< The SSP coefficient is computed from the Shu–Osher form (`ssp_coefficient`), never typed in. The Butcher literals
!< of `euler`, `ssp-22` and `ssp-33` and the Shu–Osher digits of `ssp-54` are written as in ADAM's `ssp_tableau`, so
!< the two give the same tableaux bitwise. An unknown name returns a form of 0 stages.
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_tableaux
implicit none
private
public :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54, SCHEME_SSP_104
public :: SCHEME_RK4, SCHEME_LSRK_33, SCHEME_LSRK_54
public :: RK_SCHEMES
public :: butcher
public :: shu_osher
public :: shu_osher_2r
public :: williamson_2n

character(*), parameter :: SCHEME_EULER   = 'euler'   !< Forward Euler.
character(*), parameter :: SCHEME_SSP_22  = 'ssp-22'  !< SSP, 2 stages, order 2 (Heun).
character(*), parameter :: SCHEME_SSP_33  = 'ssp-33'  !< SSP, 3 stages, order 3.
character(*), parameter :: SCHEME_SSP_54  = 'ssp-54'  !< SSP, 5 stages, order 4.
character(*), parameter :: SCHEME_SSP_104 = 'ssp-104' !< SSP, 10 stages, order 4, low storage.
character(*), parameter :: SCHEME_RK4     = 'rk4'     !< Classic 4-stage, order 4.
character(*), parameter :: SCHEME_LSRK_33 = 'lsrk-33' !< Williamson 2N, 3 stages, order 3.
character(*), parameter :: SCHEME_LSRK_54 = 'lsrk-54' !< Carpenter–Kennedy 2N, 5 stages, order 4.
character(7), parameter :: RK_SCHEMES(8) = [character(7) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, &
                                            SCHEME_SSP_54, SCHEME_SSP_104, SCHEME_RK4, SCHEME_LSRK_33, &
                                            SCHEME_LSRK_54]                             !< Every RK scheme.

contains
   pure function butcher(scheme) result(bt)
   !< Butcher tableau of a scheme.
   character(*), intent(in) :: scheme !< Scheme name.
   type(butcher_tableau)    :: bt     !< Tableau (0 stages if the name is unknown).

   select case(trim(adjustl(scheme)))
   case(SCHEME_EULER)
      call literal(1_I4P, 1_I4P)
      bt%b(1) = 1._R_P
      bt%ssp_coefficient = ssp_coefficient(shu_osher(scheme))
   case(SCHEME_SSP_22)
      call literal(2_I4P, 2_I4P)
      bt%b(1) = 0.5_R_P
      bt%b(2) = 0.5_R_P
      bt%a(2,1) = 1._R_P
      bt%c(2) = 1._R_P
      bt%ssp_coefficient = ssp_coefficient(shu_osher(scheme))
   case(SCHEME_SSP_33)
      call literal(3_I4P, 3_I4P)
      bt%b(1) = 1._R_P/6._R_P
      bt%b(2) = 1._R_P/6._R_P
      bt%b(3) = 2._R_P/3._R_P
      bt%a(2,1) = 1._R_P
      bt%a(3,1) = 0.25_R_P ; bt%a(3,2) = 0.25_R_P
      bt%c(2) = 1._R_P
      bt%c(3) = 0.5_R_P
      bt%ssp_coefficient = ssp_coefficient(shu_osher(scheme))
   case(SCHEME_SSP_54, SCHEME_SSP_104)
      bt = butcher_from_shu_osher(shu_osher(scheme), order=4_I4P)
   case(SCHEME_RK4)
      call literal(4_I4P, 4_I4P)
      bt%a(2,1) = 0.5_R_P
      bt%a(3,2) = 0.5_R_P
      bt%a(4,3) = 1._R_P
      bt%b = [1._R_P/6._R_P, 1._R_P/3._R_P, 1._R_P/3._R_P, 1._R_P/6._R_P]
      bt%c = [0._R_P, 0.5_R_P, 0.5_R_P, 1._R_P]
   case(SCHEME_LSRK_33)
      bt = butcher_from_williamson_2n(williamson_2n(scheme), order=3_I4P)
   case(SCHEME_LSRK_54)
      bt = butcher_from_williamson_2n(williamson_2n(scheme), order=4_I4P)
   endselect
   contains
      pure subroutine literal(s, p)
      !< Allocate a zero tableau of s stages and declared order p.
      integer(I4P), intent(in) :: s !< Stages.
      integer(I4P), intent(in) :: p !< Order.

      bt%stages = s
      bt%order = p
      allocate(bt%a(s,s), bt%b(s), bt%c(s))
      bt%a = 0._R_P
      bt%b = 0._R_P
      bt%c = 0._R_P
      endsubroutine literal
   endfunction butcher

   pure function shu_osher(scheme) result(so)
   !< Shu–Osher form of an SSP scheme (euler, ssp-22, ssp-33, ssp-54, ssp-104).
   character(*), intent(in) :: scheme !< Scheme name.
   type(shu_osher_form)     :: so     !< Shu–Osher form (0 stages if the scheme has none).
   integer(I4P)             :: i      !< Counter.

   select case(trim(adjustl(scheme)))
   case(SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33)
      so = shu_osher_from_2r(shu_osher_2r(scheme))
   case(SCHEME_SSP_54)
      ! Spiteri & Ruuth (2002, SIAM J. Numer. Anal. 40, 469), SSP coefficient 1.508; the 15 digits printed by
      ! Isherwood, Grant & Gottlieb (2018, arXiv:1810.04800, section 4), as in ADAM's ssp_tableau. Every coefficient is
      ! non-negative and each row of alpha sums to 1.
      call zero(5_I4P)
      so%alpha(2,1) = 1._R_P
      so%beta(2,1) = 0.391752226571890_R_P
      so%alpha(3,1) = 0.444370493651235_R_P ; so%alpha(3,2) = 0.555629506348765_R_P
      so%beta(3,2) = 0.368410593050371_R_P
      so%alpha(4,1) = 0.620101851488403_R_P ; so%alpha(4,3) = 0.379898148511597_R_P
      so%beta(4,3) = 0.251891774271694_R_P
      so%alpha(5,1) = 0.178079954393132_R_P ; so%alpha(5,4) = 0.821920045606868_R_P
      so%beta(5,4) = 0.544974750228521_R_P
      so%alpha(6,3) = 0.517231671970585_R_P ; so%alpha(6,4) = 0.096059710526147_R_P
      so%alpha(6,5) = 0.386708617503268_R_P
      so%beta(6,4) = 0.063692468666290_R_P ; so%beta(6,5) = 0.226007483236906_R_P
   case(SCHEME_SSP_104)
      ! Ketcheson (2008, SIAM J. Sci. Comput. 30, 2113), SSPRK(10,4), SSP coefficient 6. Its low-storage recipe, with
      ! u^n = Y1: five forward-Euler steps of dt/6; q2 = (u^n + 9 q1)/25 and q1 = 15 q2 - 5 q1; four more Euler steps;
      ! u^{n+1} = q2 + 3/5 q1 + dt/10 F(q1). Written as Shu–Osher stages, the stage after the fifth Euler step is
      ! Y6 = 3/5 Y1 + 2/5 Y5 + dt/15 F(Y5) and u^{n+1} = Y11 = 1/25 Y1 + 9/25 Y5 + 3/5 Y10 + dt (3/50 F(Y5) + 1/10 F(Y10)).
      call zero(10_I4P)
      do i=2, 5
         so%alpha(i,i-1) = 1._R_P
         so%beta(i,i-1) = 1._R_P/6._R_P
      enddo
      so%alpha(6,1) = 3._R_P/5._R_P ; so%alpha(6,5) = 2._R_P/5._R_P
      so%beta(6,5) = 1._R_P/15._R_P
      do i=7, 10
         so%alpha(i,i-1) = 1._R_P
         so%beta(i,i-1) = 1._R_P/6._R_P
      enddo
      so%alpha(11,1) = 1._R_P/25._R_P ; so%alpha(11,5) = 9._R_P/25._R_P ; so%alpha(11,10) = 3._R_P/5._R_P
      so%beta(11,5) = 3._R_P/50._R_P ; so%beta(11,10) = 1._R_P/10._R_P
   endselect
   contains
      pure subroutine zero(s)
      !< Allocate a zero Shu–Osher form of s stages.
      integer(I4P), intent(in) :: s !< Stages.

      so%stages = s
      allocate(so%alpha(s+1,s), so%beta(s+1,s))
      so%alpha = 0._R_P
      so%beta = 0._R_P
      endsubroutine zero
   endfunction shu_osher

   pure function shu_osher_2r(scheme) result(r)
   !< Two-register Shu–Osher form (euler, ssp-22, ssp-33): ADAM's `runge-kutta-1/2/3` ark, brk, crk.
   character(*), intent(in) :: scheme !< Scheme name.
   type(shu_osher_2r_form)  :: r      !< Two-register form (0 stages if the scheme has none).

   select case(trim(adjustl(scheme)))
   case(SCHEME_EULER)
      r = shu_osher_2r_form(stages=1_I4P, a=[1._R_P], b=[0._R_P], c=[1._R_P])
   case(SCHEME_SSP_22)
      r = shu_osher_2r_form(stages=2_I4P, a=[1._R_P, 0.5_R_P], b=[0._R_P, 0.5_R_P], c=[1._R_P, 0.5_R_P])
   case(SCHEME_SSP_33)
      r = shu_osher_2r_form(stages=3_I4P, a=[1._R_P, 0.75_R_P, 1._R_P/3._R_P], b=[0._R_P, 0.25_R_P, 2._R_P/3._R_P], &
                            c=[1._R_P, 0.25_R_P, 2._R_P/3._R_P])
   endselect
   endfunction shu_osher_2r

   pure function williamson_2n(scheme) result(w)
   !< Williamson 2N form (lsrk-33, lsrk-54).
   character(*), intent(in) :: scheme !< Scheme name.
   type(williamson_2n_form) :: w      !< 2N form (0 stages if the scheme has none).

   select case(trim(adjustl(scheme)))
   case(SCHEME_LSRK_33)
      ! Williamson (1980, J. Comput. Phys. 35, 48)
      w = williamson_2n_form(stages=3_I4P, a=[0._R_P, -5._R_P/9._R_P, -153._R_P/128._R_P], &
                             b=[1._R_P/3._R_P, 15._R_P/16._R_P, 8._R_P/15._R_P])
   case(SCHEME_LSRK_54)
      ! Carpenter & Kennedy (1994, NASA TM-109112), the five-stage fourth-order 2N scheme (solution 3)
      w = williamson_2n_form(stages=5_I4P,                                                                     &
                             a=[0._R_P,                                       -567301805773._R_P/1357537059087._R_P, &
                                -2404267990393._R_P/2016746695238._R_P,      -3550918686646._R_P/2091501179385._R_P, &
                                -1275806237668._R_P/842570457699._R_P],                                      &
                             b=[1432997174477._R_P/9575080441755._R_P,       5161836677717._R_P/13612068292357._R_P, &
                                1720146321549._R_P/2090206949498._R_P,       3134564353537._R_P/4481467310338._R_P,  &
                                2277821191437._R_P/14882151754819._R_P])
   endselect
   endfunction williamson_2n
endmodule selah_schemes
