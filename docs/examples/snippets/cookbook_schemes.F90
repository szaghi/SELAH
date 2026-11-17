program cookbook_schemes
!< Cookbook, schemes: choosing a scheme and its form; the SSP time step.
use selah
implicit none
character(16) :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('choose') ; call choose
case('ssp-dt') ; call ssp_dt
case default   ; error stop 'usage: cookbook_schemes choose|ssp-dt'
endselect

contains
   subroutine choose
   !< For each scheme: its order, the fewest registers any of its forms needs, and in which form.
   character(14), parameter :: FORMS(4) = [character(14) :: FORM_BUTCHER, FORM_SHU_OSHER_2R, FORM_WILLIAMSON_2N, &
                                                             FORM_KETCHESON_2R]
   type(butcher_tableau)    :: bt
   integer(I4P)             :: s, f, n, best
   character(14)            :: form

   print '(A)', 'scheme   order  SSP   fewest registers  in form'
   do s=1, size(RK_SCHEMES)
      bt = butcher(trim(RK_SCHEMES(s)))
      best = huge(1_I4P)
      do f=1, size(FORMS)
         n = registers_needed(trim(RK_SCHEMES(s)), trim(FORMS(f)))
         if (n > 0 .and. n < best) then
            best = n ; form = FORMS(f)
         endif
      enddo
      print '(A8,I6,A6,I12,6X,A)', RK_SCHEMES(s), bt%order, merge('  yes', '   no', bt%ssp_coefficient > 0._R_P), &
            best, trim(form)
   enddo
   endsubroutine choose

   subroutine ssp_dt
   !< The largest strong-stability-preserving step of each SSP scheme for a spatial operator whose forward-Euler
   !< step is TVD up to dt_FE, and the residual evaluations it costs per unit time.
   real(R_P), parameter  :: DT_FE = 1.e-3_R_P   ! the forward-Euler limit of your spatial operator (CFL_FE dx / a)
   type(butcher_tableau) :: bt
   integer(I4P)          :: s

   print '(A)', 'scheme   C (SSP coeff)  dt_max = C dt_FE  residuals per unit time  C / stages'
   do s=1, size(RK_SCHEMES)
      bt = butcher(trim(RK_SCHEMES(s)))
      if (bt%ssp_coefficient <= 0._R_P) cycle   ! not SSP: no strong-stability guarantee
      print '(A8,F12.4,ES17.3,F21.1,F15.3)', RK_SCHEMES(s), bt%ssp_coefficient, bt%ssp_coefficient * DT_FE, &
            bt%stages / (bt%ssp_coefficient * DT_FE), bt%ssp_coefficient / bt%stages
   enddo
   endsubroutine ssp_dt
endprogram cookbook_schemes
