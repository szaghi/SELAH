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
