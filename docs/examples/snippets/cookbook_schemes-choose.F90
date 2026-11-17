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
