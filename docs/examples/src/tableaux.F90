!run tableaux tableaux
program tableaux
!< The Butcher tableau of every scheme of the catalogue, printed from the library (the Schemes page includes it).
use selah
implicit none
type(butcher_tableau)       :: bt
type(splitting_composition) :: sc
integer(I4P)                :: s, i

do s=1, size(RK_SCHEMES)
   bt = butcher(trim(RK_SCHEMES(s)))
   print '(A,I0,A,I0,A)', trim(RK_SCHEMES(s))//': ', bt%stages, ' stages, order ', bt%order, ' (c | A, then b)'
   do i=1, bt%stages
      print '(F10.6,A,*(F10.6))', bt%c(i), ' |', bt%a(i,1:max(i-1,1))
   enddo
   print '(A,*(F10.6))', '           ', bt%b
   print '(A)', ''
enddo
do s=1, size(SPLITTING_SCHEMES)
   sc = splitting(trim(SPLITTING_SCHEMES(s)))
   print '(A,I0,A,I0,A)', trim(SPLITTING_SCHEMES(s))//': ', sc%flows, ' flows, order ', sc%order, &
                          ' (partition, fraction of the step)'
   print '(*(I3,F10.6,:,/))', (sc%partition(i), sc%coefficient(i), i=1, sc%flows)
   print '(A)', ''
enddo
endprogram tableaux
