program catalogue
!< Tutorial, chapter 2: the catalogue of schemes, their storage forms and their coefficients.
use selah
implicit none
type(butcher_tableau)       :: bt, bt2
type(splitting_composition) :: sc
integer(I4P)                :: s, i

print '(A)', 'scheme   stages order  SSP coeff  butcher two-reg 2N  ketcheson'
do s=1, size(RK_SCHEMES)
   bt = butcher(trim(RK_SCHEMES(s)))                             ! the names are blank-padded: trim them
   print '(A8,I6,I6,F11.4,4I8)', RK_SCHEMES(s), bt%stages, bt%order, bt%ssp_coefficient, &
         registers_needed(trim(RK_SCHEMES(s)), FORM_BUTCHER), registers_needed(trim(RK_SCHEMES(s)), FORM_SHU_OSHER_2R), &
         registers_needed(trim(RK_SCHEMES(s)), FORM_WILLIAMSON_2N), registers_needed(trim(RK_SCHEMES(s)), FORM_KETCHESON_2R)
enddo
bt = butcher(SCHEME_SSP_33)
print '(/,A)', 'ssp-33, Butcher tableau (c | A, then b):'
do i=1, bt%stages
   print '(F6.3,A,3F7.3)', bt%c(i), ' |', bt%a(i,:)
enddo
print '(A,3F7.3)', '       ', bt%b
bt2 = butcher_from_shu_osher_2r(shu_osher_2r(SCHEME_SSP_33))      ! two-register coefficients -> Butcher
print '(A,ES9.2)', 'two-register form of ssp-33, back to Butcher: max |difference| = ', &
      max(maxval(abs(bt2%a - bt%a)), maxval(abs(bt2%b - bt%b)))
print '(A,L1)', 'explicit: ', is_explicit(bt)
print '(/,A)', 'splitting        flows order  sum of the coefficients of partitions 1, 2'
do s=1, size(SPLITTING_SCHEMES)
   sc = splitting(trim(SPLITTING_SCHEMES(s)))
   print '(A15,I7,I6,2F9.4)', SPLITTING_SCHEMES(s), sc%flows, sc%order, &
         sum(sc%coefficient, mask=sc%partition == 1), sum(sc%coefficient, mask=sc%partition == 2)
enddo
endprogram catalogue
