print '(A)', 'scheme   stages order  SSP coeff  butcher two-reg 2N  ketcheson'
do s=1, size(RK_SCHEMES)
   bt = butcher(trim(RK_SCHEMES(s)))                             ! the names are blank-padded: trim them
   print '(A8,I6,I6,F11.4,4I8)', RK_SCHEMES(s), bt%stages, bt%order, bt%ssp_coefficient, &
         registers_needed(trim(RK_SCHEMES(s)), FORM_BUTCHER), registers_needed(trim(RK_SCHEMES(s)), FORM_SHU_OSHER_2R), &
         registers_needed(trim(RK_SCHEMES(s)), FORM_WILLIAMSON_2N), registers_needed(trim(RK_SCHEMES(s)), FORM_KETCHESON_2R)
enddo
