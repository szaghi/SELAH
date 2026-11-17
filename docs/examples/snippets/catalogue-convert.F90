bt2 = butcher_from_shu_osher_2r(shu_osher_2r(SCHEME_SSP_33))      ! two-register coefficients -> Butcher
print '(A,ES9.2)', 'two-register form of ssp-33, back to Butcher: max |difference| = ', &
      max(maxval(abs(bt2%a - bt%a)), maxval(abs(bt2%b - bt%b)))
print '(A,L1)', 'explicit: ', is_explicit(bt)
