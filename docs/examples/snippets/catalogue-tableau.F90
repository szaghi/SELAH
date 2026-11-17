bt = butcher(SCHEME_SSP_33)
print '(/,A)', 'ssp-33, Butcher tableau (c | A, then b):'
do i=1, bt%stages
   print '(F6.3,A,3F7.3)', bt%c(i), ' |', bt%a(i,:)
enddo
print '(A,3F7.3)', '       ', bt%b
