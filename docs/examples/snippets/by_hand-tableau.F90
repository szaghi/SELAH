bt = butcher(SCHEME_SSP_54)
allocate(regs(N, bt%stages))
print '(A,I0,A,I0)', 'ssp-54: stages ', bt%stages, ', order ', bt%order
