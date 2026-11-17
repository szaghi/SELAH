print '(/,A)', 'splitting        flows order  sum of the coefficients of partitions 1, 2'
do s=1, size(SPLITTING_SCHEMES)
   sc = splitting(trim(SPLITTING_SCHEMES(s)))
   print '(A15,I7,I6,2F9.4)', SPLITTING_SCHEMES(s), sc%flows, sc%order, &
         sum(sc%coefficient, mask=sc%partition == 1), sum(sc%coefficient, mask=sc%partition == 2)
enddo
