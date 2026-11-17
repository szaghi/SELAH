do b=1, NB
   print '(A,I0,A,I2,A,I0,A)', 'block ', b, ': ', count(q(1,1:NI,1,1,b) /= q0(1,1:NI,1,1,b)), ' of ', NI, &
         ' cells changed'//trim(merge(' (solid)', '        ', b == 3))
enddo
