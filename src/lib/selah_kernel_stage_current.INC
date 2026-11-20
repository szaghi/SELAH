   ! Loop of the box kernel stage_current of SELAH (selah_kernels_box): y = (h - sum_{r<n} coef(r) regs(r)) / coef(n),
   ! the stage current of a charge-conserving scheme (selah_stepper `stage_current`): the private sum, then one store.
   !$acc parallel loop independent collapse(5) DEVICEVAR(coef, h, regs, y) private(d)
   !$omp OMPLOOP collapse(5) DEVICEPTR(coef, h, regs, y) private(d)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
                  d = h(i1,i2,i3,i4,i5)
                  !$acc loop seq
                  do k=1, n-1
                     d = d - coef(k) * regs(i1,i2,i3,i4,i5,k)
                  enddo
                  y(i1,i2,i3,i4,i5) = d / coef(n)
               enddo
            enddo
         enddo
      enddo
   enddo
