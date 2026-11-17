   ! Loop of the box kernel copy of SELAH (selah_kernels_box): y = x.
   ! Included unmasked, and once per layout with SELAH_ACTIVE the mask condition; SELAH_MASK_ARG adds the mask to the
   ! data clauses. Only the cells of the box are visited: dimension 1 innermost, the contiguous one.
   !$acc parallel loop independent collapse(5) DEVICEVAR(x, y SELAH_MASK_ARG)
   !$omp OMPLOOP collapse(5) DEVICEPTR(x, y SELAH_MASK_ARG)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
#ifdef SELAH_ACTIVE
                  if (SELAH_ACTIVE) then
#endif
                  y(i1,i2,i3,i4,i5) = x(i1,i2,i3,i4,i5)
#ifdef SELAH_ACTIVE
                  endif
#endif
               enddo
            enddo
         enddo
      enddo
   enddo
