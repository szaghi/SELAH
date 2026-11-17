   ! Loop of the box kernel step_update of SELAH (selah_kernels_box): d = sum_k coef(k) regs(k); q = q + dt d; with
   ! SELAH_WITH_DQ also dq = d, and dq = 0 where the mask is inactive (ADAM FNL's update_q order).
   ! Included unmasked, and once per layout with SELAH_ACTIVE the mask condition; SELAH_MASK_ARG adds the mask to the
   ! data clauses. Only the cells of the box are visited: dimension 1 innermost, the contiguous one.
   !$acc parallel loop independent collapse(5) DEVICEVAR(coef, regs, q SELAH_DQ_ARG SELAH_MASK_ARG) private(d)
   !$omp OMPLOOP collapse(5) DEVICEPTR(coef, regs, q SELAH_DQ_ARG SELAH_MASK_ARG) private(d)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
#ifdef SELAH_ACTIVE
                  if (SELAH_ACTIVE) then
#endif
                  d = 0._R_P
                  !$acc loop seq
                  do k=1, n
                     d = d + coef(k) * regs(i1,i2,i3,i4,i5,k)
                  enddo
                  q(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5) + dt * d
#ifdef SELAH_WITH_DQ
                  dq(i1,i2,i3,i4,i5) = d
#endif
#ifdef SELAH_ACTIVE
#ifdef SELAH_WITH_DQ
                  else
                     dq(i1,i2,i3,i4,i5) = 0._R_P
#endif
                  endif
#endif
               enddo
            enddo
         enddo
      enddo
   enddo
