   ! Loop of the box kernel williamson_2n_stage of SELAH (selah_kernels_box): dqr = a dqr + dt r; q = q + b dqr: a
   ! stage of a Williamson 2N form.
   ! Included unmasked, and once per layout with SELAH_ACTIVE the mask condition; SELAH_MASK_ARG adds the mask to the
   ! data clauses. Only the cells of the box are visited: dimension 1 innermost, the contiguous one. With SELAH_A_ZERO
   ! (a = 0, the first stage of every 2N scheme) dqr is not read: dqr = dt r, so that a register holding garbage or NaN
   ! (0 * NaN = NaN) cannot leak into the step, and one field read is saved.
   !$acc parallel loop independent collapse(5) DEVICEVAR(r, dqr, q SELAH_MASK_ARG)
   !$omp OMPLOOP collapse(5) DEVICEPTR(r, dqr, q SELAH_MASK_ARG)
   do i5=lo5, hi5
      do i4=lo4, hi4
         do i3=lo3, hi3
            do i2=lo2, hi2
               do i1=lo1, hi1
#ifdef SELAH_ACTIVE
                  if (SELAH_ACTIVE) then
#endif
#ifdef SELAH_A_ZERO
                  dqr(i1,i2,i3,i4,i5) = dt * r(i1,i2,i3,i4,i5)
#else
                  dqr(i1,i2,i3,i4,i5) = a * dqr(i1,i2,i3,i4,i5) + dt * r(i1,i2,i3,i4,i5)
#endif
                  q(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5) + b * dqr(i1,i2,i3,i4,i5)
#ifdef SELAH_ACTIVE
                  endif
#endif
               enddo
            enddo
         enddo
      enddo
   enddo
