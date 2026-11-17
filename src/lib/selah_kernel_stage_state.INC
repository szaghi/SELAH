   ! Loop of the box kernel stage_state of SELAH (selah_kernels_box): u = q + sum_k dt coef(k) regs(k): the stage
   ! state of a Butcher form (ADAM FNL's compute_stage: the private sum, then one store).
   ! Included unmasked, and once per layout with SELAH_ACTIVE the mask condition; SELAH_MASK_ARG adds the mask to the
   ! data clauses. Only the cells of the box are visited: dimension 1 innermost, the contiguous one.
   !$acc parallel loop independent collapse(5) DEVICEVAR(coef, q, regs, u SELAH_MASK_ARG) private(d)
   !$omp OMPLOOP collapse(5) DEVICEPTR(coef, q, regs, u SELAH_MASK_ARG) private(d)
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
                     d = d + dt * coef(k) * regs(i1,i2,i3,i4,i5,k)
                  enddo
                  u(i1,i2,i3,i4,i5) = q(i1,i2,i3,i4,i5) + d
#ifdef SELAH_ACTIVE
                  endif
#endif
               enddo
            enddo
         enddo
      enddo
   enddo
