   ! Setup of every box kernel of SELAH (selah_kernels_box): the loop range, and the mask component.
   lo1 = bx%lo(1) ; lo2 = bx%lo(2) ; lo3 = bx%lo(3) ; lo4 = bx%lo(4) ; lo5 = bx%lo(5)
   hi1 = bx%hi(1) ; hi2 = bx%hi(2) ; hi3 = bx%hi(3) ; hi4 = bx%hi(4) ; hi5 = bx%hi(5)
   mc = 1_I4P
   if (present(mask)) then
      mc = ubound(mask, dim=component_dim(bx))
      if (present(mask_comp)) mc = mask_comp
   endif
