   ! Locals of every box kernel of SELAH (selah_kernels_box): the loop range copied out of the box (the loops never
   ! reference the derived type, which keeps them device-friendly), the counters and the mask component.
   integer(I4P) :: lo1, lo2, lo3, lo4, lo5 !< First index of the loop range.
   integer(I4P) :: hi1, hi2, hi3, hi4, hi5 !< Last index of the loop range.
   integer(I4P) :: i1, i2, i3, i4, i5      !< Counters.
   integer(I4P) :: mc                      !< Mask component: mask_comp, or the last one (ADAM's all-solids summary).
