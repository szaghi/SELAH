! the mask has the field's bounds and one component; a cell is active (fluid) where mask < 0
allocate(mask(1:1, 0:NI+1, 0:2, 0:2, 1:NB))
mask = -1._R_P
mask(:,:,:,:,3) = 1._R_P                                          ! block 3 is solid
