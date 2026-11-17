!< SELAH, L1: the box, the part of a block-structured field a kernel updates.
module selah_box
!< SELAH, L1: the box, the part of a block-structured field a kernel updates.
!<
!< A field is a rank-5 array: components, three cell indices with ghost cells, and blocks, in one of two layouts:
!<
!<+ `LAYOUT_V_FIRST`: `(v, i, j, k, b)`, the components contiguous (ADAM's CPU fields);
!<+ `LAYOUT_V_LAST`: `(b, i, j, k, v)`, the blocks contiguous (ADAM's FNL device fields).
!<
!< A box holds the lower bounds of the allocation (`lb`, ghosts included) and the range the kernels loop over (`lo:hi`):
!< the interior cells of the active blocks `1..blocks_number` of a pool that may hold more, and all components or a
!< range of them (a partition of a splitting scheme). The kernels never write outside the box: not the ghost cells, not
!< the pool slots past `blocks_number`. A box is passed at every call, never stored, because the active blocks change at
!< every regrid.
use penf, only : I4P
implicit none
private
public :: box
public :: LAYOUT_V_FIRST, LAYOUT_V_LAST
public :: box_interior
public :: box_components
public :: component_dim

integer(I4P), parameter :: LAYOUT_V_FIRST = 1_I4P !< Layout (v, i, j, k, b).
integer(I4P), parameter :: LAYOUT_V_LAST  = 2_I4P !< Layout (b, i, j, k, v).

type :: box
   !< The part of a rank-5 field a kernel updates.
   integer(I4P) :: layout = LAYOUT_V_FIRST !< Layout, LAYOUT_V_FIRST or LAYOUT_V_LAST.
   integer(I4P) :: lb(5)  = 1_I4P          !< Lower bounds of the allocation, ghosts included.
   integer(I4P) :: lo(5)  = 1_I4P          !< First index of the loop range, per dimension.
   integer(I4P) :: hi(5)  = 0_I4P          !< Last index of the loop range, per dimension.
endtype box

contains
   pure function box_interior(layout, nv, ni, nj, nk, ngc, blocks_number) result(bx)
   !< The interior cells of the active blocks, all components; the arrays are allocated with ghosts 1-ngc..n+ngc on the
   !< cell indices and from 1 on the component and block indices.
   integer(I4P), intent(in) :: layout        !< Layout.
   integer(I4P), intent(in) :: nv            !< Components.
   integer(I4P), intent(in) :: ni, nj, nk    !< Interior cells per direction.
   integer(I4P), intent(in) :: ngc           !< Ghost cells.
   integer(I4P), intent(in) :: blocks_number !< Active blocks.
   type(box)                :: bx            !< Box.

   bx%layout = layout
   bx%lb = [1_I4P, 1_I4P-ngc, 1_I4P-ngc, 1_I4P-ngc, 1_I4P]
   bx%lo = 1_I4P
   if (layout == LAYOUT_V_LAST) then
      bx%hi = [blocks_number, ni, nj, nk, nv]
   else
      bx%hi = [nv, ni, nj, nk, blocks_number]
   endif
   endfunction box_interior

   pure function component_dim(bx) result(d)
   !< The dimension of the components: 1 (V_FIRST) or 5 (V_LAST).
   type(box), intent(in) :: bx !< Box.
   integer(I4P)          :: d  !< Dimension.

   d = merge(5_I4P, 1_I4P, bx%layout == LAYOUT_V_LAST)
   endfunction component_dim

   pure function box_components(bx, v1, v2) result(sub)
   !< The same box restricted to the components v1..v2 (the variables of one partition).
   type(box),    intent(in) :: bx     !< Box.
   integer(I4P), intent(in) :: v1, v2 !< Component range.
   type(box)                :: sub    !< Restricted box.

   sub = bx
   sub%lo(component_dim(bx)) = v1
   sub%hi(component_dim(bx)) = v2
   endfunction box_components
endmodule selah_box
