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
!<
!< **Fields of rank 1 to 4.** A field of lower rank is a rank-5 field with singleton dimensions: the kernels and the
!< stepper keep one rank-5 implementation. Dimension d of a rank-k field is dimension d of the rank-5 view when the
!< field's first dimension holds its components, d + 1 when the field is scalar (`box_dim`):
!<
!<| Field | As rank 5 `(v, i, j, k, b)` |
!<|-------|-----------------------------|
!<| `(i)`, scalar | `(1, i, 1, 1, 1)` |
!<| `(v, i)` | `(v, i, 1, 1, 1)` |
!<| `(v, i, j)` | `(v, i, j, 1, 1)` |
!<| `(v, i, j, k)` | `(v, i, j, k, 1)` |
!<| `(i, j)`, `(i, j, k)`, `(i, j, k, b)`, scalar | `(1, i, j, 1, 1)`, `(1, i, j, k, 1)`, `(1, i, j, k, b)` |
!<
!< The rank-k constructors (`box_interior(lb, ub, ngc, components)`, `box_range(lb, lo, hi, components)`) take the
!< `components` flag explicitly: the mapping is never guessed. Fields of rank below 5 have the layout
!< `LAYOUT_V_FIRST`.
use penf, only : I4P
implicit none
private
public :: box
public :: LAYOUT_V_FIRST, LAYOUT_V_LAST
public :: box_interior
public :: box_range
public :: box_components
public :: box_dim
public :: component_dim

integer(I4P), parameter :: LAYOUT_V_FIRST = 1_I4P !< Layout (v, i, j, k, b).
integer(I4P), parameter :: LAYOUT_V_LAST  = 2_I4P !< Layout (b, i, j, k, v).

type :: box
   !< The part of a rank-5 field a kernel updates.
   integer(I4P) :: layout = LAYOUT_V_FIRST !< Layout, LAYOUT_V_FIRST or LAYOUT_V_LAST.
   integer(I4P) :: lb(5)  = 1_I4P          !< Lower bounds of the allocation, ghosts included.
   integer(I4P) :: lo(5)  = 1_I4P          !< First index of the loop range, per dimension.
   integer(I4P) :: hi(5)  = 0_I4P          !< Last index of the loop range, per dimension.
   integer(I4P) :: rank   = 5_I4P          !< Rank of the field (1 to 5), viewed as rank 5.
   logical      :: components = .true.     !< The field's first dimension holds its components (always, for rank 5).
endtype box

interface box_interior
   !< The interior cells of a field: of the active blocks of a rank-5 field, or of a field of rank 1 to 4.
   module procedure box_interior_blocks, box_interior_rank
endinterface box_interior

contains
   pure function box_interior_blocks(layout, nv, ni, nj, nk, ngc, blocks_number) result(bx)
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
   endfunction box_interior_blocks

   pure function box_interior_rank(lb, ub, ngc, components) result(bx)
   !< The interior cells of a field of rank k = size(lb), allocated with the bounds lb:ub (ghosts included): ngc ghost
   !< cells off each end of every dimension but the components, all the components.
   integer(I4P), intent(in) :: lb(:)      !< Lower bounds of the field, ghosts included (size k, 1 to 5).
   integer(I4P), intent(in) :: ub(:)      !< Upper bounds of the field, ghosts included (size k).
   integer(I4P), intent(in) :: ngc        !< Ghost cells.
   logical,      intent(in) :: components !< The field's first dimension holds its components.
   type(box)                :: bx         !< Box.
   integer(I4P)             :: lo(size(lb)), hi(size(lb)) !< Loop range of the field.

   lo = lb + ngc
   hi = ub - ngc
   if (components) then
      lo(1) = lb(1)
      hi(1) = ub(1)
   endif
   bx = box_range(lb, lo, hi, components)
   endfunction box_interior_rank

   pure function box_range(lb, lo, hi, components) result(bx)
   !< The loop range lo:hi of a field of rank k = size(lb) whose lower bounds are lb, as a box.
   integer(I4P), intent(in) :: lb(:)      !< Lower bounds of the field, ghosts included (size k, 1 to 5).
   integer(I4P), intent(in) :: lo(:)      !< First index of the loop range (size k).
   integer(I4P), intent(in) :: hi(:)      !< Last index of the loop range (size k).
   logical,      intent(in) :: components !< The field's first dimension holds its components.
   type(box)                :: bx         !< Box.
   integer(I4P)             :: k, d       !< Rank, dimension.

   k = size(lb, kind=I4P)
   if (size(lo) /= k .or. size(hi) /= k) error stop 'selah box: lb, lo and hi must have the rank of the field'
   if (k < 1_I4P .or. k + merge(0_I4P, 1_I4P, components) > 5_I4P) &
      error stop 'selah box: a field has rank 1 to 5 (rank 5: with components)'
   bx%layout = LAYOUT_V_FIRST
   bx%rank = k
   bx%components = components
   bx%lb = 1_I4P ; bx%lo = 1_I4P ; bx%hi = 1_I4P
   do d=1_I4P, k
      bx%lb(box_dim(bx, d)) = lb(d)
      bx%lo(box_dim(bx, d)) = lo(d)
      bx%hi(box_dim(bx, d)) = hi(d)
   enddo
   endfunction box_range

   elemental function box_dim(bx, d) result(m)
   !< The dimension of the rank-5 view that holds dimension d of the field: d, or d + 1 for a scalar field.
   type(box),    intent(in) :: bx !< Box.
   integer(I4P), intent(in) :: d  !< Dimension of the field (1 to its rank).
   integer(I4P)             :: m  !< Dimension of the rank-5 view.

   m = d + merge(0_I4P, 1_I4P, bx%components)
   endfunction box_dim

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
