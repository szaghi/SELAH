---
title: 3. A field and a box
---

# 3. A field and a box

A block-structured solver stores a field as a pool of blocks. Each block has interior cells and a frame of ghost cells,
and each cell has `NV` components. SELAH's kernels take such an array in one of two layouts, and a **box**: the range
of cells they loop over. They write those cells, and nothing else.

<<< @/examples/snippets/box.F90{fortran}

- **The array.** Ghost cells sit on the cell indices, `1-NGC : N+NGC`; components and blocks start at 1. The pool
  (`NB = 4` blocks) may be larger than the active set.
- **`box_interior(layout, nv, ni, nj, nk, ngc, blocks)`** builds the box of the interior cells of the first `blocks`
  blocks, every component:
  - `lb` holds the array's lower bounds, which the kernels need to address it;
  - `lo:hi` is the loop range.
- **`box_components(bx, v1, v2)`** narrows a box to the components `v1..v2`: a partition of a splitting scheme, or the
  passive scalars of a solver.
- **The layouts.**
  - `LAYOUT_V_FIRST` is `(v,i,j,k,b)`, components innermost: the CPU layout.
  - `LAYOUT_V_LAST` is `(b,i,j,k,v)`, blocks innermost: ADAM's device layout, where consecutive device threads read
    consecutive blocks.
  - `component_dim(bx)` says which dimension holds the components.

## Running it

<<< @/examples/output/box.txt{text}

The kernel updated exactly the 144 interior values of the 3 active blocks. The ghost cells and the fourth block of the
pool kept the sentinel bit for bit, and so did component 1 when the box was narrowed to component 2. This is the
first clause of SELAH's stage contract: **a kernel never writes a ghost cell**. Ghost cells belong to the caller (see
[Concepts](/guide/concepts)).

::: warning Explicit bounds
A field with ghost cells must reach the kernels with its own lower bounds. An assumed-shape dummy without bounds
silently renumbers it from 1, and the interior is then off by `NGC`. Every SELAH dummy is declared
`q(bx%lb(1):, ..., bx%lb(5):)`.
:::

::: tip What you learned
The field, the box (`lb`, `lo`, `hi`), the two layouts, component ranges; a kernel writes the box and nothing else.
Reference: [Concepts](/guide/concepts), [Kernels](/guide/kernels).
:::

Next: [4. The stepper](./04-stepper).
