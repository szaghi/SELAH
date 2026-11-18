---
title: Any rank
---

# Fields of any rank

[[toc]]

SELAH's stepper takes a field of rank 1 to 5 at its own rank. A field of rank below 5 is viewed as rank 5 with
singleton dimensions, so the kernels and the stage contract are the same at every rank, and nothing is copied. The
mapping and the preconditions: [Kernels](/guide/kernels#fields-of-rank-1-to-4),
[Stepper](/guide/stepper#fields-of-rank-1-to-4).

## A 1-D, 2-D or 3-D solver

Build the box and the stepper from the field's own bounds. Say whether its first dimension holds components:

| Field | Box | Stepper |
|---|---|---|
| `q(1-g:n+g)`, scalar | `box_interior([1-g], [n+g], g, components=.false.)` | `initialize(..., lb=lbound(q), ub=ubound(q), components=.false.)` |
| `q(nv, 1-g:n+g)`, components | `box_interior([1, 1-g], [nv, n+g], g, components=.true.)` | `components=.true.` |
| `q(1-g:n+g, 1-g:m+g)`, scalar 2-D | `box_interior([1-g, 1-g], [n+g, m+g], g, components=.false.)` | `components=.false.` |
| `q(nv, 1-g:n+g, 1-g:m+g)`, components 2-D | `box_interior([1, 1-g, 1-g], [nv, n+g, m+g], g, components=.true.)` | `components=.true.` |
| `q(nv, i, j, k)`, components 3-D | `box_interior(lbound(q), ubound(q), g, components=.true.)` | `components=.true.` |

The four step procedures take `q` at its rank. The stage state and the residual register are
`call st%stage_state_view(u)` and `call st%rhs_register_view(r)`, pointers of the field's rank with its bounds, as in
the [tutorial](/tutorial/04-stepper#a-field-of-any-rank):

<<< @/examples/snippets/ranks-step.F90{fortran}

## A periodic field with no ghost cells

A solver that wraps its indices (`u(0:n-1)`, the neighbours of cell `i` taken modulo `n` by the residual) has no ghost
cells. Its box is the whole field, with `ngc = 0`:

```fortran
bx = box_interior(lb=[0], ub=[n-1], ngc=0, components=.false.)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[0], ub=[n-1], components=.false.)
```

## What to watch

- **Contiguous fields.** The field reaches the rank-5 procedures without a copy only when it is contiguous. A section
  such as `q(1, :)` of a larger array makes the compiler copy it in and out at every call. In a device build that
  copy would run on the host over device memory, so pass whole arrays.
- **The rank is checked at every call.** The box, the stepper and `q` must agree in rank, in `components` and in
  bounds; a mismatch stops with a message.
- **Caller registers** (`registers=`) are for rank-5 fields only. At rank 1 to 4 the stepper allocates its own.
