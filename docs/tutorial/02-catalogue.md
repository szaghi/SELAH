---
title: 2. The catalogue
---

# 2. The catalogue

SELAH ships a fixed catalogue: eight Runge–Kutta schemes and four splittings. Each comes with its coefficients, its
order, its SSP coefficient and the storage forms it can run in. This chapter prints them.

<<< @/examples/snippets/catalogue.F90{fortran}

- `RK_SCHEMES` and `SPLITTING_SCHEMES` list the names. They are fixed-length strings padded with blanks, so `trim` them
  before passing them on.
- `registers_needed(scheme, form)` is the number of field-sized buffers a form needs, or −1 when the scheme has no
  such form. The Butcher form needs $s+1$: the stage state and one residual per stage.
- `butcher_from_shu_osher_2r` converts the two-register coefficients back to a Butcher tableau. The difference is
  exactly zero: the conversion uses the summation order of ADAM's code.
- `splitting(name)` returns the flows of a composition: which partition each flow advances (`partition`) and by what
  fraction of the step (`coefficient`). Both partitions add up to one whole step.

## Running it

<<< @/examples/output/catalogue.txt{text}

Read the table by rows. `ssp-104` needs 11 registers in Butcher form but 2 in Ketcheson's form. `rk4` and the 2N
schemes have no SSP coefficient: 0 means "not SSP" (or "not computed", for a tableau you build yourself).

::: warning The SSP coefficient is a lower bound
`ssp_coefficient` is computed from the Shu–Osher coefficients SELAH stores. A different Shu–Osher representation of
the same scheme can give a larger value. For the catalogue the values match the published optimal ones
([P2](/verification/properties/ssp-coefficients)).
:::

::: tip What you learned
The catalogue, the storage forms and their register counts, and the conversions between coefficient forms.
Reference: [Schemes](/guide/schemes), [API: L0](/api/#l0-coefficients).
:::

Next: [3. A field and a box](./03-box).
