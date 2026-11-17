---
title: Sources
---

# Sources

The schemes, the forms and the properties SELAH implements and verifies, with the work they come from.

## Runge–Kutta schemes and forms

- J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed., Wiley, 2016. The Butcher form,
  the order conditions as rooted trees.
- W. Kutta, *Beitrag zur näherungsweisen Integration totaler Differentialgleichungen*, Z. Math. Phys. 46 (1901)
  435–453. The classic `rk4`.
- C.-W. Shu, S. Osher, *Efficient implementation of essentially non-oscillatory shock-capturing schemes*, J. Comput.
  Phys. 77 (1988) 439–471. The Shu–Osher form, `ssp-22` and `ssp-33`.
- S. Gottlieb, C.-W. Shu, *Total variation diminishing Runge–Kutta schemes*, Math. Comp. 67 (1998) 73–85.
- S. Gottlieb, C.-W. Shu, E. Tadmor, *Strong stability-preserving high-order time discretization methods*, SIAM Rev.
  43 (2001) 89–112. The SSP property and coefficient, the theorem of [P1](/verification/properties/ssp-tvd).
- R. J. Spiteri, S. J. Ruuth, *A new class of optimal high-order strong-stability-preserving time discretization
  methods*, SIAM J. Numer. Anal. 40 (2002) 469–491. `ssp-54`.
- L. Isherwood, Z. J. Grant, S. Gottlieb, *Strong stability preserving integrating factor Runge–Kutta methods*, SIAM J.
  Numer. Anal. 56 (2018) 3276–3307. The 15-digit coefficients of `ssp-54`.
- D. I. Ketcheson, *Highly efficient strong stability-preserving Runge–Kutta methods with low-storage
  implementations*, SIAM J. Sci. Comput. 30 (2008) 2113–2136. `ssp-104` and its two-register recipe.
- J. H. Williamson, *Low-storage Runge–Kutta schemes*, J. Comput. Phys. 35 (1980) 48–56. The 2N form, `lsrk-33`.
- M. H. Carpenter, C. A. Kennedy, *Fourth-order 2N-storage Runge–Kutta schemes*, NASA TM-109112, 1994. `lsrk-54`.

## Splitting schemes

- G. Strang, *On the construction and comparison of difference schemes*, SIAM J. Numer. Anal. 5 (1968) 506–517.
- H. Yoshida, *Construction of higher order symplectic integrators*, Phys. Lett. A 150 (1990) 262–268. `yoshida-4`.
- S. Blanes, P. C. Moan, *Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods*, J. Comput.
  Appl. Math. 142 (2002) 313–330. `blanes-moan-s6`, `blanes-moan-s10`.
- E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer, 2006. Symplecticity and the
  long-time energy behaviour of [P3](/verification/properties/splitting-energy).

## Software

- [ADAM](https://github.com/szaghi/adam): the code SELAH was carved from, and the oracle of [C1](/verification/contracts/adam-equivalence).
- [FUNDAL](https://github.com/szaghi/FUNDAL): device memory and the portable directives.
- [PENF](https://github.com/szaghi/PENF): the kinds.
- [SETH](https://github.com/szaghi/SETH): the documentation and verification tooling SELAH's follows.
