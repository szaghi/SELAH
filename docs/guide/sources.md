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

## Limiters, additive splitting and particles

- X. Zhang, C.-W. Shu, *On positivity-preserving high order discontinuous Galerkin schemes for compressible Euler
  equations on rectangular meshes*, J. Comput. Phys. 229 (2010) 8918–8934. A limiter inside SSP stages: the
  forward-Euler step each stage must make admissible (`stage_euler_fraction`).
- B. Sportisse, *An analysis of operator splitting techniques in the stiff case*, J. Comput. Phys. 161 (2000)
  140–168. Strang splitting with a stiff operator: its order reduction, the stiff operator outside.
- S. Blanes, F. Casas, *On the necessity of negative coefficients for operator splitting schemes of order higher than
  two*, Appl. Numer. Math. 54 (2005) 23–37. Why `additive_splitting` refuses the high-order compositions for an
  irreversible operator.
- R. H. D. Townsend, *An exact integration scheme for radiative cooling in hydrodynamical simulations*, Astrophys. J.
  Suppl. 181 (2009) 391–397. The exact cooling flow of the model problem of [A6](/verification/accuracy/additive-splitting).
- T. Quinn, R. P. Perrine, D. C. Richardson, R. Barnes, *A symplectic integrator for Hill's equations*, Astron. J. 139
  (2010) 803–807. Kick–drift–kick for sink particles.
- J. P. Boris, *Relativistic plasma simulation: optimization of a hybrid code*, Proc. Fourth Conf. Numerical
  Simulation of Plasmas, NRL, 1970, 3–67. `PUSHER_BORIS`.
- H. Qin, S. Zhang, J. Xiao, J. Liu, Y. Sun, W. M. Tang, *Why is Boris algorithm so good?*, Phys. Plasmas 20 (2013)
  084503. Volume preservation.
- J.-L. Vay, *Simulation of beams or plasmas crossing at relativistic velocity*, Phys. Plasmas 15 (2008) 056701.
  `PUSHER_VAY`, the force-free test.
- A. V. Higuera, J. R. Cary, *Structure-preserving second-order integration of relativistic charged particle
  trajectories in electromagnetic fields*, Phys. Plasmas 24 (2017) 052104. `PUSHER_HIGUERA_CARY`.
- T. Zh. Esirkepov, *Exact charge conservation scheme for Particle-in-Cell simulation with an arbitrary form-factor*,
  Comput. Phys. Commun. 135 (2001) 144–153. The charge-conserving current the stage currents rest on.

## Perfectly matched layers and exponential integrators

- J.-P. Bérenger, *A perfectly matched layer for the absorption of electromagnetic waves*, J. Comput. Phys. 114 (1994)
  185–200.
- M. Kuzuoglu, R. Mittra, *Frequency dependence of the constitutive parameters of causal perfectly matched
  anisotropic absorbers*, IEEE Microw. Guided Wave Lett. 6 (1996) 447–449. The complex-frequency-shifted PML.
- J. A. Roden, S. D. Gedney, *Convolution PML (CPML): an efficient FDTD implementation of the CFS-PML for arbitrary
  media*, Microw. Opt. Technol. Lett. 27 (2000) 334–339. The recursion whose Runge–Kutta analogue `decay=` is.
- M. Hochbruck, A. Ostermann, *Exponential integrators*, Acta Numer. 19 (2010) 209–286. The exponential stage format
  and its equilibrium conditions.

## Software

- [ADAM](https://github.com/szaghi/adam): the code SELAH was carved from, and the oracle of [C1](/verification/contracts/adam-equivalence).
- [FUNDAL](https://github.com/szaghi/FUNDAL): device memory and the portable directives.
- [PENF](https://github.com/szaghi/PENF): the kinds.
- [SETH](https://github.com/szaghi/SETH): the documentation and verification tooling SELAH's follows.
