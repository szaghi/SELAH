# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.1.0] — 2026-11-17
### Added
- **tableaux**: The coefficient forms of a Runge-Kutta scheme and their conversion to Butcher

- **schemes**: The Runge-Kutta catalogue, verified scheme by scheme

- **splittings**: Symmetric composition schemes, Strang to Blanes-Moan S10

- **storage**: The buffers each storage form needs; use selah exports the whole of L0

- **kernels**: L1, the stage operations over a box or a flat array

- **stepper**: L2, a time step the caller drives one stage at a time


### Documentation
- **readme**: Authors, without Giacomo, in an ironic tone

- **guide**: The schemes and the L0 API, with every verification number of the night

- **guide**: The kernels, the L1 API, and status M2

- **guide**: Performance on the A30, the kernels on the device, status M3

- **guide**: The stepper, its contract and its verification; status M4

- **guide**: Whole steps through the stepper on the A30; status M5


### Fixed
- **ci**: Fetch the dependencies before the matrix build

- **tableaux**: Allocate the tableau before assigning it, for nvfortran

- **schemes**: No internal procedure defines a pure function's result

- **test**: Copy the coefficients into a variable before sending them to the device


### Performance
- **bench**: SELAH's kernels against ADAM's on the A30: parity per kernel, a 20% faster SSP(5,4) step

- **bench**: Whole steps through the stepper against ADAM's FNL step on the A30

- **kernels**: The stage state with no register to sum is a copy



