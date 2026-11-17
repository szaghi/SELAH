# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.2.0] — 2026-11-17
### Added
- **stepper**: Stages on q, an opt-in that spares the Butcher form its stage-1 copy


### Documentation
- **guide**: Stages on q, the option, the query and their measurement

- The theme, figure tooling and navigation of the documentation refactor

- **tutorial**: Ten chapters, each a real program whose code and output the pages include

- **cookbook**: Recipes with real code and output, porting and performance references

- **figures**: The verification figures, tables and casts, drawn from the report data

- **verification**: The verification section, case by case, with its catalogue and regression page

- **guide**: The guide rewritten from the audit: concepts, every procedure's preconditions, the home page

- **verification**: The final A30 run of the documented tree

- README in SETH's layout, ADAM out of the landing pages

- New identity, Hamon in Tsuyukusa (Sumi family)


### Fixed
- **stepper,kernels**: The defects of the documentation audit, and the tests it found missing


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



