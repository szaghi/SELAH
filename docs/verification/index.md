---
title: Verification
---

# Verification

Every claim these pages make about SELAH is backed by a verification case: a Fortran program in `src/tests/` that runs
the library and asserts the result against something that shares no code with it. The references are:
- an exact solution;
- a closed form;
- a theorem;
- a published coefficient;
- a golden file;
- ADAM's own code, copied byte for byte.

The same programs are the test suite. `bash scripts/run_tests.sh` runs all of them, with gfortran, nvfortran and ifx,
in double and single precision, and on an NVIDIA A30. With `SELAH_VV_DIR` set, the same programs also write the data
these pages are drawn from (`scripts/run_vv.sh`). The figures and tables are SELAH's real results, regenerated after
any change.

Every case page has the same sections: case and scope, setup, the program and how to run it, results, findings. A
card's tag says what kind of check it is:

- **physics**: an order of accuracy, a stability region, a strong-stability or conservation property, against a
  closed form or a theorem;
- **contract**: a bitwise identity, an equivalence (with ADAM, between storage forms), the stage contract, a pinned
  regression value, host–device agreement, or the suite's own power to detect defects.

<div class="selah-casts">
<figure class="selah-cast" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-light.png" aria-label="SSP(5,4) advecting a square wave at the SSP limit and above it"><source src="/verification/casts/ssp-tvd-light.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-dark.png" aria-label="SSP(5,4) advecting a square wave at the SSP limit and above it"><source src="/verification/casts/ssp-tvd-dark.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-dark.mp4" type="video/mp4"></video>
<figcaption>SSP(5,4) with first-order upwind at ν = C: the total variation never grows. Above its threshold it does. <a href="./properties/ssp-tvd">Case P1</a>.</figcaption>
</figure>
<figure class="selah-cast" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-light.png" aria-label="Phase portrait of the oscillator: Strang splitting against rk4"><source src="/verification/casts/splitting-phase-light.webm" type="video/webm"><source src="/verification/casts/splitting-phase-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-dark.png" aria-label="Phase portrait of the oscillator: Strang splitting against rk4"><source src="/verification/casts/splitting-phase-dark.webm" type="video/webm"><source src="/verification/casts/splitting-phase-dark.mp4" type="video/mp4"></video>
<figcaption>The oscillator at Δt = 0.6. The Strang splitting stays on a closed curve; rk4 spirals inwards. <a href="./properties/splitting-energy">Case P3</a>.</figcaption>
</figure>
</div>

## What a user of SELAH should know

The findings that change how SELAH should be used, each documented on its page:

- **The SSP coefficient is a sufficient bound.** Every SSP scheme keeps a forward-Euler TVD operator TVD up to
  $\nu = C\,\nu_{FE}$ ([P1](./properties/ssp-tvd)). On the linear upwind operator the actual threshold is higher
  (ssp-54: 1.25 C), but only $C$ is guaranteed for every TVD operator.
- **SSP(10,4) is the most economical SSP scheme of the catalogue.** Its effective coefficient $C/s$ is 0.6, against
  0.30 for ssp-54, and it runs in two registers ([P2](./properties/ssp-coefficients)).
- **`stage_weight(k)` is the only correct weight for a flux register.** Uniform weights violate conservation at the
  size of the time error for every scheme whose weights are not uniform ([P4](./properties/conservation)).
- **A splitting reaches its order only with exact sub-flows.** The stepper's sub-flow is forward Euler on one
  partition: exact when that partition's residual does not depend on itself (separable Hamiltonians), first order
  otherwise ([A4](./accuracy/splittings)).
- **In single precision some schemes lose their exact stability properties at round-off.** rk4, for instance, keeps
  its imaginary-axis interval only to within about 1e-8 of $|R| = 1$ ([A5](./accuracy/stability),
  [C7](./contracts/precision)).

## Accuracy

<SelahIcon name="order-conditions" :size="40" label="accuracy" float />

Order conditions, observed orders on scalar problems, through the stepper and for the splittings, and the linear
stability regions. [Group page](./accuracy/).

<div class="selah-gallery selah-catalog">

[![Order conditions](/verification/order-conditions.png)<strong>A1 Order conditions</strong><span>the eight conditions up to order 4, in quadruple precision, for every scheme</span><em>physics</em>](./accuracy/order-conditions)

[![Scalar orders](/verification/scalar-orders-2.png)<strong>A2 Orders on scalar ODEs</strong><span>every scheme on a linear and a nonlinear ODE, four halvings</span><em>physics</em>](./accuracy/scalar-orders)

[![Stepper orders](/verification/stepper-orders.png)<strong>A3 Orders through the stepper</strong><span>every scheme in every storage form, on a multi-block advection with an exact semi-discrete solution</span><em>physics</em>](./accuracy/stepper-orders)

[![Splitting orders](/verification/splitting-orders.png)<strong>A4 Splittings</strong><span>the four compositions on the oscillator and the pendulum, both partition orders</span><em>physics</em>](./accuracy/splittings)

[![Stability](/verification/stability.png)<strong>A5 Linear stability</strong><span>the regions |R(z)| ≤ 1 from the library's tableaux, against closed forms</span><em>physics</em>](./accuracy/stability)

[<strong>A6 Additive splitting</strong><span>the field step and the caller's flow: orders with exact sub-flows, through the stepper, time-dependent, stiff</span><em>physics</em>](./accuracy/additive-splitting)

[<strong>A7 Particle sub-flows</strong><span>KDK's discrete solution, Kepler orders and bounded energy, gas and sinks coupled</span><em>physics</em>](./accuracy/particle-flows)

</div>

## Properties

<SelahIcon name="ssp" :size="40" label="properties" float />

Strong stability, the SSP coefficients, long-time energy behaviour and conservation. [Group page](./properties/).

<div class="selah-gallery selah-catalog">

[![SSP TVD](/verification/ssp-tvd.png)<strong>P1 Strong stability (TVD)</strong><span>a CFL sweep of a TVD operator: no growth of the total variation up to C</span><em>physics</em>](./properties/ssp-tvd)

[<strong>P2 SSP coefficients</strong><span>computed from the Shu–Osher forms, against the published optimal values</span><em>physics</em>](./properties/ssp-coefficients)

[![Splitting energy](/verification/splitting-energy.png)<strong>P3 Splitting energy</strong><span>1000 periods of the oscillator: bounded energy error for splittings, drift for Runge–Kutta</span><em>physics</em>](./properties/splitting-energy)

[![Telescoping](/verification/telescoping.png)<strong>P4 Conservation</strong><span>the stage fluxes weighted by stage_weight telescope to the step's mass change</span><em>physics</em>](./properties/conservation)

</div>

## Contracts

<SelahIcon name="selah" :size="40" label="contracts" float />

Equivalences, the stage contract, pinned values, the device and the suite's own sensitivity. [Group page](./contracts/).

<div class="selah-gallery selah-catalog">

[<strong>C1 Equivalence with ADAM</strong><span>bitwise ADAM's FNL kernels and coefficients, on the host and on the A30</span><em>contract</em>](./contracts/adam-equivalence)

[![Forms](/verification/forms.png)<strong>C2 Storage forms</strong><span>every low-storage form computes its Butcher form's solution to round-off</span><em>contract</em>](./contracts/forms)

[<strong>C3 Stage contract</strong><span>ghost cells, aliasing, stages on q, masks, caller registers, the error stops</span><em>contract</em>](./contracts/stage-contract)

[<strong>C4 Golden regression</strong><span>1312 values pinned in a file written once, compared on every compiler</span><em>contract</em>](./contracts/golden)

[<strong>C5 Device</strong><span>the kernels and the stepper on an NVIDIA A30, against the host and ADAM</span><em>contract</em>](./contracts/device)

[<strong>C6 Mutation testing</strong><span>28 deliberate defects: does the suite catch every one?</span><em>contract</em>](./contracts/mutation)

[<strong>C7 Precision</strong><span>what runs and what is asserted in single precision</span><em>contract</em>](./contracts/precision)

</div>

## Running it

```bash
fobis fetch                                   # PENF, FUNDAL and GIDEON
fobis build --mode tests-gnu && bash scripts/run_tests.sh
bash scripts/run_vv.sh                        # report mode: data, figures, tables, casts
python3 scripts/mutation_test.py              # C6: about half an hour
```

The list of programs, with what each asserts, is on [Test programs](./unit-tests). What is pinned between releases,
and the rules for changing it, are on [Regression](./regression).
