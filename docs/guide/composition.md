---
title: Composition
---

# Composing the stepper with your own operators

The stepper integrates the residual you give it. Two common needs fall outside one residual:
- a term you integrate yourself, exactly or with a method of your own: radiative cooling with Townsend's exact flow, a
  stiff relaxation, a particle step;
- particles whose kick and drift are exact sub-flows.

SELAH owns the composition in both cases: which operator runs, for what fraction of the step, from what time. You own
the operators. The theory, with its sources, is in `notes/limiters-splitting-particles.md` in the repository.

## Additive splitting: the field step and your flow

For $u' = F(u) + G(u)$, with $F$ the residual the stepper integrates and $G$ a term you integrate yourself:

```fortran
use selah
type(additive_splitting) :: sp
type(stepper)            :: st

call sp%initialize('strang')                       ! flow h/2, field h, flow h/2
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb, ub)
do n=1, nsteps
   do m=1, sp%substeps()
      h  = sp%substep_fraction(m) * dt
      tm = t + sp%substep_time(m) * dt
      if (sp%substep_operator(m) == OPERATOR_FIELD) then
         ! one step of the stepper of length h from time tm: begin_step, the stages (residual at
         ! tm + st%stage_time(k)*h), finish_step
      else
         ! your flow of G over h from time tm (OPERATOR_FLOW)
      endif
   enddo
   t = t + dt
enddo
```

| Scheme | Substeps | Order (exact sub-flows) | Composition, flow 1 first | Default `first` |
|---|---|---|---|---|
| `lie` (`SCHEME_LIE`) | 2 | 1 | field $h$, flow $h$ | `OPERATOR_FIELD` |
| `strang` | 3 | 2 | flow $h/2$, field $h$, flow $h/2$ | `OPERATOR_FLOW` |
| `yoshida-4`, `blanes-moan-s6`, `blanes-moan-s10` | 7, 13, 21 | 4, 4, 6 | the palindromes of [Schemes](./schemes) | `OPERATOR_FLOW`; `reversible=.true.` required |

**Which operator outside.** Put the cheap flow outside (the default): the field step, a Runge–Kutta step, runs once
per step, and an exact autonomous flow merges across steps ($\varphi_a \circ \varphi_b = \varphi_{a+b}$), so `strang`
then costs what `lie` does. With a stiff flow, outside is also the more accurate placement (Verwer & Sportisse 1998).

**Orders.** With both sub-flows exact the composition has its order; with the stepper as the field step, of order
$p \ge 2$, `strang` stays second order and `lie` first order. A stiff flow may reduce the order (to 1 for `strang` in
the stiff limit, in the literature); on the cooling model $u' = a u - k u^2$ with Townsend's form of exact flow,
`selah_test_additive_splitting` measures no reduction up to $k u_0 \Delta t = 125$, which is a measurement on that
model, not a guarantee.

**Backward substeps.** Every splitting of order 3 or more runs each operator backward over part of the step. That is
ill-posed for an irreversible operator: an upwind or limited field step, cooling, diffusion, accretion. `initialize`
stops the program on such a scheme unless you pass `reversible=.true.`, declaring both operators reversible;
`backward()` tells.

**Time.** Each operator integrates its own time over its substep: substep $m$ starts at `substep_time(m)`, the sum of
the same operator's fractions before $m$ (strang: the flow halves start at $t$ and $t + \Delta t/2$, the field step at
$t$). This is the composition of $(u, t_1, t_2)' = (F(u,t_1), 1, 0) + (G(u,t_2), 0, 1)$, and it keeps the order
whether $F$, $G$ or both depend on time (measured). It differs from the stepper's own splitting form, whose
`stage_time` freezes time through one partition: there the sub-flows are Euler steps of one partition, which cannot
carry time.

| Query | Returns |
|---|---|
| `initialize(scheme [, first] [, reversible])` | sets the composition; stops on an unknown scheme or on backward substeps without `reversible` |
| `substeps()` | the substeps of a step |
| `substep_operator(m)` | `OPERATOR_FIELD` (a step of the stepper) or `OPERATOR_FLOW` (your flow) |
| `substep_fraction(m)` | the substep's length as a fraction of the step (negative: backward) |
| `substep_time(m)` | the substep's start as a fraction of the step |
| `order()`, `backward()` | the order with exact sub-flows; whether a substep runs backward |

## Particles: kick and drift

A particle array whose first dimension holds positions and velocities, say `p(6, np)` with $x$ = 1:3 and $v$ = 4:6,
is a rank-2 field for the stepper's splitting form. The kick $v' = a(x)$ and the drift $x' = v$ are exact sub-flows
(each partition's residual does not depend on its own variables), so the splittings reach their orders, and `strang`
with the velocities as partition 1 is kick–drift–kick leapfrog:

```fortran
call st%initialize(SCHEME_STRANG, FORM_SPLITTING, lb=[1, 1], ub=[6, np], components=.true., &
                   partitions=reshape([4, 6, 1, 3], [2, 2]))   ! partition 1: v (kick), 2: x (drift)
bx = box_interior(lb=[1, 1], ub=[6, np], ngc=0, components=.true.)
call st%begin_step(bx, p)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, p)
   call st%rhs_register_view(r)
   if (st%stage_partition(k) == 1) then
      ! r(4:6,:) = the acceleration at the positions p(1:3,:)
   else
      r(1:3,:) = p(4:6,:)
   endif
   call st%end_stage(k, bx, dt, p)
enddo
call st%finish_step(bx, dt, p)
```

The half kick is SELAH's fraction, not a typed 1/2: `stage_euler_fraction(1)` is 0.5. The higher-order compositions
(yoshida-4, blanes-moan-s6, -s10) run the kick and the drift backward over part of the step, which is harmless for a
particle Hamiltonian. With a fixed step the energy error stays bounded (no drift: `selah_test_particle_flows`, 50
Kepler orbits); a step that changes from step to step loses that property.

**Particles and a field.** Gas and sink particles, for instance: the gas is the field step, one step of the particle
stepper over the substep is the flow, and `additive_splitting` composes them (`lie` first order, `strang` second
order: measured). Accretion is irreversible: `lie` or `strang` only.

## Verification

| What | Case |
|---|---|
| lie, strang and the palindromes reach their orders with exact sub-flows, with the stepper as the field step, and with both terms time-dependent; the cooling model | [A6](/verification/accuracy/additive-splitting) |
| KDK's discrete solution on the oscillator, the orders on Kepler, bounded energy; gas and sinks coupled | [A7](/verification/accuracy/particle-flows) |
