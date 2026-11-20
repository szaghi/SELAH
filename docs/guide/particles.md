---
title: Particles
---

# Particles

Charged particles in a field solver (particle-in-cell) are advanced in one of two ways, and SELAH provides both:
- **with the field's Runge–Kutta scheme**, in the same stage loop: the particles are a
  [companion](./stepper#companion-states) of the field stepper, their residual the Lorentz force at each stage;
- **with a two-level pusher** (Boris and its relativistic successors), staggered in time with the fields: a kernel,
  `selah_push`.

The theory, with its sources: `notes/companions-pushers-pml.md` in the repository.

## With the field's scheme

The particles' residual is $(v, (q/m)(E + v\times B))$ with $E$ and $B$ gathered at the stage positions, and the
field's residual gets the current deposited from the stage particles. Each stage of the coupled step is then the
scheme's. What this gives:

| Property | Value |
|---|---|
| stability of the gyration $\omega_c$ (and of every oscillatory mode $\omega$) | $\omega\,\Delta t \le$ `st%imag_axis_bound()`: 0 for euler and ssp-22 (unstable at any step), $\sqrt 3$ for ssp-33, $2\sqrt 2$ for rk4, 3.2784 for ssp-54, 4.9215 for ssp-104 |
| gyration energy | decays by $\lvert R(i\omega_c\Delta t)\rvert^2$ per step, every stable scheme (rk4: about $(\pi/36)(\omega_c\Delta t)^5$ per gyration) |
| phase-space volume | not preserved |
| $E\times B$ drift | exact: $v_E$ is an equilibrium of every explicit scheme |

## Two-level pushers

```fortran
! p(8, np): x at 1:3, v (or u = gamma v) at 4:6, q at 7, m at 8; e, b: the fields at x^n, (3, np)
call selah_push(PUSHER_HIGUERA_CARY, LAYOUT_V_FIRST, n, dt, 1, 4, 7, 8, p, e, b, c=c_light)
```

The positions are at $t^n$, the velocities (momenta) at $t^{n-1/2}$; one call takes the first `n` particles to
$v^{n+1/2}$, $x^{n+1}$. `LAYOUT_V_LAST` takes `p(np, 8)`, `e(np, 3)`, `b(np, 3)`. The relativistic pushers hold
$u = \gamma v$ and need $c$; a Gaussian-units caller passes $B/c$.

| Pusher | Use | Volume | $E + v\times B = 0$ | $E\times B$ drift |
|---|---|---|---|---|
| `PUSHER_BORIS` | non-relativistic | preserved | — | exact |
| `PUSHER_BORIS_RELATIVISTIC` | relativistic, the classic | preserved | **not** kept (1.2e-7 per step at $\gamma = 10$, $qB\Delta t/m = 0.1$) | wrong velocity (2 % at $E/cB = 0.5$, $qB\Delta t/m = 1$) |
| `PUSHER_VAY` | relativistic beams | **not** preserved ($\det - 1 = -1.2\cdot10^{-4}$ at $q\Delta t/m = 1$) | kept | exact |
| `PUSHER_HIGUERA_CARY` | relativistic, the structure-preserving choice | preserved | kept | exact |

Boris turns the velocity by $2\arctan(\omega_c\Delta t/2)$ per step (the exact angle is $\omega_c\Delta t$), keeps
$|v|$ exactly in a pure magnetic field, and is stable at any step (aliased above $\omega_c\Delta t \approx 2$).

## Charge-conserving stage currents

A particle-in-cell code keeps Gauss's law exactly when its current is the one that moves the deposited charge (a
charge-conserving deposit, Esirkepov 2001). With the particles in the field's Runge–Kutta step, the field at stage
$k+1$ sees $\sum_{r\le k} a_{k+1,r} J_r$: Gauss's law holds at every stage when that weighted sum is the
charge-conserving current of the displacement from $q^n$ to the particles' state at stage $k+1$. Within the window of
stage $k$, after the particles' residual:

```fortran
call sp%lookahead(k, bp, dt, p, pt)            ! the particles' state at stage k+1 (the step's result at the last stage)
! ... deposit H, the charge-conserving current of the displacement p -> pt over dt (the caller's: Esirkepov)
call st%stage_current(k, bx, h, jall, jk)      ! J_k = (H - sum_{r<k} a_{k+1,r} J_r) / a_{k+1,k}
! ... jk into jall(..., k); the field's residual at stage k with the current jk
```

Then $\sum_{r\le k} a_{k+1,r} J_r = H$ at every stage and $\sum_k b_k J_k = H(q^n \to q^{n+1})$ at the step. It
needs $a_{k+1,k} \ne 0$ and $b_s \ne 0$, true for the whole catalogue (`stage_current_admissible()`). The individual
$J_k$ are not physical currents: the inversion amplifies them (up to 25 times for ssp-104); only their weighted sums
are. This is PRISM's construction; no paper on it was found (the theory note gives its proof). In two and three
dimensions Esirkepov's current follows the straight chord between the end points, so the coupling is second order in
time whatever the scheme's order.

## Verification

| What | Case |
|---|---|
| companions, the gyration on a companion, the axis bounds, the drift | [C8](/verification/contracts/companion), [A5](/verification/accuracy/stability) |
| the four pushers: defining equations, angle, speed, drift, force-free, volume, the non-relativistic limit | [A8](/verification/accuracy/pushers) |
| Gauss's law at every stage and step with the stage currents; the look-ahead state | [C9](/verification/contracts/stage-currents) |
