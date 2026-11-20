---
title: P5 Exact decay
icon: ssp
---

# P5 Exact decay

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P5 | Properties | P | CPU (gfortran, nvfortran, ifx; R8P and R4P, the orders in R8P); GPU | no | `selah_test_decay` |

## Case and scope

`decay=lambda` integrates $\psi' = g - \lambda\psi$'s decay exactly, the caller's residual the source alone
([guide](/guide/stepper#a-decay-integrated-exactly)). The case checks the properties the theory gives the variant
(`notes/companions-pushers-pml.md` §5): exactness for a frozen source at any $\lambda\Delta t$, the explicit step at
$\lambda = 0$, exact static equilibria and stability of a coupled field–auxiliary pair at a stiff $\lambda\Delta t$,
and its order.

## Setup

| Item | Value |
|---|---|
| Frozen source | 5 cells, $\lambda\Delta t$ = 0.1, 1, 10, $10^3$, $10^6$, $g = 1$, $\psi^0 = 3$, 8 steps; exact $g/\lambda + (\psi^0 - g/\lambda)e^{-\lambda t}$ |
| $\lambda = 0$ | 3 steps against the same steps without `decay` |
| Coupled pair | $u' = -\psi + J$ (the field stepper), $\psi' = \sigma u - \lambda\psi$ (a companion with `decay`), $\sigma = 3$, $J = 1/2$ |
| Equilibrium | from $u^* = \lambda J/\sigma$, $\psi^* = J$, $\lambda\Delta t = 50$, 4 steps |
| Stability | from $(1, 0)$, $\lambda\Delta t = 50$, 20 steps: with `decay`, and with $-\lambda\psi$ in the residual (explicit) |
| Order | $\lambda = 2$, from $(1, 0)$ to $t = 1$, 16 to 128 steps; reference rk4, $10^5$ steps |

## Program

```bash
./exe/selah_test_decay
```

## Results

| Check | Result | Bound |
|---|---|---|
| frozen source, the largest relative error over every $\lambda\Delta t$ and step | 3.4e-16 … 1.0e-15 (ssp-54), every scheme | $10^{-13}$ |
| $\lambda = 0$ | bitwise the explicit step, every scheme | exact |
| coupled equilibrium at $\lambda\Delta t = 50$ | kept, every scheme | 64 ε |
| stability at $\lambda\Delta t = 50$ | bounded with `decay`; the explicit stages blow up, every scheme | — |
| coupled order | 1.01 (euler), 1.98 (ssp-22, ssp-33, lsrk-33), 2.00 (ssp-54, ssp-104, rk4, lsrk-54) | $\min(p, 2) \pm 0.15$ |

Four mutants each fail it: the stage exponent without its abscissa, $\varphi_1 = 1$, the large-$\lambda$ branch
without the $-1$, every abscissa 1.

## Findings

- **Second order, whatever the scheme.** The exponential stage caps the order of the auxiliary field at 2: an
  established, stable construction chosen over higher-order variants that are new (the row-normalised Lawson variant
  of the research, order $\min(p, 3)$, measured but unpublished) or have the wrong steady state (Lawson). Whether the
  PML's second-order error reaches the interior is a question for a reflection test in the field solver.
