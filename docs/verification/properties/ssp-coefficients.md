---
title: P2 SSP coefficients
icon: ssp
---

# P2 SSP coefficients

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| P2 | Properties | P | CPU (gfortran, nvfortran, ifx; R8P and R4P) | no | `selah_test_schemes`, `selah_test_tableau_forms`, `selah_vv_stability` |

## Case and scope

SELAH never types an SSP coefficient: `butcher(scheme)%ssp_coefficient` is computed from the scheme's Shu–Osher form,
$C = \min_{\beta_{ik} > 0} \alpha_{ik}/\beta_{ik}$, or 0 if any coefficient is negative. The case checks the computed
values against the published optimal ones. It also checks the function on hand-made forms, and the disc property
[A5](../accuracy/stability) derives from $C$.

## Setup

| Scheme | Published $C$ | Source | Effective $C/s$ |
|---|---|---|---|
| `euler`, `ssp-22`, `ssp-33` | 1 | Shu & Osher (1988); Gottlieb & Shu (1998) | 1, 0.5, 0.333 |
| `ssp-54` | 1.508 | Spiteri & Ruuth (2002) | 0.302 |
| `ssp-104` | 6 | Ketcheson (2008) | 0.6 |
| `rk4`, `lsrk-33`, `lsrk-54` | not SSP (0) | | |

Bounds: $|C - 1| \le 64\,\varepsilon$ for euler, ssp-22, ssp-33; $|C - 6| \le 64\,\varepsilon \cdot 6$ for ssp-104;
$|C - 1.508| < 5\cdot10^{-4}$ for ssp-54 (Spiteri & Ruuth publish $C$ to four digits; SELAH computes 1.5082 from their
15-digit coefficients); exactly 0 for rk4 and the 2N schemes.

## Program

```bash
./exe/selah_test_schemes ; ./exe/selah_test_tableau_forms
```

## Results

The values are those of the [stability table](../accuracy/stability) (column $C$), and every one matches.

## Findings

- **$C$ is a property of a representation.** The minimum is taken over the Shu–Osher coefficients SELAH stores, and
  another Shu–Osher form of the same Butcher tableau can give a larger value. For the catalogue the stored forms are the
  optimal ones, so $C$ is the scheme's SSP coefficient. For a tableau you build yourself, 0 means "not computed" as well
  as "not SSP".
- SSP(10,4) has the best effective coefficient in the catalogue: per residual evaluation it advances twice as far as
  SSP(5,4) at its strong-stability limit ([Cookbook](/cookbook/schemes#the-largest-strong-stability-preserving-time-step)).
