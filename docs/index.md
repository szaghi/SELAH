---
layout: home

hero:
  name: SELAH
  text: Stage-wise Evolution Library for Accelerators and HPC
  tagline: Time integrators as stage-wise operators for block-structured field solvers, on CPUs and GPUs, in modern Fortran. Every claim verified.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: SELAH, one drop falling on still water and four rings spreading from it
  actions:
    - theme: brand
      text: Quick start
      link: /guide/quickstart
    - theme: alt
      text: Tutorial
      link: /tutorial/
    - theme: alt
      text: Verification
      link: /verification/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/SELAH

features:
  - icon:
      light: /icons/tableaux.svg
      dark: /icons/tableaux-dark.svg
      width: 44
      height: 44
    title: A verified catalogue
    details: Eight Runge–Kutta schemes and four splittings. Order conditions in quadruple precision, stability regions, SSP coefficients, all checked against closed forms and published values.
    link: /guide/schemes
    linkText: Schemes
  - icon:
      light: /icons/kernels.svg
      dark: /icons/kernels-dark.svg
      width: 44
      height: 44
    title: Kernels over your arrays
    details: Stage operations on block-structured fields with ghost cells, in two layouts, masked or not. They write the box and nothing else.
    link: /guide/kernels
    linkText: Kernels
  - icon:
      light: /icons/stepper.svg
      dark: /icons/stepper-dark.svg
      width: 44
      height: 44
    title: A stepper you drive
    details: begin_stage, end_stage, finish_step, around your own ghost fill and residual, under a stage contract an AMR orchestrator can rely on.
    link: /guide/stepper
    linkText: Stepper
  - icon:
      light: /icons/low-storage.svg
      dark: /icons/low-storage-dark.svg
      width: 44
      height: 44
    title: Low storage, same answer
    details: Two-register, Williamson 2N and Ketcheson forms in 2 registers instead of s + 1, equal to the Butcher form to round-off.
    link: /guide/concepts#storage-forms
    linkText: Storage forms
  - icon:
      light: /icons/cpu-gpu.svg
      dark: /icons/cpu-gpu-dark.svg
      width: 44
      height: 44
    title: CPU and GPU, one source
    details: One kernel source for the host and OpenACC through FUNDAL. On an A30 the stepper costs nothing over its kernels, and the step residual comes in the update's pass.
    link: /guide/performance
    linkText: Performance
  - icon:
      light: /icons/order-conditions.svg
      dark: /icons/order-conditions-dark.svg
      width: 44
      height: 44
    title: Evidence, not adjectives
    details: 21 verification programs, a golden file, 28 of 28 mutants killed, an A30, three compilers, two precisions. Every figure regenerated from the programs.
    link: /verification/
    linkText: Verification
---

<script setup>
import { withBase } from 'vitepress'
import { VPTeamMembers } from 'vitepress/theme'

const authors = [
  {
    avatar: 'https://github.com/szaghi.png',
    name: 'Stefano Zaghi',
    title: 'Chief Yak Shaver & Keeper of the Time Step',
    desc: 'Believes every problem is one more library away from being solved, and keeps proving it.',
    links: [{ icon: 'github', link: 'https://github.com/szaghi' }],
  },
  {
    avatar: withBase('/avatars/claude.svg'),
    name: 'Claude',
    org: 'Anthropic',
    orgLink: 'https://claude.ai',
    title: 'Tireless Rubber Duck & Professional Pedant',
    desc: "AI pair programmer; checks the order conditions so humans don't have to, and has an opinion on every comma.",
  },
]
</script>

<div class="selah-home">

## See it run

<div class="selah-casts">
<figure class="selah-cast" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-light.png" aria-label="SSP(5,4) at the SSP limit and above its threshold"><source src="/verification/casts/ssp-tvd-light.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/ssp-tvd-dark.png" aria-label="SSP(5,4) at the SSP limit and above its threshold"><source src="/verification/casts/ssp-tvd-dark.webm" type="video/webm"><source src="/verification/casts/ssp-tvd-dark.mp4" type="video/mp4"></video>
<figcaption>Strong stability: SSP(5,4) at ν = C keeps the total variation from growing; above its threshold it does not. <a href="./verification/properties/ssp-tvd">P1</a></figcaption>
</figure>
<figure class="selah-cast" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-light.png" aria-label="Phase portrait, Strang against rk4"><source src="/verification/casts/splitting-phase-light.webm" type="video/webm"><source src="/verification/casts/splitting-phase-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/splitting-phase-dark.png" aria-label="Phase portrait, Strang against rk4"><source src="/verification/casts/splitting-phase-dark.webm" type="video/webm"><source src="/verification/casts/splitting-phase-dark.mp4" type="video/mp4"></video>
<figcaption>Symplecticity: the Strang splitting stays on a closed orbit, rk4 spirals in. <a href="./verification/properties/splitting-energy">P3</a></figcaption>
</figure>
</div>

## Time integrators

<p class="selah-lead">Every scheme the public API exports: pass the constant to <code>stepper%initialize</code> with one of its storage forms (registers in brackets). <a href="./guide/schemes">Schemes</a> has the tableaux and their derivations.</p>

| Runge–Kutta | Name | Stages | Order | SSP coefficient | Storage forms | Source |
|---|---|---|---|---|---|---|
| `SCHEME_EULER`   | `euler`   | 1  | 1 | 1           | Butcher (2), two-register (2) | forward Euler |
| `SCHEME_SSP_22`  | `ssp-22`  | 2  | 2 | 1           | Butcher (3), two-register (2) | Shu & Osher (1988) |
| `SCHEME_SSP_33`  | `ssp-33`  | 3  | 3 | 1           | Butcher (4), two-register (2) | Shu & Osher (1988) |
| `SCHEME_SSP_54`  | `ssp-54`  | 5  | 4 | 1.508180049 | Butcher (6) | Spiteri & Ruuth (2002) |
| `SCHEME_SSP_104` | `ssp-104` | 10 | 4 | 6           | Butcher (11), Ketcheson (2) | Ketcheson (2008) |
| `SCHEME_RK4`     | `rk4`     | 4  | 4 | —           | Butcher (5) | Kutta (1901) |
| `SCHEME_LSRK_33` | `lsrk-33` | 3  | 3 | —           | Butcher (4), Williamson 2N (2) | Williamson (1980) |
| `SCHEME_LSRK_54` | `lsrk-54` | 5  | 4 | —           | Butcher (6), Williamson 2N (2) | Carpenter & Kennedy (1994) |

| Splitting | Name | Flows | Order | Storage form | Source |
|---|---|---|---|---|---|
| `SCHEME_STRANG`          | `strang`          | 3  | 2 | splitting (1) | Strang (1968) |
| `SCHEME_YOSHIDA_4`       | `yoshida-4`       | 7  | 4 | splitting (1) | Yoshida (1990) |
| `SCHEME_BLANES_MOAN_S6`  | `blanes-moan-s6`  | 13 | 4 | splitting (1) | Blanes & Moan (2002), S6 |
| `SCHEME_BLANES_MOAN_S10` | `blanes-moan-s10` | 21 | 6 | splitting (1) | Blanes & Moan (2002), S10 |

## Quick start

<p class="selah-lead">Build the library and run the verification programs, or pull SELAH into your own project.</p>

::: code-group

```bash [Build and test]
git clone https://github.com/szaghi/SELAH && cd SELAH
fobis fetch                                                   # PENF and FUNDAL, into src/third_party/
fobis build --mode static-gnu                                 # lib/gnu/static/libselah.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh     # the verification programs
```

```ini [FoBiS dependency]
# fobos of your project; then: fobis fetch && fobis build
[dependencies]
deps_dir = src/third_party
PENF     = https://github.com/szaghi/PENF :: branch=master    # SELAH's own dependencies: fobis fetch is not transitive
FUNDAL   = https://github.com/szaghi/FUNDAL :: branch=main
SELAH    = https://github.com/szaghi/SELAH :: branch=main
```

```toml [fpm dependency]
# fpm.toml of your project
[dependencies]
SELAH = { git = "https://github.com/szaghi/SELAH", branch = "main" }
```

:::

<p class="selah-lead">One step, driven stage by stage (<a href="./guide/quickstart">the whole program</a>):</p>

```fortran
use selah
type(stepper) :: st

call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q))
call st%begin_step(bx, q)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, q)
   u => st%stage_state() ; r => st%rhs_register()
   call fill_ghosts(u)            ! yours
   call residual(u, r)            ! yours
   call st%end_stage(k, bx, dt, q)
enddo
call st%finish_step(bx, dt, q)
```

## Authors

</div>

<VPTeamMembers size="medium" :members="authors" />

<div class="selah-home">

## License

<p class="selah-lead">SELAH is available under your choice of license.</p>

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](https://github.com/szaghi/SELAH/blob/main/LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](https://github.com/szaghi/SELAH/blob/main/LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](https://github.com/szaghi/SELAH/blob/main/LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](https://github.com/szaghi/SELAH/blob/main/LICENSE.mit.md) |

</div>
