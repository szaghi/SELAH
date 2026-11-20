# Limiters, additive splitting, particle sub-flows and masks: the theory behind SELAH #5

*2026-11-20, night selah45. The theory check of szaghi/SELAH#5 (Step 1): what each feature must guarantee, traced to
the llm-wiki bundles and the primary literature. Labels: **[wiki: bundle/path]** an llm-wiki concept page
(`/home/stefano/llm-wiki`; status `established` unless stated); **[paper: DOI]** the primary source (DOIs checked on
Crossref); **[computed]** derived or computed for this note and reproduced by a SELAH test; **[unverified]** general
knowledge not checked against a source. The research draft was written by an agent from the wiki and the papers; every
number below that a test can check is checked by one (named at each section).*

The consumer is ADAM's FLUME (szaghi/adam#134, I3-I5, I7-I8): its positivity limiter, its radiative-cooling split, its
sink particles and its adapter's masked stage-state copy.

## 1. The forward-Euler step of each stage (`stage_euler_fraction`, `ssp_coefficient`)

**Shu–Osher form.** $Y_1 = u^n$, $Y_i = \sum_{k<i} (\alpha_{ik} Y_k + \Delta t\, \beta_{ik} F(Y_k))$, rows of
$\alpha$ summing to 1. With $\alpha, \beta \ge 0$ and $\beta_{ik} > 0 \Rightarrow \alpha_{ik} > 0$, stage $i$ is a
convex combination of forward-Euler steps $Y_k + \Delta t (\beta_{ik}/\alpha_{ik}) F(Y_k)$
[paper: 10.1016/0021-9991(88)90177-5; wiki: ode-integration/methods/strong-stability-preserving-runge-kutta.md].

**The step of stage k.** The residual of stage $k$ enters forward-Euler steps of $\Delta t\,\beta_{ik}/\alpha_{ik}$,
one per later stage $i$. Let $f_k = \max_{i>k,\ \beta_{ik}>0} \beta_{ik}/\alpha_{ik}$. For a convex admissible set
$\mathcal G$ (positive density and pressure), if $Y_k \in \mathcal G$ and $Y_k + f_k \Delta t F_k \in \mathcal G$, then
$Y_k + \tau F_k = (1 - \tau/\tau_k) Y_k + (\tau/\tau_k)(Y_k + \tau_k F_k) \in \mathcal G$ for every
$0 \le \tau \le \tau_k = f_k\Delta t$; by induction every stage and $u^{n+1}$ are admissible [computed; the convexity
step of wiki: high-order-capturing-schemes/properties/positivity-preserving.md]. The max is sufficient and minimal: a
smaller step leaves the longest forward-Euler step uncertified. Three conditions:

1. $F_k$ is **one** vector, the limited residual, used unchanged by every later stage;
2. what the limiter computes from the step (a flux limiter's low-order backbone, a time-centred source) is computed with
   $\tau_k$: a source time-centred within the forward-Euler step sits at $\tau_k/2$ from the stage state, not at
   $\Delta t/2$;
3. the stages are evaluated as the convex combinations they are. In floating point the Butcher form
   $u^n + \Delta t \sum_j a_{ij} F_j$ sums signed increments and can round below a floor
   (`selah_test_stage_euler`: $-2\cdot10^{-19}$ in ssp-104's Butcher form, exactly 0 in its Ketcheson form)
   [computed]. A limiter is safest in the Shu–Osher, two-register or Ketcheson form.

**Relation to C.** $\max_k f_k = 1/C$, with $C = \min \alpha_{ik}/\beta_{ik}$ the SSP coefficient of the stored
representation; for an optimal representation it is the method's (Kraaijevanger) coefficient
[paper: 10.1007/BF01933264]. The individual $f_k$ depend on the representation: the canonical form of SSP(5,4) puts
every residual into steps of $\Delta t/C = 0.663\,\Delta t$; SELAH's stored form gives stage 1 only $0.392\,\Delta t$
[computed]. The limiter must use the representation SELAH evaluates, so the tests pin that one.

| Scheme | $C$ | $f_k$ | Source |
|---|---|---|---|
| euler, ssp-22, ssp-33 | 1 | 1 every stage | [paper: 10.1016/0021-9991(88)90177-5] |
| ssp-54 | 1.50818004918983 | 0.391752226571890, 0.663050807850946, 0.663050807850949, 0.663050807850948, 0.584438703993960 | Spiteri & Ruuth [paper: 10.1137/S0036142901389025], the digits of Isherwood, Grant & Gottlieb [paper: 10.1137/17M1143290]; [computed] |
| ssp-104 | 6 | 1/6 every stage | Ketcheson [paper: 10.1137/07070485X] |
| rk4, lsrk-33, lsrk-54 | 0 | 0 | not SSP: no forward-Euler decomposition |

In the whole catalogue every ratio $\beta_{ik}/\alpha_{ik}$ of a stage is the same (ssp-54's stage 4 enters stages 5
and 6 with 0.663051 to 15 digits) [computed]: the max is attained by every step of the stage.

**FLUME's choice** (every stage a forward-Euler step of $\Delta t$) is exact for SSP-22/33 and safe but over-limiting
for SSP-54 (by $1/f_k$ = 1.5 to 2.55) and SSP-104 (6): it demands a backbone admissible at $\Delta t$ instead of
$C\,\Delta t_{FE}$, losing SSP-54's larger step, and centres its gravity term at $\Delta t/2$ instead of
$\tau_k/2 = (0.196, 0.332, 0.332, 0.332, 0.292)\,\Delta t$ [computed].

**Limiters.** The Zhang–Shu scaling limiter [paper: 10.1016/j.jcp.2010.08.016 (Euler equations; the scalar
maximum-principle paper is 10.1016/j.jcp.2009.12.030); wiki: high-order-capturing-schemes/methods/positivity-preserving-scaling-limiter.md]
needs $\tau_k$ only in its CFL condition; the flux (blending) limiter of FLUME
[wiki: high-order-capturing-schemes/methods/positivity-preserving-flux-limiter.md] needs it in the update it makes
admissible. Whether limiting inside the stages reduces the order is an **open contradiction** of the wiki
[wiki: high-order-capturing-schemes/properties/positivity-preserving.md ## Contradictions]: SELAH claims positivity,
not order, under limiting.

**SELAH's choice.** `stage_euler_fraction(k)` returns $f_k$, and 0 for a scheme without the decomposition;
`ssp_coefficient()` returns $C$, 0 for those schemes and for a splitting. The draft recommended refusing the query
for a non-SSP scheme; SELAH returns 0 and documents that a caller checks `ssp_coefficient() > 0` before trusting a
limiter (no program stop in a query). Verified by `selah_test_stage_euler` (values, every form, planted forms, a
positivity limiter that is positive with $f_k\Delta t$ and negative with half of it).

## 2. Additive splitting of the field step and a caller's flow (`selah_additive_splitting`)

**Orders.** For $u' = F(u) + G(u)$ with exact sub-flows, Lie is order 1 (local error $\tfrac12 h^2 [L_1, L_2]$ in the
linear case) and Strang order 2 (local error $\tfrac{h^3}{24}([L_1,[L_1,L_2]] + 2[L_2,[L_1,L_2]])$); when the parts
commute both are exact, so a test needs non-commuting parts [wiki: ode-integration/methods/operator-splitting.md;
paper: 10.1137/0705041]. With an RK field step of order $p \ge 2$ in place of the exact flow, Strang keeps order 2
[unverified; measured by `selah_test_additive_splitting`].

**Which flow outside.** The cheap flow (cooling) outside: $\varphi_{h/2} \circ \mathrm{RK}_h \circ \varphi_{h/2}$.
RK is not a flow, so RK outside doubles the field cost; an exact autonomous flow merges across steps
($\varphi_{a}\circ\varphi_{b} = \varphi_{a+b}$), so Strang then costs what Lie does (Strang is conjugate to Lie
[wiki: ode-integration/methods/processing-and-effective-order.md]); and with a stiff operator, putting it outside is
the more accurate ordering (Verwer & Sportisse, CWI report MAS-R9830, 1998; Sportisse [paper: 10.1006/jcph.2000.6495],
read through its abstract and a citing summary). SELAH's default for the palindromes: `first=OPERATOR_FLOW`.

**Order reduction.** In the stiff limit Strang drops to order 1 for the stiff components, in both orderings (Verwer &
Sportisse 1998 §4.3; conditions in Descombes & Massot [paper: 10.1007/s00211-003-0496-3]); boundary conditions reduce
it too [wiki: ode-integration/methods/operator-splitting.md §Order reduction]; splitting shifts steady states by
$O(h)$ (balanced splitting, Speth et al. [paper: 10.1137/120878641]). Townsend's exact cooling
[paper: 10.1088/0067-0049/181/2/391] removes the cooling sub-step's time error, not the splitting error
[wiki: astrophysical-simulation-codes/methods/radiative-cooling-integration.md]. **Measured:** on the cooling model
$u' = a u - k u^2$ (the field step $a u$ by rk4, the exact flow of $-k u^2$, Townsend's form), Lie and Strang keep
orders 1 and 2 with $k u_0$ up to $10^3$ ($k u_0 \Delta t$ up to 125), the error constant bounded as $k$ grows.
SELAH does not generalise this: the guide states that the order may drop with a stiff flow.

**The negative-coefficient barrier.** Every real splitting of order $\ge 3$ has a negative fraction in each operator
[wiki: ode-integration/results/negative-coefficient-barrier.md; papers: Sheng 10.1093/imanum/9.2.199, Suzuki
10.1063/1.529425, Goldman & Kaper 10.1137/0733018, Blanes & Casas 10.1016/j.apnum.2004.10.005]. A backward upwind or
limited field step is unstable; a backward cooling flow heats gas. **SELAH's choice:** `initialize` stops on yoshida-4,
blanes-moan-s6 or -s10 unless the caller declares both operators reversible (`reversible=.true.`);
`selah_test_xfail_additive_backward` checks the stop.

**Time.** Two consistent rules exist. (a) Time as a variable with $t' = 1$ carried by one part, frozen through the
other (the stepper's splitting form uses this: its Euler sub-flows of one partition cannot carry time). (b) Each
operator integrates its own time over its substep, the composition of $(u, t_1, t_2)' = (F(u,t_1),1,0) +
(G(u,t_2),0,1)$, whose exact solution keeps $t_1 = t_2$ [computed]. SELAH's additive splitting takes (b): substep $m$
starts at the sum of the same operator's fractions before $m$; the field step evaluates its stages at
$t_m + c_k h$, the caller's flow integrates from $t_m$. It keeps the order with both terms time-dependent: measured
(`selah_test_additive_splitting`, lie 1, strang 2).

**Production codes** [wiki: astrophysical-simulation-codes]: most apply cooling as a Lie step after the hydro step
(Townsend) or before it (TIGRESS, sub-cycled); Athena++ and IDEFIX put super-time-stepped diffusion in Strang form; no
source in the bundles demonstrates second order for the cooled hydro system
[wiki: methods/radiative-cooling-integration.md, methods/super-time-stepping.md].

## 3. Particle sub-flows (the splitting form on a particle array)

**KDK is Strang** with exact sub-flows, kick $v' = a(x)$ and drift $x' = v$: Störmer–Verlet, symplectic,
time-symmetric, preserving linear and quadratic invariants such as angular momentum under central forces
[wiki: ode-integration/methods/operator-splitting.md "Names"; ode-integration/properties/symplecticity.md]. With a
fixed step the energy error is bounded over exponentially long times; **variable steps destroy this** [wiki: same
pages]. The stepper's splitting form already is this integrator: partition 1 the velocities, partition 2 the
positions, `strang` (the half kick a fraction of SELAH's, not a typed 1/2). No new code; the tests prove it.

**Harmonic oscillator** $x'' = -\omega^2 x$: stable for $\omega h < 2$; the KDK map has trace $2 - (\omega h)^2$ and
determinant 1, so from $x_0 = 1$, $v_0 = 0$ the discrete solution is $x_n = \cos(n\theta)$ with
$\cos\theta = 1 - (\omega h)^2/2$; the phase error per period is $\approx \pi(\omega h)^2/12 = \pi^3/(3N^2)$ with $N$
steps per period [computed; `selah_test_particle_flows`: $x_n$ to $3\cdot10^{-13}$ over 100 periods].

**Kepler** ($GM = 1$, $e = 0.5$): observed orders 2, 4, 4, 6 for strang, yoshida-4, blanes-moan-s6, blanes-moan-s10;
energy error bounded with no drift (measured as the largest error over every step of an orbit: sampled once per
orbit, it follows the precessing perihelion) [computed]. Higher-order compositions need backward flows, harmless for a
particle Hamiltonian (reversible) [paper: 10.1016/0375-9601(90)90092-3; wiki: ode-integration/methods/composition-methods.md].
**Open:** whether SELAH's blanes-moan-s6 digits are Blanes & Moan's general S6 or their Runge–Kutta–Nyström SRKN6
cannot be told on Kepler (an RKN problem); a non-RKN two-part test is needed.

**Gas and sinks.** Federrath et al. sub-cycle the sinks with the gas frozen; TIGRESS advances positions before the MHD
step (Quinn et al.'s KDK [paper: 10.1088/0004-6256/139/2/803]); Cholla's step order is hydro, sources, kick and drift,
deposit, Poisson solve, final kick [wiki: astrophysical-simulation-codes/methods/sink-particles.md,
methods/particle-mesh-gravity.md]. Coupled through the additive splitting (the gas the field step, a KDK step of the
sinks the flow), Lie is order 1 and Strang order 2 (`selah_test_particle_flows`). Accretion is irreversible: only
non-negative fractions, Lie or Strang.

## 4. Masked cells in the stage states (`masked_q`)

**No source addresses what an RK stage register holds in masked cells** (checked: the immersed-boundary-methods,
ode-integration and high-order-capturing-schemes bundles). The nearest: the ghost-fluid method re-extrapolates its
ghost values after every stage [wiki: immersed-boundary-methods/methods/ghost-fluid-method.md]; halo values frozen
across stages limit the accuracy there to second order (Goswami et al.)
[wiki: ode-integration/methods/runge-kutta-stage-boundary-treatment.md]; cells uncovered by a moving body have no
valid history [wiki: immersed-boundary-methods/definitions/fresh-cell-problem.md], so the mask is constant within a
step.

**SELAH's policy** [computed]: the masked cells of every stage state hold $q^n$, which is the Runge–Kutta solution of
$u' = 0$ there in every form. The forms other than Butcher keep it automatically (the state is $q$, untouched in
masked cells); the Butcher form's register did not (its masked cells held the previous stage's values), and
`masked_q=.true.` makes the stage-state kernel write $q^n$ there in the same pass. Values the active stencils read
from masked cells that must evolve within the step (immersed-boundary ghost values rebuilt from the fluid) remain the
caller's to refresh at every stage: $q^n$ is a default, not a boundary condition. Verified by
`selah_test_stepper_options` (bitwise $q^n$ in every form, bitwise ADAM's own copy) and `selah_test_kernels`,
`selah_test_device_kernels` (both layouts).
