# Companion states, particle pushers, charge-conserving stage currents and PML auxiliary fields: the theory behind SELAH #4

*2026-11-20, night selah45. The theory check of szaghi/SELAH#4 (Step 1). Labels: **[wiki: bundle/path]** an llm-wiki
concept page (`/home/stefano/llm-wiki`); **[paper: DOI]** the primary source (DOIs checked on Crossref; full texts not
re-read unless stated); **[computed]** derived or computed for this note, at 40 digits by the research agent's scripts
(kept with the night's scratch files) and, where a test exists, reproduced by a SELAH test; **[unverified]** general
knowledge. The research draft was written by an agent; SELAH's tests are the verification.*

The consumer is ADAM's PRISM (szaghi/adam#134, L3, I1, I2, I6): its PIC particles and PML auxiliary fields stepped in
lockstep with the field stepper, its charge-conserving current, its imaginary-axis constant.

## 0. What the sources cover

| Topic | Coverage |
|---|---|
| Boris (non-relativistic), volume preservation; RK4's secular drift against Boris's bounded error | [wiki: fusion-mhd-simulation/methods/particle-pushers.md] (Liu et al. 2022, Hoelzl et al. 2021); [wiki: plasma-propulsion-simulation/methods/particle-pushers.md] |
| Vay (2008), Higuera–Cary (2017) | **not in the wiki**: primary literature |
| Esirkepov, Villasenor–Buneman | bibliography entries only (div-free-amr) |
| Charge conservation with multi-stage RK PIC | **absent from every bundle; no paper found** (wiki, Crossref, web: not exhaustive) |
| PML | **not covered** (bibliography entries only, no ingestion) |
| Stability function, imaginary-axis inclusion | [wiki: ode-integration/definitions/stability-function.md] |
| Lawson, exponential RK, ETD | [wiki: ode-integration/methods/lawson-methods.md], [wiki: ode-integration/methods/exponential-integrators.md] (Hochbruck & Ostermann [paper: 10.1017/S0962492910000048]) |

## 1. Companion states

PRISM's particles (`(8, np)` on the host, `(np, 8)` on the device, components 1:6 evolving, `np` changing every step)
and PML auxiliary fields (six face arrays, rebuilt at a regrid) are states stepped with the field's scheme, step and
stage loop, their residuals built inside the field residual (two-way coupling within each stage). A Runge–Kutta step
of the coupled system $(y, z)' = (f(y,z), g(y,z))$ is the same scheme applied to each part with the same stage loop
[unverified; the method of lines]: so a companion is a **second stepper with the field stepper's scheme and form**,
driven through the same calls. The stepper's box already gives what the issue asks: an active extent changing between
steps (the box's upper bound, inside registers allocated for a capacity), a component sub-range (the box's component
range), either layout (the device's `(np, 8)` is a scalar rank-2 field whose second index is the component). Stepped
together the parts are bitwise the system stepped as one state when the kernels sum in the same order (they do: the
same kernels, element by element) [computed; `selah_test_companion`].

## 2. PIC particles in the field's Runge–Kutta step

**Stability.** Linearised in uniform $B$ the velocity has eigenvalues $0, \pm i\omega_c$ (and the coupled modes
$\pm i\omega_{UH}$, $\pm i c k$) [unverified]: every one must satisfy $|\lambda|\Delta t \le \beta_I$, the
imaginary-axis bound (the largest $\beta$ with $|R(iy)| \le 1$ on $[0, \beta]$) [wiki:
ode-integration/definitions/stability-function.md]. From SELAH's tableaux, as the smallest positive root of odd
multiplicity of $P(y) = |R(iy)|^2 - 1$ (0 when the lowest-order term of $P$ is positive) [computed; identical to
`selah_vv_stability`'s bisection]:

| scheme | leading term of $\lvert R(iy)\rvert^2 - 1$ | $\beta_I$ | $\beta_R$ |
|---|---|---|---|
| euler, ssp-22 | $+y^2$, $+y^4/4$ | 0 (unstable for any $\omega\Delta t > 0$) | 2 |
| ssp-33, lsrk-33 | $-y^4/12$ | $\sqrt 3$ = 1.7320508 | 2.512745 |
| rk4 | $-y^6/72$ | $2\sqrt 2$ = 2.8284271 | 2.785294 |
| ssp-54 | $-0.0049335\,y^6$ | 3.2783556 | 5.331473 |
| lsrk-54 | $-7y^6/1800$ | 3.3407180 | 4.656757 |
| ssp-104 | $-y^6/3240$ | 4.9214531 | 13.917047 |

PRISM's hard-coded 3.278 for SSP(5,4) is right to four digits; $\sqrt 3$ for SSP(3,3) is exact.

**Gyration.** Per step the amplitude is $|R(i\omega_c\Delta t)|$ and the phase $\arg R$: every stable scheme damps the
gyration energy (rk4: $\Delta E/E \approx -(\pi/36)(\omega_c\Delta t)^5$ per gyration) and none preserves phase-space
volume; Boris keeps $|v|$ exactly [computed; wiki: fusion-mhd-simulation/methods/particle-pushers.md].

**E×B drift.** In uniform $E \perp B$ the drift velocity $v_E = E\times B/B^2$ is an equilibrium of the velocity ODE,
and any explicit Runge–Kutta step maps an equilibrium onto itself (every stage residual is 0): the drift is exact for
every scheme and form [computed].

## 3. Two-level pushers: Boris, Vay, Higuera–Cary

Staggered leapfrog: $x$ at $t^n$, $v$ (or $u = \gamma v$) at $t^{n-1/2}$; $E$, $B$ at $x^n$; $\varepsilon = q\Delta t/(2m)$,
$\gamma(u) = \sqrt{1 + |u|^2/c^2}$, SI units. Each formula is the explicit solution of its implicit defining equation
(residual $\le 10^{-40}$ at 40 digits) [computed]; the relativistic ones are reconstructed from the papers' defining
equations, not checked against their printed formulas.

- **Boris** (Boris 1970, NRL proceedings; Birdsall & Langdon): $v^- = v^{n-1/2} + \varepsilon E$, $t = \varepsilon B$,
  $s = 2t/(1+|t|^2)$, $v' = v^- + v^-\times t$, $v^+ = v^- + v'\times s$, $v^{n+1/2} = v^+ + \varepsilon E$,
  $x^{n+1} = x^n + \Delta t\,v^{n+1/2}$. Solves $(v^{n+1/2} - v^{n-1/2})/\Delta t = (q/m)[E + \bar v\times B]$.
  Volume-preserving, not symplectic [paper: 10.1063/1.4818428]; $|v|$ exact in pure $B$; rotation angle
  $2\arctan(\omega_c\Delta t/2)$ per step; $v_E$ an exact fixed point; second order.
- **Relativistic Boris**: the same on $u$, with $t = (\varepsilon/\gamma(u^-))B$ and
  $x^{n+1} = x^n + \Delta t\,u^{n+1/2}/\gamma(u^{n+1/2})$. Volume-preserving; fails the force-free test
  $E + v\times B = 0$ and drifts at the wrong $E\times B$ velocity ($v^*/v_E - 1$ = 2.28e-2 at $E/cB = 0.5$,
  $qB\Delta t/m = 1$).
- **Vay** [paper: 10.1063/1.2837054]: the average velocity in the magnetic term;
  $u' = u^{n-1/2} + 2\varepsilon E + \varepsilon(u^{n-1/2}/\gamma)\times B$, $\tau = \varepsilon B$,
  $u^* = u'\cdot\tau/c$, $\sigma = \gamma(u')^2 - |\tau|^2$,
  $\gamma^{n+1/2} = \sqrt{(\sigma + \sqrt{\sigma^2 + 4(|\tau|^2 + u^{*2})})/2}$, $t = \tau/\gamma^{n+1/2}$,
  $s = 1/(1+|t|^2)$, $u^{n+1/2} = s[u' + (u'\cdot t)t + u'\times t]$. Exact in the force-free case; not
  volume-preserving ($\det - 1 = -1.2\cdot10^{-4}$ at $q\Delta t/m = 1$).
- **Higuera–Cary** [paper: 10.1063/1.4979989]: $\gamma$ at the average momentum; $u^- = u^{n-1/2} + \varepsilon E$,
  $\tau$, $u^*$, $\sigma$ from $u^-$, $\gamma_{new}$ as Vay's, $t = \tau/\gamma_{new}$, $s = 1/(1+|t|^2)$,
  $\bar u = s[u^- + (u^-\cdot t)t + u^-\times t]$, $u^{n+1/2} = \bar u + \bar u\times t + \varepsilon E$.
  Volume-preserving and exact in the force-free case and for the $E\times B$ drift.

These are not Runge–Kutta stage schemes: a separate kind of integrator in SELAH, a kernel over a particle array, with
no tableau. ADAM's deleted PIC leapfrog (three-level, $E$ not divided by $\varepsilon_0$) is not ported.

## 4. Charge-conserving stage currents

With a compatible pair ($\mathrm{div}_h\,\mathrm{curl}_h = 0$) and the Butcher form,
$D_{s+1} - D^n = \Delta t \sum_{r\le s} a_{s+1,r}(\mathrm{curl}_h H_r - J_r)$, so Gauss's law holds at stage $s+1$ iff
$\sum_{r \le s} a_{s+1,r} J_r$ is a current that conserves the charge from $q^n$ to $Y_{s+1}$. Esirkepov's current
$H(q^n \to X)$ does that for one straight displacement [paper: 10.1016/S0010-4655(00)00228-9]. With the look-ahead
$\tilde X_s = x^n + \Delta t\sum_{r\le s} a_{s+1,r} V_r$ (the position of $Y_{s+1}$; $b$ in place of $a_{s+1,\cdot}$ at
the last stage, giving $x^{n+1}$), computable at stage $s$ before the field residual (the position residual is the
stage velocity, a state value), the stage current

$$J_s = \Big(\tilde H_s - \sum_{r<s} a_{s+1,r} J_r\Big) / a_{s+1,s}, \qquad \tilde H_s = H(q^n \to \tilde X_s),$$

gives $\sum_{r\le s} a_{s+1,r} J_r = \tilde H_s$, hence Gauss's law at every stage, and $\sum_k b_k J_k = H(q^n \to
q^{n+1})$ at the step [computed]. It needs $a_{s+1,s} \ne 0$ ($s < S$) and $b_S \ne 0$: true for the whole catalogue;
the inversion amplifies the individual $J_s$ (by up to $\|L^{-1}\|_\infty$ = 25 for ssp-104), so only their weighted
sums are physical currents. **No paper constructs this**: PRISM's construction is original, §4 is its proof. In 2-D
and 3-D Esirkepov's current follows the straight chord, not the curved orbit: the coupling is second order in time
whatever the scheme's order (magnetisation deficit $1 - \sin\theta/\theta$ per gyration) [computed: analytic, not
measured in a PIC run]. Villasenor & Buneman [paper: 10.1016/0010-4655(92)90169-Y], Umeda et al.
[paper: 10.1016/S0010-4655(03)00437-5] give other single-displacement currents; Jacobs & Hesthaven
[paper: 10.1016/j.jcp.2005.09.008] pair RK fields and particles with divergence cleaning instead.

## 5. PML auxiliary fields

$\psi' = g(F) - \lambda\psi$, $g = (\hat\sigma/\kappa^2)\partial_n F$, $\lambda = \alpha + \hat\sigma/\kappa$ (CLASSIC:
$\kappa = 1$, $\alpha = 0$; CFS [paper: 10.1109/75.544545]; Bérenger [paper: 10.1006/jcph.1994.1159]; Bermúdez
[paper: 10.1016/j.jcp.2006.09.018]). **The source is stiff too** ($g \propto \hat\sigma$): the stiffness is not
confined to the decay [computed].

**Explicit stages (PRISM today).** Order $p$, exact equilibria; but the coupled stable set is **not** the rectangle
$\{k\Delta t \le \beta_I\}\times\{\lambda\Delta t \le \beta_R\}$: for rk4 at $\hat\sigma\Delta t = 2.7$ (with $\kappa = 1$,
$\alpha = 0$) the largest stable $k\Delta t$ is 0.89, not 2.83, and 0 at 3 [computed; single Fourier mode of the 1-D
CFS-PML system]. PRISM's two separate limits are therefore insufficient where its typical profiles sit
($\hat\sigma_{max}\Delta t \approx 1.6$–3.2) [computed].

**Exponential treatment.** In Hochbruck & Ostermann's format with $z = -\lambda\Delta t$:
$\psi_i = e^{c_i z}\psi^n + \Delta t\sum_{j<i} a_{ij}(z) g_j$, $\psi^{n+1} = e^z\psi^n + \Delta t\sum_j b_j(z) g_j$;
equilibria (and a frozen source) are exact iff $\sum_j b_j(z) = \varphi_1(z)$ and
$\sum_j a_{ij}(z) = c_i\varphi_1(c_i z)$ [wiki: ode-integration/methods/exponential-integrators.md]. Candidates
[computed on the 4×4 single-mode system]:

| | coefficients | order with the field tableau | frozen source / equilibria | large $\lambda\Delta t$ |
|---|---|---|---|---|
| Lawson (integrating factor) [wiki: lawson-methods.md] | $a_{ij}e^{(c_i-c_j)z}$, $b_je^{(1-c_j)z}$ | $p$ (non-stiff) | **wrong** (moves a static field by up to 0.2) | overflows for non-monotone $c$ (ssp-33, ssp-54, ssp-104) |
| F: exponential Euler per stage on the RK-averaged source (the RK analogue of Roden–Gedney's CPML) | $a_{ij}\varphi_1(c_iz)$, $b_j\varphi_1(z)$ | 2 | exact | stable for every catalogue scheme |
| N: row-normalised Lawson (**new**, proposed by the research) | $c_i\varphi_1(c_iz)\,a_{ij}e^{-c_jz}/\sum_k a_{ik}e^{-c_kz}$, same for $b$ | $\min(p, 3)$ | exact | stable; ψ transfer function 10–100× better than F |

Every exponential variant lowers the stable CFL inside a stiff PML (rk4: to 0.72 $\beta_I$ with F, 0.43 with N), so the
caller takes $\Delta t$ from the coupled limit, not from two separate ones [computed].

**SELAH's choice.** F, implemented first: it is the established construction (Roden & Gedney's recursion
[paper: 10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A] is its one-stage form), exact for a frozen source and
static equilibria, stable for every scheme, a pointwise map at stage formation; its order 2 bounds the PML's own time
error, not the interior's. N is new (measured on one family, not proved): it is recorded here, not implemented
tonight; the choice between F, N and PRISM's explicit stages is the owner's (szaghi/SELAH#4 Step 5.2). The
plane-wave reflection test needs a stencil: it belongs to ADAM; SELAH tests the stencil-free surrogates (a frozen
source, a decay, the order with a smooth source).
