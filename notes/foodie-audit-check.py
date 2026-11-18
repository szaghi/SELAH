#!/usr/bin/env python3
"""Check the coefficients of FOODIE's integrators, read from its Fortran sources, in 60-digit arithmetic.

Why this script exists: notes/foodie-audit.md lists what FOODIE implements and SELAH does not. Before a scheme is
ported, its coefficients must be right, checked independently of FOODIE's own tests (its RK tableau audit and its
constant-step convergence and exactness tests). This script parses the coefficients from the sources, with no
transcription by hand, then:

- explicit Runge-Kutta (SSP, embedded, low-storage 2N): every order condition over the rooted trees up to the declared
  order plus one, c = A 1, and for the 2N form the stated C against the c of its Butcher form;
- linear multistep (Adams-Bashforth, Adams-Moulton, BDF, LMM-SSP): the exact order conditions on polynomials;
- LMM-SSP with variable step size, multistep-multistage RK-SSP, linear-SSP RK: the update formulas as FOODIE codes them,
  run on y' = -y^2 (exact solution 1/(1+t)) with exact starting values: observed orders under step halving, with
  constant steps and, for the variable-step schemes, with steps that change;
- literal hygiene: coefficient literals without a real kind (default real) and integer divisions.

Usage: python3 notes/foodie-audit-check.py [--foodie ~/fortran/FOODIE] > notes/foodie-audit-check.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60
TINY = mp.mpf(10) ** -40  # "zero" for exactly typed coefficients


# --------------------------------------------------------------------------------------------------------- parsing
LIT = re.compile(r'(?<![\w.])(\d+\.\d*|\.\d+|\d+)([eEdD][+-]?\d+)?(_[A-Za-z0-9_]+)?(?![\w.])')
ASSIGN = re.compile(r'^\s*(?:self%)?(\w+)\s*(?:\(([^)]*)\))?\s*=\s*(.+)$')


def case_blocks(path: Path) -> dict[str, list[str]]:
    """The statements of each `case('name')` block, continuation lines joined, split on ';', comments removed."""
    lines = path.read_text().splitlines()
    joined: list[str] = []
    buf = ''
    for raw in lines:
        code = raw.split('!')[0].rstrip()
        if code.endswith('&'):
            buf += code[:-1] + ' '
            continue
        joined.append(buf + code)
        buf = ''
    blocks: dict[str, list[str]] = {}
    current = None
    for line in joined:
        m = re.match(r"\s*case\s*\('([^']+)'\)", line)
        if m:
            current = m.group(1)
            blocks[current] = []
            continue
        if re.match(r'\s*(case\s*(default|\()|endselect|end\s+select)', line):
            current = None
            continue
        if current is not None:
            blocks[current].extend(s for s in line.split(';') if s.strip())
    return blocks


def to_python(expr: str, flags: list[str], where: str) -> str:
    """A Fortran constant expression as Python over mpf, flagging default-real literals and integer divisions."""
    e = re.sub(r'real\(\s*([^,()]+?)\s*,\s*kind\s*=\s*\w+\s*\)', r'mpf(\1)', expr, flags=re.I)

    def lit(m: re.Match) -> str:
        num, exp, kind = m.group(1), m.group(2) or '', m.group(3) or ''
        is_real = '.' in num or exp
        if is_real and not kind:
            flags.append(f'{where}: default-real literal {m.group(0)} (single precision)')
        text = num + exp.replace('d', 'e').replace('D', 'e')
        return f"mpf('{text}')" if is_real else f'I({num})'

    e = LIT.sub(lit, e)
    if re.search(r'I\(\d+\)\s*/\s*I\(\d+\)', e):
        flags.append(f'{where}: integer division in {expr.strip()}')
    return e


class FInt(int):
    """A Fortran integer literal: '/' between two of them is an integer division, as in Fortran."""

    def __truediv__(self, other):  # noqa: D105
        if isinstance(other, FInt):
            return FInt(int(self) // int(other))
        return mp.mpf(int(self)) / other

    def __rtruediv__(self, other):  # noqa: D105
        return other / mp.mpf(int(self))


def coefficients(path: Path, scheme: str, flags: list[str]) -> dict[str, dict]:
    """name -> {index tuple (or ()) -> mpf} for the assignments of `scheme`'s case block."""
    out: dict[str, dict] = {}
    for stmt in case_blocks(path)[scheme]:
        m = ASSIGN.match(stmt)
        if not m or m.group(1) in ('stages', 'steps', 'description_', 'registers', 'iterations'):
            continue
        name, idx, expr = m.group(1), m.group(2), m.group(3)
        if 'allocate' in stmt or '%' in expr or not re.search(r'\d', expr):
            continue
        try:
            value = eval(to_python(expr, flags, f'{path.name}:{scheme}:{name}'), {'mpf': mp.mpf, 'I': FInt})  # noqa: S307
        except Exception:  # noqa: BLE001, S112 (an assignment that is not a constant: not a coefficient)
            continue
        key = tuple(int(i) for i in idx.split(',')) if idx else ()
        out.setdefault(name, {})[key] = mp.mpf(value)
    return out


def vector(c: dict, n: int, base: int = 1) -> list:
    return [c.get((i,), mp.mpf(0)) for i in range(base, base + n)]


def matrix(c: dict, n: int) -> list[list]:
    return [[c.get((i, j), mp.mpf(0)) for j in range(1, n + 1)] for i in range(1, n + 1)]


# ------------------------------------------------------------------------------------------- Runge-Kutta, trees
def trees(order: int) -> list[tuple]:
    """Every rooted tree with at most `order` vertices, a tree being the sorted tuple of its subtrees."""
    by_n: dict[int, set] = {1: {()}}
    for n in range(2, order + 1):
        found: set = set()
        # children multiset with total size n-1
        def build(remaining: int, maxsize: int, acc: list, found: set = found):  # noqa: ANN202
            if remaining == 0:
                found.add(tuple(sorted(acc)))
                return
            for size in range(min(remaining, maxsize), 0, -1):
                for t in by_n[size]:
                    build(remaining - size, size, [*acc, t], found)
        build(n - 1, n - 1, [])
        by_n[n] = found
    return [t for n in range(1, order + 1) for t in by_n[n]]


def size(t: tuple) -> int:
    return 1 + sum(size(u) for u in t)


def gamma(t: tuple) -> int:
    g = size(t)
    for u in t:
        g *= gamma(u)
    return g


def order_residuals(a: list[list], b: list, order: int) -> dict[int, mp.mpf]:
    """Max |b . Phi(t) - 1/gamma(t)| over the trees of each order 1..order."""
    s = len(b)
    memo: dict[tuple, list] = {}

    def phi(t: tuple) -> list:
        if t not in memo:
            v = [mp.mpf(1)] * s
            for u in t:
                w = phi(u)
                aw = [mp.fsum(a[i][j] * w[j] for j in range(s)) for i in range(s)]
                v = [v[i] * aw[i] for i in range(s)]
            memo[t] = v
        return memo[t]

    res: dict[int, mp.mpf] = {}
    for t in trees(order):
        r = abs(mp.fsum(b[i] * phi(t)[i] for i in range(s)) - mp.mpf(1) / gamma(t))
        res[size(t)] = max(res.get(size(t), mp.mpf(0)), r)
    return res


def observed_rk_order(a: list[list], b: list, declared: int) -> tuple[int, dict]:
    res = order_residuals(a, b, declared + 1)
    p = 0
    for n in range(1, declared + 2):
        if res[n] < mp.mpf(10) ** -14:
            p = n
        else:
            break
    return p, res


def low_storage_to_butcher(A: list, B: list) -> tuple[list[list], list]:
    """FOODIE's 2N recurrence (dq = A_s dq + dt F(u); u = u + B_s dq) as a Butcher tableau."""
    s = len(A)
    d = [mp.mpf(0)] * s
    w = [mp.mpf(0)] * s
    rows = []
    for k in range(s):
        rows.append(list(w))  # stage k+1 is evaluated at u_{k}
        e = [mp.mpf(1) if j == k else mp.mpf(0) for j in range(s)]
        d = [A[k] * d[j] + e[j] for j in range(s)]
        w = [w[j] + B[k] * d[j] for j in range(s)]
    return rows, w


def fmt(x) -> str:
    return mp.nstr(x, 3, min_fixed=0, max_fixed=0) if x != 0 else '0'


# --------------------------------------------------------------------------------------------- linear multistep
def lmm_order(residual, pmax: int = 20) -> int:
    """Largest p with residual(q) = 0 for q = 0..p (residual(q): the scheme applied to y = tau^q)."""
    p = -1
    for q in range(pmax + 1):
        if abs(residual(q)) < TINY:
            p = q
        else:
            break
    return p


def pw(x, q):
    return mp.mpf(1) if q == 0 else mp.mpf(x) ** q


# ---------------------------------------------------------------------------------------- the y' = -y^2 runs
def f(y):
    return -y * y


def exact(t):
    return 1 / (1 + t)


def observed_orders(step, k: int, t_end, ns: list[int], ratios=None) -> list:
    """Observed orders of a k-step method `step(history_u, history_dt, dt) -> u_new` on y' = -y^2."""
    errs = []
    for n in ns:
        h = mp.mpf(t_end) / n
        dts = [h * (ratios[i % len(ratios)] if ratios else 1) for i in range(n + k)]
        times = [mp.mpf(0)]
        for i in range(n + k):
            times.append(times[-1] + dts[i])
        hist_u = [exact(times[i]) for i in range(k)]
        hist_dt = [dts[i - 1] if i > 0 else dts[0] for i in range(k)]  # the step that produced each history value
        for i in range(k - 1, n + k - 1):
            u = step(hist_u, hist_dt, dts[i])
            hist_u = [*hist_u[1:], u]
            hist_dt = [*hist_dt[1:], dts[i]]
        errs.append(abs(hist_u[-1] - exact(times[n + k - 1])))
    return [mp.log(errs[i] / errs[i + 1]) / mp.log(mp.mpf(ns[i + 1]) / ns[i]) for i in range(len(ns) - 1)]


# ----------------------------------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--foodie', default=str(Path('~/fortran/FOODIE').expanduser()))
    args = ap.parse_args()
    lib = Path(args.foodie) / 'src/lib'
    flags: list[str] = []
    print(f'FOODIE sources: {lib}  (mpmath {mp.__version__}, {mp.mp.dps} digits)\n')

    # ---- explicit RK in Butcher form
    print('## Explicit Runge-Kutta: order conditions (max residual per tree order)\n')
    for fname, prefix in (('foodie_integrator_runge_kutta_ssp.F90', 'runge_kutta_ssp_stages_'),
                          ('foodie_integrator_runge_kutta_embedded.F90', 'runge_kutta_emd_stages_')):
        path = lib / fname
        for scheme in case_blocks(path):
            if not scheme.startswith(prefix):
                continue
            s, declared = (int(x) for x in re.findall(r'\d+', scheme[len(prefix) - 7:])[:2])
            c = coefficients(path, scheme, flags)
            a = matrix(c.get('alph', {}), s)
            gam = vector(c.get('gamm', {}), s)
            cerr = max(abs(mp.fsum(a[i]) - gam[i]) for i in range(s))
            cols = [1, 2] if 'emd' in scheme else [None]
            for col in cols:
                if col is None:
                    b = vector(c.get('beta', {}), s)
                else:
                    b = [c.get('beta', {}).get((i, col), mp.mpf(0)) for i in range(1, s + 1)]
                p, res = observed_rk_order(a, b, declared)
                tag = scheme + ('' if col is None else f' b(:,{col})')
                print(f'{tag:48s} declared {declared:2d} observed {p:2d}   c-A1 {fmt(cerr)}   '
                      + ' '.join(f'{n}:{fmt(res[n])}' for n in sorted(res)))
    # ---- low storage 2N
    path = lib / 'foodie_integrator_runge_kutta_low_storage.F90'
    for scheme in case_blocks(path):
        s, declared = (int(x) for x in re.findall(r'\d+', scheme)[:2])
        c = coefficients(path, scheme, flags)
        A, B, C = (vector(c.get(n, {}), s) for n in ('A', 'B', 'C'))
        a, b = low_storage_to_butcher(A, B)
        cerr = max(abs(mp.fsum(a[i]) - C[i]) for i in range(s))
        p, res = observed_rk_order(a, b, declared)
        print(f'{scheme:48s} declared {declared:2d} observed {p:2d}   C-A1 {fmt(cerr)}   '
              + ' '.join(f'{n}:{fmt(res[n])}' for n in sorted(res)))

    # ---- linear multistep, tabulated
    print('\n## Linear multistep: exact order on polynomials '
          '(declared = steps for AB, steps+1 for AM, steps for BDF)\n')
    path = lib / 'foodie_integrator_adams_bashforth.F90'
    for scheme in case_blocks(path):
        k = int(scheme.rsplit('_', 1)[1])
        b = vector(coefficients(path, scheme, flags).get('b', {}), k)
        tau = [s - k for s in range(1, k + 1)]

        def res_ab(q, b=b, tau=tau, k=k):
            return 1 - q * mp.fsum(b[i] * pw(tau[i], q - 1) for i in range(k)) if q else mp.mpf(0)
        p = lmm_order(res_ab)
        print(f'{scheme:32s} declared {k:2d} observed {p:2d}')
    path = lib / 'foodie_integrator_adams_moulton.F90'
    for scheme in case_blocks(path):
        k = int(scheme.rsplit('_', 1)[1])
        b = vector(coefficients(path, scheme, flags).get('b', {}), k + 1, base=0)
        tau = [(s + 1) - k for s in range(k)] if k > 0 else []

        def res_am(q, b=b, tau=tau, k=k):
            if q == 0:
                return mp.mpf(0)
            return 1 - q * (mp.fsum(b[s] * pw(tau[s], q - 1) for s in range(k)) + b[k])
        print(f'{scheme:32s} declared {k + 1:2d} observed {lmm_order(res_am):2d}')
    path = lib / 'foodie_integrator_backward_differentiation_formula.F90'
    for scheme in case_blocks(path):
        k = int(scheme.rsplit('_', 1)[1])
        c = coefficients(path, scheme, flags)
        a = vector(c.get('a', {}), k)
        bb = c.get('b', {}).get((), mp.mpf(0))
        tau = [s - k for s in range(1, k + 1)]

        def res_bdf(q, a=a, bb=bb, tau=tau, k=k):
            return pw(1, q) + mp.fsum(a[i] * pw(tau[i], q) for i in range(k)) - (q * bb if q else 0)
        print(f'{scheme:32s} declared {k:2d} observed {lmm_order(res_bdf):2d}')
    path = lib / 'foodie_integrator_lmm_ssp.F90'
    for scheme in case_blocks(path):
        k, declared = (int(x) for x in re.findall(r'\d+', scheme)[:2])
        c = coefficients(path, scheme, flags)
        a, b = vector(c.get('a', {}), k), vector(c.get('b', {}), k)
        tau = [s - k for s in range(1, k + 1)]

        def res_ssp(q, a=a, b=b, tau=tau, k=k):
            return pw(1, q) - mp.fsum(a[i] * pw(tau[i], q) for i in range(k)) \
                - (q * mp.fsum(b[i] * pw(tau[i], q - 1) for i in range(k)) if q else 0)
        print(f'{scheme:32s} declared {declared:2d} observed {lmm_order(res_ssp):2d}   '
              f'sum a = {fmt(mp.fsum(a))}, a >= 0: {all(x >= 0 for x in a)}, b >= 0: {all(x >= 0 for x in b)}')

    # ---- LMM-SSP with variable step size: FOODIE's formula against the published one
    print("\n## LMM-SSP-VSS on y' = -y^2, t in [0, 1]: observed orders (step halving: 32 -> 64 -> 128 -> 256 steps)\n")
    ns = [32, 64, 128, 256]

    def vss(order: int, foodie: bool):
        def step(hu, hdt, dt):
            k = len(hu)
            if foodie:  # lmm_ssp_vss.F90:248-267: omega and the increment from the stored steps
                ref = hdt[-1]
                om = mp.fsum(hdt[i] / ref for i in range(k - 1))
                h = ref
            else:  # Hadjimichael et al. (2016): omega_j = dt_{n-j}/dt_n, j = 1..k-1, increment dt_n
                om = mp.fsum(hdt[i] / dt for i in range(1, k))
                h = dt
            if order == 2:
                return hu[0] / om**2 + hu[-1] * (om**2 - 1) / om**2 + f(hu[-1]) * h * (om + 1) / om
            return (hu[0] * (3 * om + 2) / om**3 + hu[-1] * (om + 1) ** 2 * (om - 2) / om**3
                    + f(hu[0]) * h * (om + 1) / om**2 + f(hu[-1]) * h * (om + 1) ** 2 / om**2)
        return step

    for k, order in ((2, 2), (3, 2), (3, 3), (4, 3), (5, 3)):
        for ratios, label in ((None, 'constant dt'), ([1, mp.mpf(3) / 2, mp.mpf(2) / 3], 'dt x 1, 1.5, 2/3')):
            o_f = observed_orders(vss(order, True), k, 1, ns, ratios)
            o_p = observed_orders(vss(order, False), k, 1, ns, ratios)
            print(f'lmm_ssp_vss_steps_{k}_order_{order:<2d} {label:17s} FOODIE: ' + ' '.join(fmt(o) for o in o_f)
                  + '   published: ' + ' '.join(fmt(o) for o in o_p))

    # ---- multistep multistage RK-SSP
    print("\n## Multistep-multistage RK-SSP on y' = -y^2, t in [0, 1], constant dt: observed orders\n")
    path = lib / 'foodie_integrator_ms_runge_kutta_ssp.F90'
    for scheme in case_blocks(path):
        k, s, declared = (int(x) for x in re.findall(r'\d+', scheme)[:3])
        c = coefficients(path, scheme, flags)
        D, Ah, A = (c.get(n, {}) for n in ('D', 'Ahat', 'A'))
        Q, Bh, B = (c.get(n, {}) for n in ('Q', 'Bhat', 'B'))

        def step(hu, hdt, dt, k=k, s=s, D=D, Ah=Ah, A=A, Q=Q, Bh=Bh, B=B, rhs=f):  # ms_runge_kutta_ssp.F90:132-155
            z = mp.mpf(0)
            y = [hu[-1]]
            for i in range(2, s + 1):
                v = mp.fsum(D.get((i, m), z) * hu[m - 1] for m in range(1, k + 1))
                v += dt * mp.fsum(Ah.get((i, m), z) * rhs(hu[m - 1]) for m in range(1, k))
                v += dt * mp.fsum(A.get((i, j), z) * rhs(y[j - 1]) for j in range(1, i))
                y.append(v)
            u = mp.fsum(Q.get((m,), z) * hu[m - 1] for m in range(1, k + 1))
            u += dt * mp.fsum(Bh.get((m,), z) * rhs(hu[m - 1]) for m in range(1, k))
            u += dt * mp.fsum(B.get((j,), z) * rhs(y[j - 1]) for j in range(1, s + 1))
            return u
        o = observed_orders(step, k, 1, [8, 16, 32, 64] if declared > 4 else ns)
        print(f'{scheme:44s} declared {declared:2d} observed ' + ' '.join(fmt(x) for x in o)
              + f'   sum Q = {fmt(mp.fsum(Q.values()))}')
        if declared > 4:  # order 8 needs errors above the floor of its typed digits: y' = -4 y^2, errors printed
            def rhs4(y):
                return -4 * y * y

            def exact4(t):
                return 1 / (1 + 4 * t)
            errs, prev = [], None
            for n in (16, 32, 64, 128, 256):
                h = mp.mpf(1) / n
                hu = [exact4(i * h) for i in range(k)]
                for _ in range(n):
                    hu = [*hu[1:], step(hu, None, h, rhs=rhs4)]
                e = abs(hu[-1] - exact4((n + k - 1) * h))
                errs.append(f'{n}: {fmt(e)}' + ('' if prev is None else f' ({fmt(mp.log(prev / e) / mp.log(2))})'))
                prev = e
            print(f"{'':44s} y' = -4y^2, steps: error (order): " + ', '.join(errs))

    # ---- linear SSP RK: the stability polynomial against exp(z), and the integer factorial
    print('\n## Linear-SSP RK: R(z) - exp(z) = O(z^(p+1)) for y\' = z y (p from the ratio at z = 1e-3, 5e-4)\n')
    for variant in ('order_s', 'order_s_1'):
        for s in (2, 3, 4, 5, 8, 12, 13, 14):
            alpha = [mp.mpf(0)] * (s + 1)

            def fact_i32(i):  # Fortran default-integer product, wrapped to 32 bits as two's complement
                v = 1
                for j in range(1, i + 1):
                    v = (v * j) & 0xFFFFFFFF
                return v - (1 << 32) if v >= 1 << 31 else v
            if variant == 'order_s':  # lssp.F90:297-307
                alpha[1] = mp.mpf(1)
                for i in range(2, s + 1):
                    alpha[i] = mp.mpf(1) / fact_i32(i)
                    for j in range(i - 1, 1, -1):
                        alpha[j] = mp.mpf(1) / (j - 1) * alpha[j - 1]
                    alpha[1] = 1 - mp.fsum(alpha[2:i + 1])
            else:  # lssp.F90:275-288
                alpha[1], alpha[2] = mp.mpf(0), mp.mpf(1)
                for i in range(3, s + 1):
                    alpha[i] = mp.mpf(2) / i * alpha[i - 1]
                    for j in range(i - 1, 1, -1):
                        alpha[j] = mp.mpf(2) / (j - 1) * alpha[j - 1]
                    alpha[1] = 1 - mp.fsum(alpha[2:i + 1])

            def R(z, s=s, alpha=alpha, variant=variant):  # lssp.F90:321-331 / 343-352
                h = mp.mpf(1) / 2 if variant == 'order_s_1' else mp.mpf(1)
                st = [mp.mpf(1)]
                for _ in range(2, s + 1):
                    st.append(st[-1] * (1 + h * z))
                st[-1] = st[-1] * (1 + h * z)
                return mp.fsum(alpha[i + 1] * st[i] for i in range(s))
            e1, e2 = (abs(R(z) - mp.exp(z)) for z in (mp.mpf('1e-3'), mp.mpf('5e-4')))
            p = mp.log(e1 / e2) / mp.log(2) - 1
            print(f'runge_kutta_lssp_{variant:9s} stages {s:2d}: order {fmt(p):>5s}   '
                  f'R(-1)-exp(-1) {fmt(R(-1) - mp.exp(-1))}   {s}! as int32: {fact_i32(s)}')

    print('\n## Literal hygiene\n')
    for path in sorted(lib.glob('foodie_integrator_*.F90')):
        for scheme in case_blocks(path):
            coefficients(path, scheme, flags)
    for line in sorted(set(flags)) or ['none']:
        print(line)
    return 0


if __name__ == '__main__':
    sys.exit(main())
