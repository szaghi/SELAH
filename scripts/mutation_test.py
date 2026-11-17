#!/usr/bin/env python3
"""Mutation testing of SELAH's test suite: does every test-worthy defect make some test fail? (SETH's tool.)

Each mutant changes one thing in a scratch copy of the tree: a coefficient of the catalogue (Butcher, Shu-Osher,
two-register, 2N, splitting), a kernel body (a sign, a factor, the mask's sign, the empty-sum or a = 0 paths), the box,
the storage counts, the SSP coefficient, the stepper's form logic (the residual slot, the stage weight, the Ketcheson
recipe, the stages on q, the copy of q^n, the splitting sub-flow). The tests are rebuilt (--force-compile) and run
with scripts/run_tests.sh; the mutant is *killed* when at least one test fails. A surviving mutant is a defect the
suite cannot see: the report lists them.

Usage (from the repository root, with fobis and the compilers of the mode on PATH, PENF and FUNDAL fetched):

    python3 scripts/mutation_test.py [--mode tests-gnu] [--workdir DIR] [--only NAME ...]

The working tree is never modified: the copy goes to --workdir (default: a new temporary directory). Exit status 0
when every mutant is killed, 1 otherwise.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

L = 'src/lib/'
# (name, file, anchor, old, new): replace the first `old` after `anchor`. Anchors in code, never in comments: a mutated
# comment is an equivalent mutant.
MUTANTS = [
    # L0: coefficients
    ('ssp22-weight', L + 'selah_schemes.F90', 'bt%b(1) = 0.5_R_P', '0.5_R_P', '0.5001_R_P'),
    ('ssp33-weight', L + 'selah_schemes.F90', 'bt%b(3) = 2._R_P/3._R_P', '2._R_P/3._R_P', '0.6667_R_P'),
    ('rk4-weight', L + 'selah_schemes.F90', 'bt%b = [1._R_P/6._R_P, 1._R_P/3._R_P', '1._R_P/3._R_P, 1._R_P/3._R_P',
     '0.3334_R_P, 1._R_P/3._R_P'),
    ('ssp54-shu-osher-beta', L + 'selah_schemes.F90', 'so%beta(4,3) = ', '0.251891774271694_R_P', '0.2519_R_P'),
    ('ssp104-alpha', L + 'selah_schemes.F90', 'so%alpha(11,1) = 1._R_P/25._R_P', '1._R_P/25._R_P', '1._R_P/26._R_P'),
    ('ssp33-two-register-c', L + 'selah_schemes.F90', 'r = shu_osher_2r_form(stages=3_I4P',
     'c=[1._R_P, 0.25_R_P, 2._R_P/3._R_P]', 'c=[1._R_P, 0.25_R_P, 0.6667_R_P]'),
    ('lsrk33-A', L + 'selah_schemes.F90', 'w = williamson_2n_form(stages=3_I4P', '-5._R_P/9._R_P', '-5.001_R_P/9._R_P'),
    ('lsrk54-B', L + 'selah_schemes.F90', 'b=[1432997174477._R_P', '1432997174477._R_P', '1432997174478._R_P'),
    ('yoshida-w1', L + 'selah_splittings.F90', 'w1 = 1._R_P / (2._R_P - cbrt2)', '(2._R_P - cbrt2)',
     '(2.0001_R_P - cbrt2)'),
    ('blanes-moan-s6-a1', L + 'selah_splittings.F90', 'a(1) =  0.0792036964311957_R_P', '0.0792036964311957_R_P',
     '0.0792_R_P'),
    ('ssp-coefficient-max', L + 'selah_tableaux.F90', 'ssp = min(ssp, so%alpha(i,k) / so%beta(i,k))', 'ssp = min(',
     'ssp = max('),
    # L1: kernels and the box
    ('axpy-sign', L + 'selah_kernel_axpy.INC', 'y(i1,i2,i3,i4,i5) = y(i1,i2,i3,i4,i5) + alpha', '+ alpha', '- alpha'),
    ('lincomb-beta', L + 'selah_kernel_lincomb.INC', '+ beta * y(i1,i2,i3,i4,i5)', '+ beta * y', '+ alpha * y'),
    ('stage-state-dt', L + 'selah_kernel_stage_state.INC', 'd = d + dt * coef(k)', 'd = d + dt * coef(k)',
     'd = d + coef(k)'),
    ('step-update-dq-zero', L + 'selah_kernel_step_update.INC', 'dq(i1,i2,i3,i4,i5) = d', '= d', '= 0.5_R_P * d'),
    ('shu-osher-2r-swap', L + 'selah_kernel_shu_osher_2r_stage.INC', 'a * q0(i1,i2,i3,i4,i5) + b * q(',
     'a * q0(i1,i2,i3,i4,i5) + b * q(', 'b * q0(i1,i2,i3,i4,i5) + a * q('),
    ('williamson-a0-reads', L + 'selah_kernel_williamson_2n_stage.INC', '#ifdef SELAH_A_ZERO',
     'dqr(i1,i2,i3,i4,i5) = dt * r(i1,i2,i3,i4,i5)', 'dqr(i1,i2,i3,i4,i5) = dqr(i1,i2,i3,i4,i5) + dt * r(i1,i2,i3,i4,i5)'),
    ('mask-sign', L + 'selah_kernels_box.F90', 'subroutine axpy_box', 'mask(mc,i2,i3,i4,i5) < 0._R_P',
     'mask(mc,i2,i3,i4,i5) > 0._R_P'),
    ('box-ghost-bound', L + 'selah_box.F90', 'bx%lb = [1_I4P, 1_I4P-ngc', '1_I4P-ngc, 1_I4P-ngc, 1_I4P-ngc',
     '2_I4P-ngc, 1_I4P-ngc, 1_I4P-ngc'),
    ('stage-state-copy-n1', L + 'selah_kernels_box.F90', 'subroutine stage_state_box', 'if (n == 0_I4P) then',
     'if (n <= 1_I4P) then'),
    # L0/L2: storage and the stepper
    ('registers-butcher', L + 'selah_storage.F90', 'n = bt%stages + 1_I4P', 'bt%stages + 1_I4P', 'bt%stages'),
    ('rhs-slot', L + 'selah_stepper.F90', 'slot = self%k_active + 1_I4P', '+ 1_I4P', '+ 0_I4P'),
    ('stage-weight-c', L + 'selah_stepper.F90', 'w = self%bt%b(k)', 'self%bt%b(k)', 'self%bt%c(k)'),
    ('ketcheson-euler-step', L + 'selah_stepper.F90', 'if (k < 10_I4P) call selah_axpy(bx, dt/6._R_P',
     'dt/6._R_P', 'dt/5._R_P'),
    ('ketcheson-combination', L + 'selah_stepper.F90', 'call selah_lincomb(bx, 9._R_P/25._R_P', '9._R_P/25._R_P',
     '9._R_P/24._R_P'),
    ('stages-on-q-all', L + 'selah_stepper.F90', 'self%on_q(k) = all(self%bt%a(k,:) == 0._R_P)',
     'all(self%bt%a(k,:) == 0._R_P)', '.true.'),
    ('two-register-no-qn', L + 'selah_stepper.F90', 'case(FORM_SHU_OSHER_2R, FORM_KETCHESON_2R) !',
     'call selah_copy(bx, q, self%regs(:,:,:,:,:,1), mask, mask_comp)', 'continue'),
    ('splitting-subflow', L + 'selah_stepper.F90', 'self%sc%coefficient(k) * dt, r1, q', 'self%sc%coefficient(k) * dt',
     '1.0001_R_P * self%sc%coefficient(k) * dt'),
]


def run(cmd: list[str], cwd: Path, timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def failing_tests(work: Path, mode: str) -> tuple[bool, list[str]]:
    """Build the tests; return (built, names of the failing test programs)."""
    # --force-compile: FoBiS does not rebuild the dependents of a changed cpp #include (the kernels' .INC bodies)
    b = run(['fobis', 'build', '--mode', mode, '--force-compile'], work)
    if b.returncode != 0 or not (work / 'exe').exists():
        return False, []
    t = run(['bash', 'scripts/run_tests.sh'], work)
    # the runner's lines '  FAIL  <program>' / '  XPASS <program>' (the programs' own '  FAIL  <check>' lines are not)
    failed = [line.split()[1].split('/')[-1] for line in t.stdout.splitlines()
              if line.split()[:1] in (['FAIL'], ['XPASS']) and line.split()[1].split('/')[-1].startswith('selah_')]
    if t.returncode != 0 and not failed:
        failed = ['(run_tests.sh failed)']
    return True, failed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mode', default='tests-gnu')
    ap.add_argument('--workdir', type=Path)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    repo = Path.cwd()
    work = args.workdir or Path(tempfile.mkdtemp(prefix='selah-mutants-'))
    if work.exists():
        shutil.rmtree(work)
    files = run(['git', 'ls-files'], repo).stdout.split() + run(['git', 'ls-files', '--others', '--exclude-standard'],
                                                                repo).stdout.split()
    for f in files:
        if not (repo / f).is_file():
            continue
        dst = work / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / f, dst)
    for dep in ('PENF', 'FUNDAL'):
        shutil.copytree(repo / 'src/third_party' / dep, work / 'src/third_party' / dep, dirs_exist_ok=True)
    built, failed = failing_tests(work, args.mode)
    if not built or failed:
        print(f'baseline is not green ({failed}): fix the suite first')
        return 1
    print(f'baseline green in {work}')
    survivors = []
    selected = [m for m in MUTANTS if not args.only or m[0] in args.only]
    for name, rel, anchor, old, new in selected:
        path = work / rel
        original = path.read_text()
        i = original.find(anchor)
        j = original.find(old, i) if i >= 0 else -1
        if i < 0 or j < 0:
            print(f'{name:28s} NOT APPLIED (anchor or text not found)')
            survivors.append(name)
            continue
        path.write_text(original[:j] + new + original[j + len(old):])
        built, failed = failing_tests(work, args.mode)
        path.write_text(original)
        if not built:
            print(f'{name:28s} does not build (invalid mutant)')
            survivors.append(name)
        elif failed:
            print(f'{name:28s} killed by {", ".join(sorted(set(failed)))}')
        else:
            print(f'{name:28s} SURVIVED')
            survivors.append(name)
        sys.stdout.flush()
    print(f'{len(selected) - len(survivors)} of {len(selected)} mutants killed')
    if survivors:
        print('survivors: ' + ', '.join(survivors))
    return 0 if not survivors else 1


if __name__ == '__main__':
    sys.exit(main())
