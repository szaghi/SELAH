#!/usr/bin/env python3
"""Draw the figures of SELAH's Performance page from the benchmark records kept in notes/.

Why this script exists: the Performance page reports A30 measurements (notes/m3-bench.md, the M3 kernel table, and
notes/m5-bench-data/job-1368116.txt, the raw output of selah_bench_stepper). The figures are drawn from those records,
not retyped, in the light and the dark theme.

Outputs (docs/public/guide/): perf-kernels.png (SELAH against ADAM, kernel by kernel), perf-steps.png (whole steps:
time and field accesses), and their -dark twins.

Usage: make_performance_figures.py [--out DIR] [--theme {light,dark,both}]
Needs numpy and matplotlib.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import selah_palette as P  # noqa: E402

mpl.rcParams['figure.facecolor'] = P.WASHI
mpl.rcParams['savefig.facecolor'] = P.WASHI
import matplotlib.pyplot as plt  # noqa: E402

from figure_theme import THEMES, theme, themes  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DPI = 110


def m3_kernels() -> list[tuple[str, float]]:
    """(operation, SELAH/ADAM) rows of the M3 table that have an ADAM counterpart, the whole step excluded."""
    rows = []
    for line in (ROOT / 'notes/m3-bench.md').read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if len(cells) == 7 and re.fullmatch(r'[0-9.]+', cells[4]) and 'whole' not in cells[0]:
            rows.append((cells[0], float(cells[4])))
    return rows


def m5_steps() -> list[tuple[str, int, float]]:
    """(operation, accesses, ms) of pass 2 of job 1368116."""
    rows, on = [], False
    for line in (ROOT / 'notes/m5-bench-data/job-1368116.txt').read_text().splitlines():
        if line.startswith('pass 2'):
            on = True
            continue
        if on and line.startswith('pass 3'):
            break
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if on and len(cells) == 4 and re.fullmatch(r'[0-9]+', cells[1]):
            rows.append((cells[0], int(cells[1]), float(cells[2])))
    return rows


def kernels_figure(out: Path, dark: bool) -> None:
    rows = m3_kernels()
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    y = np.arange(len(rows))
    ax.barh(y, [r[1] for r in rows], color=P.BLUE_D if dark else P.BLUE, height=0.6)
    ax.axvline(1.0, color='0.4', lw=0.8, ls='--')
    ax.set_yticks(y, [r[0] for r in rows], fontsize=8)
    ax.set_xlim(0.9, 1.05)
    ax.set_xlabel('time, SELAH / ADAM (A30, job 1368089)')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    fig.savefig(out / 'perf-kernels.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def steps_figure(out: Path, dark: bool) -> None:
    rows = [r for r in m5_steps() if 'copy' not in r[0] and 'stage state of stage 1' not in r[0]
            and 'control' not in r[0] and 'stage 1 on q (kernels)' not in r[0]]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw={'width_ratios': [3, 2]})
    y = np.arange(len(rows))
    colours = [(P.UMBER_D if dark else P.UMBER) if r[0].startswith('ADAM') else (P.BLUE_D if dark else P.BLUE)
               for r in rows]
    axs[0].barh(y, [r[2] for r in rows], color=colours, height=0.6)
    for i, r in enumerate(rows):
        axs[0].text(r[2] + 0.05, i, f'{r[2]:.2f} ms', va='center', fontsize=7.5)
    axs[0].set_yticks(y, [r[0] for r in rows], fontsize=7.5)
    axs[0].invert_yaxis()
    axs[0].set_xlabel('time per step [ms] (A30, 256 blocks of 16³ cells, 5 components)')
    axs[0].grid(axis='x', alpha=0.3)
    axs[1].barh(y, [r[1] for r in rows], color=colours, height=0.6)
    axs[1].set_yticks(y, [''] * len(rows))
    axs[1].invert_yaxis()
    axs[1].set_xlabel('field accesses per step')
    axs[1].grid(axis='x', alpha=0.3)
    fig.savefig(out / 'perf-steps.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=Path('docs/public/guide'))
    ap.add_argument('--theme', choices=THEMES, default='both')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for th in themes(args.theme):
        with theme(th):
            for job in (kernels_figure, steps_figure):
                print(f'{job.__name__} ({th})')
                job(args.out, th == 'dark')
    return 0


if __name__ == '__main__':
    sys.exit(main())
