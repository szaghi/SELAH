#!/usr/bin/env python3
"""Draw the figures of SELAH's conventions: the box, the two layouts, the stage contract, the register maps.

Why this script exists: the Concepts page explains what a kernel touches, how a field is laid out in memory, who writes
what during a stage and which register holds what in each storage form. These are pictures of SELAH's own rules (the
box of selah_box.F90, the slot logic of selah_stepper.F90), drawn here once and regenerated on demand, in the light
and the dark theme (figure_theme.theme("dark") writes the <name>-dark.png twins).

Outputs (docs/public/guide/): conv-box.png, conv-layouts.png, conv-contract.png, conv-registers.png and their twins.

Usage: make_convention_figures.py [--out DIR] [--theme {light,dark,both}]
Needs numpy and matplotlib.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use('Agg')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import selah_palette as P  # noqa: E402

mpl.rcParams['figure.facecolor'] = P.WASHI
mpl.rcParams['savefig.facecolor'] = P.WASHI
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

from figure_theme import THEMES, theme, themes  # noqa: E402

DPI = 110
C = {}  # the colours of the current theme (main)


def palette(dark: bool) -> None:
    C.clear()
    C.update(interior=P.BLUE_D if dark else P.BLUE, read='#5d7f9e' if dark else '#a9c3dc',
             ghost=P.UMBER_D if dark else P.UMBER,
             pool='0.45' if dark else '0.72', caller=P.GREEN_D if dark else P.GREEN, stepper=P.BLUE_D if dark else P.BLUE,
             rhs=P.PLUM_D if dark else P.PLUM, ink=P.INK_DARK if dark else P.INK_LIGHT)


def box_figure(out: Path) -> None:
    """A 2-D slice (i, j) of a pool of 4 blocks, 4 x 3 interior cells, 2 ghost layers, 3 active blocks."""
    ni, nj, ngc, nb, active = 4, 3, 2, 4, 3
    fig, ax = plt.subplots(figsize=(10.5, 3.0))
    w = ni + 2 * ngc + 1.2
    for b in range(nb):
        x0 = b * w
        for i in range(-ngc + 1, ni + ngc + 1):
            for j in range(-ngc + 1, nj + ngc + 1):
                interior = 1 <= i <= ni and 1 <= j <= nj
                if b >= active:
                    face, hatch = C['pool'], None
                elif interior:
                    face, hatch = C['interior'], None
                else:
                    face, hatch = 'none', '////'
                ax.add_patch(Rectangle((x0 + i + ngc - 1, j + ngc - 1), 1, 1, facecolor=face, edgecolor=C['ink'],
                                       lw=0.5, hatch=hatch, alpha=0.85 if face != 'none' else 1.0))
        ax.text(x0 + (ni + 2 * ngc) / 2, nj + 2 * ngc + 0.35, f'block {b + 1}' + ('' if b < active else ' (pool)'),
                ha='center', fontsize=9, color=C['ink'])
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor=C['interior'], label='the box: lo..hi, written by the kernels'))
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor='none', hatch='////', edgecolor=C['ink'],
                           label='ghost cells 1-ngc..0 and n+1..n+ngc: never written by SELAH'))
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor=C['pool'], label='blocks past blocks_number: never touched'))
    ax.set_xlim(-0.5, nb * w)
    ax.set_ylim(-0.3, nj + 2 * ngc + 1.0)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=8, frameon=False)
    fig.savefig(out / 'conv-box.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def layouts_figure(out: Path) -> None:
    """The first 16 values in memory of a field with 2 components and 4 blocks, in both layouts."""
    fig, axs = plt.subplots(2, 1, figsize=(10, 2.4))
    for ax, title, seq in ((axs[0], 'LAYOUT_V_FIRST (v,i,j,k,b): components innermost, the CPU layout',
                            [(f'v{v}', f'b1', v) for i in range(8) for v in (1, 2)][:16]),
                           (axs[1], 'LAYOUT_V_LAST (b,i,j,k,v): blocks innermost, the device layout (coalesced over b)',
                            [(f'b{b}', 'v1', b) for i in range(4) for b in (1, 2, 3, 4)][:16])):
        cols = [C['interior'], C['caller'], C['rhs'], C['ghost']]
        for n, (lab, _, k) in enumerate(seq):
            ax.add_patch(Rectangle((n, 0), 1, 1, facecolor=cols[(k - 1) % 4], edgecolor=C['ink'], lw=0.6, alpha=0.8))
            ax.text(n + 0.5, 0.5, lab, ha='center', va='center', fontsize=8, color='white')
        ax.text(16.3, 0.5, '→ address', va='center', fontsize=8, color=C['ink'])
        ax.set_title(title, fontsize=9, loc='left')
        ax.set_xlim(0, 18)
        ax.set_ylim(0, 1)
        ax.axis('off')
    fig.subplots_adjust(hspace=0.8)
    fig.savefig(out / 'conv-layouts.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def contract_figure(out: Path) -> None:
    """One step, two stages shown: who writes what and when."""
    fig, ax = plt.subplots(figsize=(12, 3.0))
    lanes = {'stepper': 2.0, 'caller': 1.0}
    x = 0.0

    def block(x0, w, lane, text, colour):
        ax.add_patch(FancyBboxPatch((x0, lanes[lane] - 0.36), w, 0.72, boxstyle='round,pad=0.02,rounding_size=0.08',
                                    facecolor=colour, edgecolor=C['ink'], lw=0.6, alpha=0.9))
        ax.text(x0 + w / 2, lanes[lane], text, ha='center', va='center', fontsize=7.5, color='white')

    block(x, 1.7, 'stepper', 'begin_step', C['stepper'])
    x += 1.9
    for k in (1, 2):
        block(x, 1.9, 'stepper', f'begin_stage({k})', C['stepper'])
        block(x + 2.1, 1.7, 'caller', 'ghost fill', C['caller'])
        block(x + 4.0, 1.9, 'caller', 'residual into\nrhs_register()', C['rhs'])
        block(x + 6.1, 1.7, 'stepper', f'end_stage({k})', C['stepper'])
        ax.add_patch(Rectangle((x + 2.0, 0.45), 4.0, 2.1, facecolor='none', edgecolor=C['ink'], lw=0.7, ls='--'))
        ax.text(x + 4.0, 2.62, f"stage {k}: the caller's window", ha='center', fontsize=7.5, color=C['ink'])
        x += 8.0
    ax.text(x + 0.3, 2.0, '...', ha='center', va='center', fontsize=12, color=C['ink'])
    x += 0.9
    block(x, 1.8, 'stepper', 'finish_step\n(q, and dq)', C['stepper'])
    for lane, y in lanes.items():
        ax.text(-0.3, y, lane, ha='right', va='center', fontsize=9, color=C['ink'])
    ax.set_xlim(-1.6, x + 2.0)
    ax.set_ylim(0.3, 3.0)
    ax.axis('off')
    fig.savefig(out / 'conv-contract.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def registers_figure(out: Path) -> None:
    """What each register holds during stage k, per form: the slot logic of selah_stepper.F90."""
    forms = [
        ('Butcher, s = 4: s + 1 registers', 4, lambda k: (['U (stage state)'] + [f'F{j}' for j in range(1, 5)], k)),
        ('two-register (ssp-33)', 3, lambda k: (['q^n', 'R'], 1)),
        ('Williamson 2N (lsrk-54)', 5, lambda k: (['ΔQ', 'R'], 1)),
        ('Ketcheson (ssp-104)', 10, lambda k: (['q2', 'R'], 1)),
    ]
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.3), gridspec_kw={'width_ratios': [5, 2, 2, 2]})
    for ax, (title, s, slots) in zip(axs, forms):
        for k in range(1, s + 1):
            names, rhs = slots(k)
            for j, nm in enumerate(names):
                if title.startswith('Butcher'):
                    if j == 0:
                        colour = C['stepper']
                    elif j == k:
                        colour = C['rhs']
                    elif j < k:
                        colour = C['read']
                    else:
                        colour = 'none'
                else:
                    colour = C['rhs'] if j == rhs else C['stepper']
                ax.add_patch(Rectangle((j, s - k), 1, 1, facecolor=colour, edgecolor=C['ink'], lw=0.5,
                                       alpha=0.8 if colour != 'none' else 1.0))
                if k == 1:
                    ax.text(j + 0.5, s + 0.2, nm, ha='center', va='bottom', fontsize=7, color=C['ink'], rotation=30)
        ax.set_xlim(0, len(slots(1)[0]))
        ax.set_ylim(0, s + 1.6)
        ax.set_yticks([s - k + 0.5 for k in range(1, s + 1)], [f'k={k}' for k in range(1, s + 1)], fontsize=6.5)
        ax.set_xticks([])
        ax.set_title(title, fontsize=8.5)
        ax.set_facecolor('none')
        for sp in ax.spines.values():
            sp.set_visible(False)
    axs[0].add_patch(Rectangle((0, 0), 0, 0, facecolor=C['stepper'], label='written by the stepper'))
    axs[0].add_patch(Rectangle((0, 0), 0, 0, facecolor=C['rhs'], label="rhs_register(): the caller's residual"))
    axs[0].add_patch(Rectangle((0, 0), 0, 0, facecolor=C['read'], label='earlier residuals, read'))
    fig.legend(loc='lower center', ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.savefig(out / 'conv-registers.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=Path('docs/public/guide'))
    ap.add_argument('--theme', choices=THEMES, default='both')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for th in themes(args.theme):
        palette(th == 'dark')
        with theme(th):
            for job in (box_figure, layouts_figure, contract_figure, registers_figure):
                print(f'{job.__name__} ({th})')
                job(args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
