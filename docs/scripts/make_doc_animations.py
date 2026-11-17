#!/usr/bin/env python3
"""Render the animations of SELAH's verification documentation from the frames of the verification programs.

Modelled on SETH's (itself on ADAM's docs/scripts/make_gallery.py): each animation is drawn twice, on the two grounds
of the docs (light washi and dark silver ink, docs/scripts/selah_palette.py), with no title in the frame and a small
time label under the panels. Per animation and theme it writes, into --out:

    <name>-<theme>.webm   VP9, constant quality (crf 36), 720 px wide
    <name>-<theme>.mp4    H.264, yuv420p, faststart
    <name>-<theme>.png    the poster: the last frame

Animations and the data they read (written by scripts/run_vv.sh into --data):

    ssp-tvd           ssp-tvd-frames-{at-c,above}.dat         SSP(5,4), first-order upwind, a square wave: at nu = C
                                                              the total variation never grows; above the threshold
                                                              oscillations appear and grow
    splitting-phase   splitting-phase-{strang,rk4}.dat         the oscillator's phase portrait at dt = 0.6: the
                                                              splitting stays on a closed curve, rk4 spirals inwards

Usage: make_doc_animations.py --data DIR [--out DIR] [--only NAME ...] [--fps N] [--width PX]
Needs numpy, matplotlib and ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import selah_palette as P  # noqa: E402


@dataclass(frozen=True)
class Theme:
    """The colours of one ground of the docs, and the line colours (selah_palette)."""

    name: str
    bg: str
    bg_alt: str
    text: str
    mesh: str
    lines: tuple[str, ...]  # exact/reference first, then the curves


THEMES = (
    Theme('light', P.WASHI, P.WASHI_ALT, P.INK_LIGHT, P.MESH_LIGHT, (P.INK_LIGHT, *P.LINES_LIGHT)),
    Theme('dark', P.SILVER, P.SILVER_ALT, P.INK_DARK, P.MESH_DARK, (P.INK_DARK, *P.LINES_DARK)),
)


def rc(th: Theme) -> dict:
    return {'figure.facecolor': th.bg, 'savefig.facecolor': th.bg, 'axes.facecolor': th.bg_alt,
            'axes.edgecolor': th.mesh, 'axes.labelcolor': th.text, 'text.color': th.text, 'xtick.color': th.mesh,
            'ytick.color': th.mesh, 'xtick.labelcolor': th.text, 'ytick.labelcolor': th.text, 'grid.color': th.mesh,
            'legend.facecolor': th.bg, 'legend.edgecolor': th.mesh}


def encode(tmp: Path, nframes: int, base: Path, fps: int, width: int) -> None:
    scale = f'scale={width}:-2:flags=lanczos'
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                    '-vf', scale, '-c:v', 'libvpx-vp9', '-crf', '36', '-b:v', '0', '-pix_fmt', 'yuv420p',
                    f'{base}.webm'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                    '-vf', scale, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                    f'{base}.mp4'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', str(tmp / f'f{nframes - 1:04d}.png'), '-vf', scale,
                    f'{base}.png'], check=True)
    print(f'  {base}.{{webm,mp4,png}}')


def read_frames(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """x, t, u per frame, exact per frame."""
    rows = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith('#')]
    x = np.array(rows[0].split(), float)
    n = x.size
    data = np.array([r.split() for r in rows[1:]], float)
    return x, data[:, 0], data[:, 1:1 + n], data[:, 1 + n:1 + 2 * n]


def ssp_tvd(data: Path, out: Path, fps: int, width: int) -> None:
    files = [data / 'ssp-tvd-frames-at-c.dat', data / 'ssp-tvd-frames-above.dat']
    if not all(p.exists() for p in files):
        print('ssp-tvd: skipped, missing frames')
        return
    x, t1, u1, e1 = read_frames(files[0])
    _, t2, u2, e2 = read_frames(files[1])
    nus = [ln.split('nu =')[1].strip() for ln in (f.read_text().splitlines()[0] for f in files)]
    n = min(t1.size, t2.size)

    def tv(u: np.ndarray) -> np.ndarray:
        return np.abs(np.diff(np.concatenate([u, u[:, :1]], axis=1), axis=1)).sum(axis=1)

    tvs = (tv(u1) / tv(u1)[0], tv(u2) / tv(u2)[0])
    tv_top = max(1.1, 1.05 * max(tvs[0].max(), tvs[1].max()))
    for th in THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'selah-ssp-tvd-{th.name}-'))
        try:
            with mpl.rc_context(rc(th)):
                for f in range(n):
                    fig, axs = plt.subplots(2, 2, figsize=(9.6, 4.8), dpi=100, gridspec_kw={'height_ratios': [2, 1]})
                    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.11, top=0.93, wspace=0.18, hspace=0.35)
                    for col, (t, u, e, title) in enumerate(((t1, u1, e1, f'ν = C = {nus[0]}'),
                                                            (t2, u2, e2, f'ν = {nus[1]}, above the threshold'))):
                        ax = axs[0, col]
                        ax.plot(x, e[f], '-', color=th.lines[0], lw=0.9, label='exact')
                        ax.plot(x, u[f], '-', color=th.lines[1], lw=1.2, label='SSP(5,4)')
                        ax.set_xlim(0, 1)
                        ax.set_ylim(-0.15, 1.15)
                        ax.grid(alpha=0.25)
                        ax.set_title(title, fontsize=9)
                        ax.legend(loc='upper right', fontsize=8)
                        ax = axs[1, col]
                        ax.axhline(1.0, color=th.lines[0], lw=0.8, ls='--')
                        ax.plot(np.arange(f + 1), tvs[col][:f + 1], '-', color=th.lines[2 + col], lw=1.3)
                        ax.set_xlim(0, n - 1)
                        ax.set_ylim(0.6, tv_top)
                        ax.set_ylabel('TV / TV$_0$', fontsize=8)
                        ax.grid(alpha=0.25)
                    fig.text(0.07, 0.02, f'step {f}', fontsize=9, color=th.text)
                    fig.savefig(tmp / f'f{f:04d}.png')
                    plt.close(fig)
            out.mkdir(parents=True, exist_ok=True)
            encode(tmp, n, out / f'ssp-tvd-{th.name}', fps, width)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def splitting_phase(data: Path, out: Path, fps: int, width: int) -> None:
    files = [data / 'splitting-phase-strang.dat', data / 'splitting-phase-rk4.dat']
    if not all(p.exists() for p in files):
        print('splitting-phase: skipped, missing data')
        return
    a = np.loadtxt(files[0], ndmin=2)
    b = np.loadtxt(files[1], ndmin=2)
    n = min(len(a), len(b))
    stride = max(n // 120, 1)
    th_circle = np.linspace(0, 2 * np.pi, 200)
    frames = list(range(stride, n + 1, stride))
    for th in THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'selah-phase-{th.name}-'))
        try:
            with mpl.rc_context(rc(th)):
                for i, f in enumerate(frames):
                    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), dpi=100)
                    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.16, top=0.9, wspace=0.25)
                    for ax, d, label, slot in ((axs[0], a, 'Strang splitting', 3), (axs[1], b, 'rk4', 4)):
                        ax.plot(np.cos(th_circle), np.sin(th_circle), '-', color=th.lines[0], lw=0.7)
                        ax.plot(d[:f, 2], d[:f, 3], '.', color=th.lines[slot], ms=2)
                        ax.plot(d[f - 1, 2], d[f - 1, 3], 'o', color=th.lines[slot], ms=5)
                        ax.set_xlim(-1.3, 1.3)
                        ax.set_ylim(-1.3, 1.3)
                        ax.set_aspect('equal')
                        ax.set_xlabel('q')
                        ax.set_ylabel('p')
                        ax.set_title(label, fontsize=9)
                        ax.grid(alpha=0.25)
                    fig.text(0.07, 0.03, f't = {a[f - 1, 0]:.1f}  (Δt = 0.6)', fontsize=9, color=th.text)
                    fig.savefig(tmp / f'f{i:04d}.png')
                    plt.close(fig)
            out.mkdir(parents=True, exist_ok=True)
            encode(tmp, len(frames), out / f'splitting-phase-{th.name}', fps, width)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=Path('docs/public/verification/casts'))
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--width', type=int, default=720)
    args = ap.parse_args()
    if shutil.which('ffmpeg') is None:
        print('ffmpeg not found: no animation rendered')
        return 1
    jobs = {'ssp-tvd': ssp_tvd, 'splitting-phase': splitting_phase}
    for name, job in jobs.items():
        if args.only and name not in args.only:
            continue
        print(name)
        job(args.data, args.out, args.fps, args.width)
    return 0


if __name__ == '__main__':
    sys.exit(main())
