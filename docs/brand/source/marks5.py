import math, random
from brush import Ink, stroke, blob, drop, arcpts, resample, catmull, kasure, poly_d
from marks4 import earc

def ring(ink, role, cx, cy, rx, ry, w, seed, gap0=-150, span=305, **kw):
    a0 = math.radians(gap0)
    kw.setdefault("t0", 0.05); kw.setdefault("t1", 0.5)
    stroke(ink, role, earc(cx, cy, rx, ry, a0, a0 + math.radians(span)), w, seed=seed, wob=kw.pop("wob", 0.12), belly=kw.pop("belly", 0.08), pw=kw.pop("pw", 0.45), minw=0.25, **kw)

def trace(ink, x, y0, y1, ws, seed):
    """The drop's fall: a faint dry line, broken twice."""
    L = y1 - y0
    for k, (a, b) in enumerate([(0.0, 0.32), (0.42, 0.68), (0.78, 0.92)]):
        stroke(ink, "ink55", [(x, y0 + L * a), (x + 0.05, y0 + L * b)], 0.45 * ws, t0=0.3, t1=0.3, seed=seed + k, wob=0.3, belly=0.0)

# ------------------------------------------------------------ 5.1 Hamon, refined
def hamon(ws=1.0):
    ink = Ink()
    drop(ink, "acc", 30.0, 9.5, 3.3 * ws)
    trace(ink, 30.0, 15.0, 36.0, ws, 10)
    cx, cy = 30.0, 41.0
    for k, (rx, w, role, g) in enumerate([(4.6, 0.95, "ink", -120), (9.6, 1.0, "ink", -60), (15.6, 0.95, "ink55", -165), (22.6, 0.85, "ink55", -30)]):
        ring(ink, role, cx, cy, rx, rx * 0.34, w * ws, 20 + k, gap0=g)
    blob(ink, "ink", cx, cy, 1.0, 0.5, seed=30)
    return ink

# ------------------------------------------------------------ 5.2 Tsukubai
def tsukubai(ws=1.0):
    ink = Ink()
    cx, cy = 31.0, 41.0
    # the stone basin: a heavy dry-brushed rim
    kasure(ink, "ink", earc(cx, cy, 21.0, 9.5, math.radians(-200), math.radians(150)), 3.4 * ws, seed=40, dry=0.45, strands=7)
    stroke(ink, "ink55", earc(cx, cy + 4.6, 22.0, 9.0, math.radians(10), math.radians(170)), 2.2 * ws, t0=0.2, t1=0.3, seed=41, wob=0.2)
    # water inside: three rings round where the drop falls
    wx, wy = 27.0, 40.0
    for k, (rx, w, role) in enumerate([(3.4, 0.85, "ink"), (7.4, 0.8, "ink"), (12.0, 0.7, "ink55")]):
        ring(ink, role, wx, wy, rx, rx * 0.36, w * ws, 42 + k, gap0=-140 + 70 * k)
    # the bamboo spout (kakehi) from the upper left
    stroke(ink, "ink", [(-1.0, 9.0), (12.0, 12.0), (24.0, 16.0)], 3.0 * ws, t0=0.0, t1=0.02, seed=50, wob=0.06, belly=0.0, minw=0.85, cap=0.9)
    for k, f in enumerate([0.35, 0.75]):
        x, y = -1.0 + 25.0 * f, 9.0 + 7.0 * f
        stroke(ink, "ground", [(x - 0.6, y - 2.0), (x + 0.6, y + 2.0)], 0.6 * ws, t0=0.05, t1=0.05, seed=51 + k, wob=0.2, track=False)
    drop(ink, "acc", 26.2, 25.0, 2.6 * ws)
    trace(ink, 25.6, 18.0, 21.5, ws, 55)
    return ink

# ------------------------------------------------------------ 5.3 Ue kara (from above)
def uekara(ws=1.0):
    ink = Ink()
    cx, cy = 29.5, 29.5
    blob(ink, "acc", cx, cy, 3.4 * ws ** 0.5, 3.4 * ws ** 0.5, seed=60, rough=0.04)
    for k, (r, w, role, g) in enumerate([(8.0, 2.1, "ink", -100), (13.4, 1.8, "ink", 40), (19.0, 1.45, "ink55", -170), (24.5, 1.15, "ink55", 120)]):
        ring(ink, role, cx, cy, r, r, w * ws, 61 + k, gap0=g, span=318)
    return ink

# ------------------------------------------------------------ 5.4 Nidan (two steps)
def nidan(ws=1.0):
    ink = Ink()
    # step n: wide rings, already pale
    for k, (rx, w) in enumerate([(8.0, 0.85), (13.5, 0.8), (19.5, 0.7)]):
        ring(ink, "ink55", 23.0, 46.0, rx, rx * 0.33, w * ws, 70 + k, gap0=-160 + 40 * k)
    # step n+1: the new drop and its first rings, dark
    drop(ink, "acc", 41.0, 11.0, 3.1 * ws)
    trace(ink, 41.0, 16.5, 31.0, ws, 75)
    for k, (rx, w) in enumerate([(4.4, 0.95), (8.6, 1.0)]):
        ring(ink, "ink", 41.0, 36.0, rx, rx * 0.34, w * ws, 78 + k, gap0=-110 + 80 * k)
    blob(ink, "ink", 41.0, 36.0, 1.0, 0.55, seed=81)
    return ink

# ------------------------------------------------------------ 5.5 Mizu no tsuki (moon in the water)
def mizunotsuki(ws=1.0):
    ink = Ink()
    cx, cy, R = 31.0, 40.0, 11.5
    rnd = random.Random(90)
    # the reflection, broken by the ripples into slices
    k = 0
    y = -R + 1.2
    while y < R - 0.6:
        half = math.sqrt(max(0.0, R * R - y * y))
        shift = 2.4 * math.sin(y * 0.9 + 0.6) + rnd.uniform(-0.4, 0.4)
        h = 1.05 if abs(y) < R * 0.7 else 0.8
        stroke(ink, "acc", [(cx - half + shift, cy + y * 0.55), (cx + half + shift, cy + y * 0.55 + 0.1)], h * ws, t0=0.15, t1=0.25, seed=91 + k, wob=0.15, belly=0.1)
        y += 3.1
        k += 1
    # the ripples that break it
    for j, (rx, w, role) in enumerate([(8.0, 0.8, "ink"), (14.0, 0.8, "ink"), (20.0, 0.7, "ink55"), (25.5, 0.65, "ink55")]):
        ring(ink, role, cx, cy, rx, rx * 0.34, w * ws, 110 + j, gap0=-170 + 55 * j)
    # the drop that broke it, small, in ink
    drop(ink, "ink", cx + 1.0, 14.0, 2.2 * ws)
    trace(ink, cx + 1.0, 18.5, 30.0, ws, 120)
    # far bank, one pale line
    stroke(ink, "ink55", [(2.0, 26.0), (20.0, 25.4), (40.0, 26.2), (62.0, 25.6)], 0.6 * ws, t0=0.3, t1=0.3, seed=125, wob=0.3)
    return ink

# ------------------------------------------------------------ 5.6 Ōkan (the crown)
def okan(ws=1.0):
    ink = Ink()
    cx, cy = 30.0, 44.0
    for k, (rx, w, role) in enumerate([(11.0, 1.1, "ink"), (17.5, 0.9, "ink55"), (24.0, 0.75, "ink55")]):
        ring(ink, role, cx, cy, rx, rx * 0.28, w * ws, 130 + k, gap0=-150 + 60 * k)
    # the crown: a rim of spikes with droplets at their tips
    n = 9
    for k in range(n):
        a = math.pi + math.pi * (k + 0.5) / n * 1.0 if False else math.radians(190 + 160 * k / (n - 1))
        bx, by = cx + 7.0 * math.cos(a), cy + 2.2 * math.sin(a)
        h = 6.0 + 2.0 * math.sin(k * 1.7)
        lean = 1.6 * math.cos(a)
        tx, ty = bx + lean, by - h
        stroke(ink, "ink", [(bx, by), (bx + lean * 0.4, by - h * 0.55), (tx, ty)], 1.4 * ws, t0=0.02, t1=0.6, seed=140 + k, wob=0.1, pw=0.6)
        blob(ink, "ink", tx + lean * 0.35, ty - 1.3, 0.75, 0.75, seed=150 + k)
    # front rim of the crown
    stroke(ink, "ink", earc(cx, cy, 7.2, 2.4, math.radians(10), math.radians(170)), 1.5 * ws, t0=0.1, t1=0.1, seed=160, wob=0.1)
    # the Worthington jet and its rebound drop
    stroke(ink, "ink", [(cx, cy - 1.0), (cx + 0.2, cy - 10.0), (cx - 0.1, cy - 18.5)], 2.6 * ws, t0=0.02, t1=0.75, seed=161, wob=0.08, pw=0.6)
    blob(ink, "acc", cx - 0.2, cy - 25.0, 2.6 * ws ** 0.5, 2.9 * ws ** 0.5, seed=162, rough=0.05)
    return ink

# ------------------------------------------------------------ 5.7 Kasure hamon (dry brush)
def kasurehamon(ws=1.0):
    ink = Ink()
    drop(ink, "acc", 30.0, 9.0, 3.2 * ws)
    trace(ink, 30.0, 14.5, 33.0, ws, 170)
    cx, cy = 30.0, 41.0
    for k, (rx, w, dry, g) in enumerate([(5.0, 1.5, 0.75, -120), (10.6, 2.0, 0.55, -50), (17.0, 2.4, 0.4, -170), (24.0, 2.4, 0.25, -20)]):
        a0 = math.radians(g)
        kasure(ink, "ink", earc(cx, cy, rx, rx * 0.36, a0, a0 + math.radians(300)), w * ws, seed=171 + k, dry=dry, strands=6, t0=0.05)
    return ink

# ------------------------------------------------------------ 5.8 Suzuri no umi (the inkstone's sea)
def suzuri(ws=1.0):
    ink = Ink()
    # the inkstone, a heavy dry rim
    x0, x1, y0, y1 = 9.0, 47.0, 21.0, 54.0
    kasure(ink, "ink", [(x0 + 1.0, y0), (x1 - 1.0, y0 - 0.3)], 2.8 * ws, seed=180, dry=0.65, strands=6)
    stroke(ink, "ink", [(x1, y0 + 1.5), (x1 - 0.3, y1 - 1.5)], 2.4 * ws, t0=0.05, t1=0.1, seed=181, wob=0.08, belly=0.0)
    kasure(ink, "ink", [(x1 - 1.0, y1), (x0 + 1.0, y1 + 0.3)], 2.8 * ws, seed=182, dry=0.65, strands=6)
    stroke(ink, "ink", [(x0, y1 - 1.5), (x0 + 0.3, y0 + 1.5)], 2.4 * ws, t0=0.05, t1=0.1, seed=183, wob=0.08, belly=0.0)
    # the sea (umi): a well of ground ink
    sea = earc(28.0, 30.0, 14.0, 5.6, 0, 2 * math.pi, n=60)
    ink.add("ink", poly_d(sea))
    for x, y in sea[::4]:
        ink.pts.append((x, y, 0.8))
    ink.pts.append((28.0, 30.0, 6.0))
    # rings cut into the ink by the drop, in ground colour
    for k, (rx, w) in enumerate([(2.6, 0.9), (6.0, 0.75), (10.0, 0.6)]):
        a0 = math.radians(-150 + 60 * k)
        stroke(ink, "ground", earc(26.0, 30.2, rx, rx * 0.38, a0, a0 + math.radians(300)), w * ws, t0=0.05, t1=0.5, seed=184 + k, wob=0.15, track=False, belly=0.2)
    # the land (oka): pale grinding strokes
    for k in range(3):
        stroke(ink, "ink55", [(15.0 + 2 * k, 41.0 + 3.6 * k), (41.0 - 2 * k, 40.6 + 3.6 * k)], 1.3 * ws, t0=0.2, t1=0.3, seed=190 + k, wob=0.3)
    # the water dropper (suiteki) above, and the drop in colour
    blob(ink, "ink", 41.0, 6.0, 5.0, 3.6, seed=195, rough=0.08)
    stroke(ink, "ink", [(37.0, 6.6), (32.5, 9.0)], 1.6 * ws, t0=0.05, t1=0.3, seed=196, wob=0.1)
    drop(ink, "acc", 26.6, 16.5, 2.4 * ws)
    return ink

MARKS5 = {"hamon": hamon, "tsukubai": tsukubai, "uekara": uekara, "nidan": nidan, "mizunotsuki": mizunotsuki, "okan": okan, "kasurehamon": kasurehamon, "suzuri": suzuri}
