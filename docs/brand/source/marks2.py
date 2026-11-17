import math, random
from brush import Ink, stroke, blob, drop, arcpts, resample, catmull, kasure, smoke, offset_path

def spiral_pts(cx, cy, r0, r1, turns, a0=0.0, step=0.4, sq=0.0):
    pts = []
    th_max = 2 * math.pi * turns
    n = int(th_max * r1 / step) + 20
    for i in range(n + 1):
        th = th_max * i / n
        r = r0 + (r1 - r0) * th / th_max
        a = a0 + th
        c, s = math.cos(a), math.sin(a)
        if sq:
            m = max(abs(c), abs(s)) ** sq
            c, s = c / m, s / m
        pts.append((cx + r * c, cy + r * s))
    return pts

def split_dense(path, step=0.3):
    P, tot = resample(path, step)
    return P

def ash(ink, P, w, seed, gap=3, seg=9):
    k = 0
    rnd = random.Random(seed)
    while k < len(P) - 4:
        j = min(len(P) - 1, k + seg + rnd.randint(-2, 3))
        stroke(ink, "ink55", P[k:j + 1], w * rnd.uniform(0.75, 1.05), t0=0.25, t1=0.3, seed=seed + k, wob=0.3, belly=0.0, smooth=False, minw=0.3)
        k = j + gap + rnd.randint(0, 2)

def roof(ink, role, x, y, s, seed, ws=1.0):
    stroke(ink, role, [(x - 1.6 * s, y + 0.65 * s), (x - 0.2 * s, y - 0.5 * s), (x + 1.6 * s, y + 0.6 * s)], 0.9 * s * ws, t0=0.15, t1=0.35, seed=seed, wob=0.12)
    stroke(ink, role, [(x - 1.0 * s, y + 1.35 * s), (x + 1.0 * s, y + 1.3 * s)], 0.62 * s * ws, t0=0.2, t1=0.25, seed=seed + 1, wob=0.1)

# =============================================================== B · incense
def senkodai(ws=1.0):
    """B2.1 One stick in a censer: ash curl, ember, thread ties ahead, smoke."""
    ink = Ink()
    # censer: rim, bowl, feet
    stroke(ink, "ink", [(13.0, 45.5), (31.0, 44.8), (49.0, 45.6)], 2.2 * ws, t0=0.06, t1=0.12, seed=1, wob=0.1)
    kasure(ink, "ink", [(14.5, 47.0), (19.0, 54.0), (31.0, 57.2), (43.0, 54.0), (47.5, 47.0)], 5.2 * ws, seed=2, dry=0.62, strands=8, t0=0.05)
    for k, (x0, x1) in enumerate([(19.5, 17.5), (42.5, 44.5)]):
        stroke(ink, "ink", [(x0, 55.5), (x1, 60.5)], 1.7 * ws, t0=0.1, t1=0.5, seed=3 + k, wob=0.1)
    # ash bed in the bowl
    stroke(ink, "ink55", [(17.0, 47.8), (31.0, 48.6), (45.0, 47.6)], 1.2 * ws, t0=0.3, t1=0.3, seed=6, wob=0.3)
    # the stick, unburnt below the ember
    stroke(ink, "ink", [(31.0, 47.6), (31.4, 34.0), (31.8, 21.6)], 1.45 * ws, t0=0.02, t1=0.02, seed=7, wob=0.05, belly=0.0, minw=0.9, cap=0.95)
    # thread ties: the stage pegs still ahead of the fire
    for k, y in enumerate([27.0, 31.5, 36.0, 42.0]):
        x = 31.0 + (47.6 - y) / (47.6 - 21.6) * -0.8 + 0.8
        stroke(ink, "ink", [(x - 2.0, y + 0.2), (x, y - 0.3), (x + 2.2, y + 0.4)], 0.8 * ws, t0=0.3, t1=0.3, seed=10 + k, wob=0.15)
        stroke(ink, "ink55", [(x + 2.8, y + 0.5), (x + 3.6, y + 2.6)], 0.7 * ws, t0=0.2, t1=0.6, seed=20 + k, wob=0.2)
    # the ember, and the ash curl drooping from it
    stroke(ink, "acc", [(31.8, 22.4), (31.9, 18.6)], 2.6 * ws, t0=0.5, t1=0.35, seed=30, wob=0.08, pw=0.4)
    stroke(ink, "ink55", [(31.9, 18.8), (33.2, 16.0), (35.6, 15.2), (37.0, 17.0)], 1.5 * ws, t0=0.05, t1=0.4, seed=31, wob=0.2)
    # smoke
    smoke(ink, "ink55", [(31.6, 16.5), (29.0, 12.0), (32.5, 8.0), (29.5, 4.0), (32.5, 0.5)], 1.6 * ws, seed=32)
    smoke(ink, "ink55", [(32.4, 14.0), (36.5, 10.5), (35.0, 6.0), (39.5, 2.0)], 1.0 * ws, seed=33, strands=1)
    return ink

def uzumaki(ws=1.0):
    """B2.2 The coiled incense of summer nights: burning from outside in."""
    ink = Ink()
    cx, cy = 29.5, 30.0
    S = spiral_pts(cx, cy, 2.6, 21.5, 3.3, a0=math.pi * 0.75)
    P = split_dense(S)
    n = len(P)
    burnt = int(n * 0.86)
    stroke(ink, "ink", P[:burnt - 19], 2.7 * ws, t0=0.01, t1=0.002, seed=40, wob=0.06, belly=0.0, smooth=False, minw=0.5, step=0.3)
    stroke(ink, "acc", P[burnt - 20:burnt + 3], 3.1 * ws, t0=0.15, t1=0.6, seed=41, wob=0.06, smooth=False, pw=0.4)
    ash(ink, P[burnt + 5:], 2.3 * ws, 42, gap=4, seg=10)
    # the stand: a short brushed pin through the centre
    stroke(ink, "ink", [(cx - 3.0, cy + 3.4), (cx + 3.2, cy - 3.2)], 1.6 * ws, t0=0.2, t1=0.3, seed=45, wob=0.1)
    # fallen ash below the coil
    for k, (x, y, r) in enumerate([(14.0, 55.5, 1.2), (19.5, 57.5, 0.9), (25.0, 56.2, 1.4), (31.5, 58.0, 0.8)]):
        blob(ink, "ink55", x, y, r * 1.5, r, seed=50 + k, rough=0.2)
    ex, ey = P[burnt - 4]
    smoke(ink, "ink55", [(ex, ey - 1.0), (ex - 3.0, ey - 5.0), (ex + 1.0, ey - 9.0), (ex - 1.5, ey - 13.0)], 1.1 * ws, seed=55)
    return ink

def koban(ws=1.0):
    """B2.3 The incense box: a square spiral of powder, hour plates ahead."""
    ink = Ink()
    # the lacquer box, double brushed frame
    kasure(ink, "ink", [(5.5, 6.0), (50.5, 5.6)], 3.4 * ws, seed=60, dry=0.7, strands=6)
    kasure(ink, "ink", [(51.0, 5.0), (50.6, 50.5)], 3.4 * ws, seed=61, dry=0.7, strands=6)
    kasure(ink, "ink", [(50.0, 51.0), (5.6, 50.6)], 3.4 * ws, seed=62, dry=0.7, strands=6)
    kasure(ink, "ink", [(5.0, 51.0), (5.4, 5.5)], 3.4 * ws, seed=63, dry=0.7, strands=6)
    S = spiral_pts(28.2, 28.2, 3.2, 17.5, 2.6, a0=math.pi * 0.25, sq=0.85)
    P = split_dense(S[::-1])
    n = len(P)
    b = int(n * 0.40)
    ash(ink, P[:b - 2], 1.8 * ws, 64, gap=3, seg=8)
    stroke(ink, "acc", P[b - 2:b + 12], 2.6 * ws, t0=0.55, t1=0.12, seed=65, wob=0.08, smooth=False, pw=0.4)
    stroke(ink, "ink", P[b + 12:], 2.0 * ws, t0=0.01, t1=0.25, seed=66, wob=0.1, smooth=False, minw=0.5, belly=0.0)
    # hour plates (jifuda) standing in the trail ahead
    for k, f in enumerate([0.52, 0.64, 0.77, 0.93]):
        x, y = P[min(n - 1, int(n * f))]
        stroke(ink, "ink", [(x, y - 4.2), (x + 0.1, y + 0.6)], 1.3 * ws, t0=0.05, t1=0.1, seed=70 + k, wob=0.1, belly=0.0, minw=0.8, cap=0.8)
        stroke(ink, "ink", [(x - 1.3, y - 4.4), (x + 1.4, y - 4.6)], 1.0 * ws, t0=0.2, t1=0.3, seed=75 + k, wob=0.1)
    return ink

def kemuri(ws=1.0):
    """B2.4 Smoke as the clock: four curls rising from one ember."""
    ink = Ink()
    # small holder and a short stick, low left
    blob(ink, "ink", 15.0, 55.0, 6.5, 2.6, seed=80, rough=0.08)
    stroke(ink, "ink", [(15.0, 53.0), (17.6, 43.0)], 1.4 * ws, t0=0.02, t1=0.02, seed=81, wob=0.05, minw=0.9, cap=0.9)
    stroke(ink, "acc", [(17.6, 43.2), (18.3, 40.4)], 2.4 * ws, t0=0.4, t1=0.3, seed=82, wob=0.05, pw=0.4)
    # the column of smoke: one ribbon, four bends (the stages), widening and fading
    pts = [(18.4, 40.0), (16.2, 34.0), (21.5, 28.5), (19.0, 22.5), (26.0, 17.0), (24.5, 11.0), (33.0, 6.5), (43.0, 8.5), (51.0, 4.5), (60.0, 5.5)]
    C = catmull(pts)
    P, tot = resample(C, 0.3)
    n = len(P)
    cuts = [0, int(n * 0.24), int(n * 0.46), int(n * 0.70), n]
    for k in range(4):
        seg = P[max(0, cuts[k] - 2):cuts[k + 1]]
        role = "ink" if k < 2 else "ink55"
        w = (1.3 + 0.75 * k) * ws
        stroke(ink, role, seg, w, t0=0.2, t1=0.3, seed=90 + k, wob=0.25, belly=0.35, smooth=False, lat=0.1)
    q = offset_path(pts, lambda t: 1.2 + 3.2 * t)
    stroke(ink, "ink55", q[int(len(q) * 0.2):int(len(q) * 0.95)], 0.7 * ws, t0=0.2, t1=0.6, seed=99, wob=0.4, smooth=False, lat=0.2)
    q = offset_path(pts, lambda t: -1.0 - 2.4 * t)
    stroke(ink, "ink55", q[int(len(q) * 0.35):int(len(q) * 0.8)], 0.5 * ws, t0=0.3, t1=0.6, seed=98, wob=0.4, smooth=False, lat=0.2)
    return ink

# =============================================================== E · the road
def enkin(ws=1.0):
    """E2.1 Perspective: the road narrows to the sun on the horizon."""
    ink = Ink()
    hz = 25.0
    vx = 40.0
    # far mountains and the horizon
    stroke(ink, "ink55", [(1.0, hz - 0.2), (63.0, hz + 0.2)], 0.8 * ws, t0=0.1, t1=0.1, seed=100, wob=0.2)
    kasure(ink, "ink55", [(2.0, hz - 1.0), (10.0, hz - 7.5), (16.0, hz - 4.5), (23.0, hz - 10.0), (31.0, hz - 1.5)], 2.0 * ws, seed=101, dry=0.4, strands=5)
    # the sun, half risen at the vanishing point: t^{n+1}
    blob(ink, "acc", vx, hz - 5.6, 6.0, 6.0, seed=102, rough=0.02)
    # the road: a pale wash between two dry edges, converging on the sun
    from brush import poly_d
    L = catmull([(-1.0, 64.0), (20.0, 45.0), (vx - 0.7, hz + 1.0)])
    R = catmull([(46.0, 64.0), (43.5, 45.0), (vx + 0.7, hz + 1.0)])
    ink.add("ink55", poly_d(L + R[::-1]))
    for x, y in L[::6] + R[::6]:
        ink.pts.append((x, y, 0.6))
    kasure(ink, "ink", [(-1.0, 64.0), (20.0, 45.0), (vx - 0.7, hz + 1.0)], 3.4 * ws, seed=103, dry=0.4, strands=6, wprof=lambda t: 1.0 - 0.75 * t)
    kasure(ink, "ink", [(46.0, 64.0), (43.5, 45.0), (vx + 0.7, hz + 1.0)], 2.4 * ws, seed=104, dry=0.3, strands=5, wprof=lambda t: 1.0 - 0.75 * t)
    # ruts: two pale lines down the middle
    stroke(ink, "ground", [(18.0, 64.0), (30.0, 45.0), (vx - 0.1, hz + 2.0)], 0.7 * ws, t0=0.0, t1=0.8, seed=105, wob=0.3, track=False, belly=0.0)
    stroke(ink, "ground", [(28.0, 64.0), (36.0, 45.0), (vx + 0.2, hz + 2.0)], 0.6 * ws, t0=0.0, t1=0.8, seed=106, wob=0.3, track=False, belly=0.0)
    # inns at c = 0, 1/2, 1/2, 1 (perspective depth)
    def depth(c):
        return hz + 1.5 + 36.0 / (1.0 + 2.6 * c) - 36.0 / 3.6 * c * 0.0
    y0, yh, y1 = depth(0.0), depth(0.5), depth(1.0)
    roof(ink, "ink", 9.0, 49.5, 2.3, 110, ws)
    roof(ink, "ink", 20.5, yh - 1.0, 1.35, 112, ws)
    roof(ink, "ink", 47.0, yh - 1.2, 1.35, 114, ws)
    roof(ink, "ink", 44.6, y1 - 1.6, 0.75, 116, ws)
    # a pine in the foreground
    stroke(ink, "ink", [(56.0, 48.0), (57.0, 34.0)], 1.2 * ws, t0=0.05, t1=0.3, seed=120, wob=0.1)
    for k, (yy, l) in enumerate([(35.0, 5.0), (39.0, 7.0), (43.0, 9.0)]):
        kasure(ink, "ink", [(56.5 - l / 2, yy + 0.6), (56.6, yy - 0.8), (56.5 + l / 2, yy + 0.7)], 1.9 * ws, seed=121 + k, dry=0.5, strands=4)
    return ink

def sugoroku(ws=1.0):
    """E2.2 The Tōkaidō board game: a spiral of stations, Kyoto in the centre."""
    ink = Ink()
    path = [(7.0, 49.0), (7.0, 7.0), (49.0, 7.0), (49.0, 49.0), (15.5, 49.0), (15.5, 15.5), (40.5, 15.5), (40.5, 40.5), (24.0, 40.5), (24.0, 24.0), (32.0, 24.0), (32.0, 30.0)]
    # straight segments with brushed corners
    P = []
    for i in range(len(path) - 1):
        seg, _ = resample([path[i], path[i + 1]], 0.3)
        P += seg if i == 0 else seg[1:]
    n = len(P)
    stroke(ink, "ink", P, 5.2 * ws, t0=0.004, t1=0.004, seed=130, wob=0.05, belly=0.0, smooth=False, minw=0.85, cap=0.85, lat=0.01)
    # cell cuts in ground colour: every station
    L, cut, k = 0.0, 4.6, 0
    nxt = cut
    for i in range(1, n - 6):
        L += math.hypot(P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1])
        if L >= nxt:
            nxt += cut
            tx, ty = P[i + 1][0] - P[i - 1][0], P[i + 1][1] - P[i - 1][1]
            l = math.hypot(tx, ty) or 1e-9
            nx, ny = -ty / l, tx / l
            x, y = P[i]
            # skip cuts right on a corner
            if any(math.hypot(x - cx, y - cy) < 2.2 for cx, cy in path):
                continue
            k += 1
            stroke(ink, "ground", [(x - nx * 3.3, y - ny * 3.3), (x + nx * 3.3, y + ny * 3.3)], 0.62 * ws, t0=0.05, t1=0.05, seed=140 + k, wob=0.3, track=False, belly=0.4)
    # agari: Kyoto in colour at the centre
    blob(ink, "acc", 29.0, 31.0, 4.3, 4.0, rot=0.1, seed=150, rough=0.06)
    # the start: a small roof at Nihonbashi
    roof(ink, "ink", 7.0, 57.0, 1.7, 152, ws)
    return ink

def tsuzuraori(ws=1.0):
    """E2.3 Hakone pass: the road winds up the mountain, the sun behind the top."""
    ink = Ink()
    blob(ink, "acc", 41.5, 13.0, 7.0, 7.0, seed=160, rough=0.02)
    # mountain mass: kasure strokes building a dark slope
    for k in range(7):
        y0 = 54.0 - k * 0.2
        kasure(ink, "ink", [(1.0 + k * 3.6, 60.0 - k * 0.3), (13.0 + k * 2.2, 38.0 - k * 1.5), (29.0 + k * 0.5, 19.0 + k * 1.3)],
               4.0 * ws, seed=170 + k, dry=0.55 - 0.03 * k, strands=6, t0=0.05)
    kasure(ink, "ink", [(29.5, 19.5), (38.0, 28.0), (44.0, 38.0), (51.0, 48.0)], 3.6 * ws, seed=178, dry=0.4, strands=6)
    # the switchback road, cut in ground colour, inns at the turns
    road = [(9.0, 60.0), (27.0, 51.5), (16.0, 45.0), (30.5, 38.0), (22.0, 31.5), (32.0, 25.0), (28.0, 19.0)]
    C = catmull(road)
    P, _ = resample(C, 0.3)
    stroke(ink, "ground", P, 1.5 * ws, t0=0.0, t1=0.6, seed=179, wob=0.2, smooth=False, track=False, belly=0.0, minw=0.3)
    roof(ink, "ink", 32.5, 49.5, 1.25, 180, ws)
    roof(ink, "ink", 35.5, 35.0, 1.1, 182, ws)
    roof(ink, "ink", 36.0, 34.0 - 0.0, 0.0, 184, ws) if False else None
    roof(ink, "ink", 13.0, 42.0, 1.1, 184, ws)
    roof(ink, "ink", 35.8, 22.5, 0.85, 186, ws)
    # pale far ridge on the right
    stroke(ink, "ink55", [(46.0, 32.0), (52.0, 27.0), (58.0, 30.0), (63.5, 26.5)], 1.0 * ws, t0=0.2, t1=0.3, seed=190, wob=0.25)
    return ink

def konomichi(ws=1.0):
    """E2.4 Incense road: the Tōkaidō laid in incense, burnt to the midpoint inns."""
    ink = Ink()
    road = [(3.0, 60.5), (24.0, 55.5), (37.0, 49.0), (30.0, 42.0), (17.0, 37.0), (22.0, 30.5), (36.0, 27.0), (44.0, 21.0), (41.5, 15.0), (36.0, 11.5)]
    C = catmull(road)
    P, tot = resample(C, 0.3)
    n = len(P)
    wprof = lambda t: 2.8 - 1.9 * t
    def seg(i0, i1, role, seed, **kw):
        part = P[i0:i1]
        if len(part) < 3:
            return
        t = (i0 + i1) / 2 / (n - 1)
        stroke(ink, role, part, wprof(t) * ws, seed=seed, smooth=False, **kw)
    b = int(n * 0.5)
    # ash behind
    k, s = 0, 200
    while k < b - 8:
        j = min(b - 3, k + 10)
        seg(k, j + 1, "ink55", s, t0=0.2, t1=0.3, wob=0.3, belly=0.0, minw=0.3); s += 1
        k = j + 4
    seg(b - 3, b + 12, "acc", 230, t0=0.5, t1=0.15, wob=0.08, pw=0.4)
    seg(b + 12, n, "ink", 231, t0=0.01, t1=0.3, wob=0.1, belly=0.0, minw=0.4)
    ex, ey = P[b + 1]
    smoke(ink, "ink55", [(ex + 0.5, ey - 1.2), (ex - 2.5, ey - 5.0), (ex + 1.0, ey - 8.5), (ex - 2.0, ey - 12.5)], 1.0 * ws, seed=232)
    # inns along the road: c = 0, 1/2, 1/2, 1
    def at(f):
        return P[min(n - 1, int(f * (n - 1)))]
    x, y = at(0.02); roof(ink, "ink", x + 4.0, y - 5.5, 1.9, 240, ws)
    x, y = at(0.5); roof(ink, "ink", x - 8.5, y + 2.5, 1.45, 242, ws); roof(ink, "ink", x - 4.0, y + 4.6, 1.45, 244, ws)
    x, y = at(0.99); roof(ink, "ink", x - 4.6, y - 0.4, 1.0, 246, ws)
    # far hills
    stroke(ink, "ink55", [(42.0, 8.0), (50.0, 3.5), (56.0, 6.0), (63.5, 2.5)], 0.9 * ws, t0=0.2, t1=0.3, seed=250, wob=0.25)
    return ink

MARKS2 = {"senkodai": senkodai, "uzumaki": uzumaki, "koban": koban, "kemuri": kemuri,
          "enkin": enkin, "sugoroku": sugoroku, "tsuzuraori": tsuzuraori, "konomichi": konomichi}
