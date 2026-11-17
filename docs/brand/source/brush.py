"""Sumi brush engine for the SELAH round-1 marks (family conventions: 64-unit
mark box, brushed 3-level AMR ground refined under ink, carved seal)."""
import math, random

def F(v):
    return f"{v:.2f}"

# ---------------------------------------------------------------- geometry
def catmull(pts, per=12):
    if len(pts) < 3:
        out = []
        (x0, y0), (x1, y1) = pts[0], pts[-1]
        n = max(2, int(math.hypot(x1 - x0, y1 - y0) / 0.35))
        return [(x0 + (x1 - x0) * i / n, y0 + (y1 - y0) * i / n) for i in range(n + 1)]
    P = [pts[0]] + list(pts) + [pts[-1]]
    out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        seg = max(per, int(math.hypot(p2[0] - p1[0], p2[1] - p1[1]) / 0.3))
        for k in range(seg):
            t = k / seg
            t2, t3 = t * t, t * t * t
            x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3)
            y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3)
            out.append((x, y))
    out.append(pts[-1])
    return out

def resample(pts, step):
    L = [0.0]
    for i in range(1, len(pts)):
        L.append(L[-1] + math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]))
    tot = L[-1]
    n = max(3, int(tot / step))
    out, j = [], 0
    for k in range(n + 1):
        s = tot * k / n
        while j < len(L) - 2 and L[j + 1] < s:
            j += 1
        seg = L[j + 1] - L[j] or 1e-9
        u = (s - L[j]) / seg
        out.append((pts[j][0] + (pts[j + 1][0] - pts[j][0]) * u, pts[j][1] + (pts[j + 1][1] - pts[j][1]) * u))
    return out, tot

def arcpts(cx, cy, r, a0, a1, n=None, rj=0.0, seed=0):
    rnd = random.Random(seed)
    ph = rnd.uniform(0, 6.28)
    n = n or max(8, int(abs(a1 - a0) * r / 1.0))
    return [(cx + (r + rj * math.sin(3 * (a0 + (a1 - a0) * i / n) + ph)) * math.cos(a0 + (a1 - a0) * i / n),
             cy + (r + rj * math.sin(3 * (a0 + (a1 - a0) * i / n) + ph)) * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]

# ---------------------------------------------------------------- strokes
class Ink:
    """Collects polygons by layer role and centre-line samples for the AMR."""
    def __init__(self):
        self.layers = {}   # role -> list of d
        self.pts = []      # (x, y, r) ink samples
        self.order = []

    def add(self, role, d, evenodd=False):
        key = (role, evenodd)
        if key not in self.layers:
            self.layers[key] = []
            self.order.append(key)
        self.layers[key].append(d)

def poly_d(poly):
    return "M" + F(poly[0][0]) + " " + F(poly[0][1]) + "L" + " ".join(F(x) + " " + F(y) for x, y in poly[1:]) + "Z"

def stroke(ink, role, pts, w, t0=0.18, t1=0.35, seed=0, wob=0.10, minw=0.10,
           belly=0.18, smooth=True, track=True, lat=0.04, step=0.35, pw=0.7, cap=0.0):
    rnd = random.Random(seed)
    ph = [rnd.uniform(0, 6.283) for _ in range(6)]
    C = catmull(pts) if smooth and len(pts) > 2 else catmull(pts[:1] + pts[-1:]) if len(pts) == 2 else pts
    P, tot = resample(C, step)
    n = len(P)
    left, right = [], []
    for i, (x, y) in enumerate(P):
        t = i / (n - 1)
        a = min(1.0, t / t0) if t0 > 0 else 1.0
        b = min(1.0, (1 - t) / t1) if t1 > 0 else 1.0
        f = max(minw, (a ** pw) * (b ** 0.9)) if (t0 > 0 or t1 > 0) else 1.0
        f = max(f, cap) if (t < 0.02 or t > 0.98) else f
        f *= 1 + belly * math.sin(math.pi * t)
        noise = 1 + wob * (0.6 * math.sin(2 * math.pi * 2.3 * t + ph[0]) + 0.4 * math.sin(2 * math.pi * 6.1 * t + ph[1]))
        hw = 0.5 * w * f * noise
        i0, i1 = max(0, i - 1), min(n - 1, i + 1)
        tx, ty = P[i1][0] - P[i0][0], P[i1][1] - P[i0][1]
        l = math.hypot(tx, ty) or 1e-9
        nx, ny = -ty / l, tx / l
        off = lat * w * math.sin(2 * math.pi * 1.7 * t + ph[2])
        cx, cy = x + nx * off, y + ny * off
        left.append((cx + nx * hw, cy + ny * hw))
        right.append((cx - nx * hw, cy - ny * hw))
        if track:
            ink.pts.append((x, y, hw))
    ink.add(role, poly_d(left + right[::-1]))

def blob(ink, role, cx, cy, rx, ry, rot=0.0, seed=0, rough=0.06, n=36, track=True):
    rnd = random.Random(seed)
    ph = [rnd.uniform(0, 6.283) for _ in range(4)]
    poly = []
    c, s = math.cos(rot), math.sin(rot)
    for k in range(n):
        a = 2 * math.pi * k / n
        r = 1 + rough * (0.6 * math.sin(2 * a + ph[0]) + 0.3 * math.sin(5 * a + ph[1]) + 0.2 * math.sin(9 * a + ph[2]))
        x, y = rx * r * math.cos(a), ry * r * math.sin(a)
        poly.append((cx + x * c - y * s, cy + x * s + y * c))
    if track:
        for k in range(0, n, 3):
            ink.pts.append((poly[k][0], poly[k][1], 0.3))
        ink.pts.append((cx, cy, max(rx, ry)))
    ink.add(role, poly_d(poly))

def drop(ink, role, cx, cy, a, seed=0, track=True):
    poly = []
    for k in range(40):
        t = 2 * math.pi * k / 40
        x = a * 0.62 * math.sin(t) * math.sin(t / 2) ** 1.1
        y = -a * math.cos(t)
        poly.append((cx + x, cy + y * 0.95 + a * 0.2))
    if track:
        ink.pts.append((cx, cy, a))
    ink.add(role, poly_d(poly))

# ---------------------------------------------------------------- family ground & seal
LEVELS = [(16, 0.30, 0.62), (8, 0.22, 0.384), (4, 0.15, 0.186)]

def meshline(out, x0, y0, x1, y1, w, seed):
    tmp = Ink()
    stroke(tmp, "m", [(x0, y0), (x1, y1)], w, t0=0.01, t1=0.01, seed=seed, wob=0.18, minw=0.7,
           belly=0.0, track=False, lat=0.25, step=0.9)
    out.append(tmp.layers[("m", False)][0])

def hit(pts, x0, y0, s):
    for x, y, r in pts:
        if x0 - r < x < x0 + s + r and y0 - r < y < y0 + s + r:
            return True
    return False

def amr(ink, levels=3, seed=7):
    """Balanced quadtree: level 0 the 16-grid, refined one level per stroke."""
    out = [[] for _ in range(3)]
    sd = seed
    for k in (16, 32, 48):
        sd += 1; meshline(out[0], -0.4, k, 64.4, k + 0.05, LEVELS[0][1], sd)
        sd += 1; meshline(out[0], k, -0.4, k - 0.05, 64.4, LEVELS[0][1], sd)
    if levels > 1:
        for i in range(4):
            for j in range(4):
                x0, y0 = 16 * i, 16 * j
                if hit(ink.pts, x0, y0, 16):
                    sd += 1; meshline(out[1], x0 + 8, y0, x0 + 8, y0 + 16, LEVELS[1][1], sd)
                    sd += 1; meshline(out[1], x0, y0 + 8, x0 + 16, y0 + 8, LEVELS[1][1], sd)
                    if levels > 2:
                        for a in range(2):
                            for b in range(2):
                                u0, v0 = x0 + 8 * a, y0 + 8 * b
                                if hit(ink.pts, u0, v0, 8):
                                    sd += 1; meshline(out[2], u0 + 4, v0, u0 + 4, v0 + 8, LEVELS[2][1], sd)
                                    sd += 1; meshline(out[2], u0, v0 + 4, u0 + 8, v0 + 4, LEVELS[2][1], sd)
    return [(" ".join(o), LEVELS[k][2]) for k, o in enumerate(out) if o]

def seal_d(x0=52.6, y0=52.6, s=9.0, seed=3):
    rnd = random.Random(seed)
    def jit(v): return v + rnd.uniform(-0.12, 0.12)
    x1, y1 = x0 + s, y0 + s
    outer = []
    for k in range(6):
        outer.append((jit(x0 + 0.5 + (s - 1) * k / 5), jit(y0)))
    outer.append((x1 - 0.2, y0 + 0.3))
    for k in range(6):
        outer.append((jit(x1), jit(y0 + 0.5 + (s - 1) * k / 5)))
    outer.append((x1 - 0.35, y1 - 0.25))
    for k in range(6):
        outer.append((jit(x1 - 0.5 - (s - 1) * k / 5), jit(y1)))
    outer.append((x0 + 0.3, y1 - 0.3))
    for k in range(6):
        outer.append((jit(x0), jit(y1 - 0.5 - (s - 1) * k / 5)))
    d = poly_d(outer)
    # carved quadrants, top-right one refined (the family seal)
    m, bar = 1.55, 0.85
    ix0, iy0, ix1, iy1 = x0 + m, y0 + m, x1 - m, y1 - m
    cx, cy = (ix0 + ix1) / 2, (iy0 + iy1) / 2
    def hole(a, b, c, e):
        return poly_d([(jit(a), jit(b)), (jit(c), jit(b)), (jit(c), jit(e)), (jit(a), jit(e))])
    h = bar / 2
    d += " " + hole(ix0, iy0, cx - h, cy - h)
    d += " " + hole(ix0, cy + h, cx - h, iy1)
    d += " " + hole(cx + h, cy + h, ix1, iy1)
    qx0, qy0, qx1, qy1 = cx + h, iy0, ix1, cy - h
    mx, my = (qx0 + qx1) / 2, (qy0 + qy1) / 2
    hh = 0.28
    d += " " + hole(qx0, qy0, mx - hh, my - hh) + " " + hole(mx + hh, qy0, qx1, my - hh)
    d += " " + hole(qx0, my + hh, mx - hh, qy1) + " " + hole(mx + hh, my + hh, qx1, qy1)
    return d

def offset_path(pts, off_fn, step=0.35):
    C = catmull(pts) if len(pts) > 2 else catmull(pts[:1] + pts[-1:])
    P, tot = resample(C, step)
    n = len(P)
    out = []
    for i, (x, y) in enumerate(P):
        i0, i1 = max(0, i - 1), min(n - 1, i + 1)
        tx, ty = P[i1][0] - P[i0][0], P[i1][1] - P[i0][1]
        l = math.hypot(tx, ty) or 1e-9
        nx, ny = -ty / l, tx / l
        o = off_fn(i / (n - 1))
        out.append((x + nx * o, y + ny * o))
    return out

def kasure(ink, role, pts, w, seed=0, dry=0.55, strands=7, t0=0.12, wob=0.1, gaps=True, wprof=None):
    """Dry brush: a loaded core that splits into bristle streaks towards the tail."""
    rnd = random.Random(seed)
    wprof = wprof or (lambda t: 1.0)
    C = catmull(pts) if len(pts) > 2 else catmull(pts[:1] + pts[-1:])
    P, tot = resample(C, 0.3)
    n = len(P)
    core_end = int(n * dry)
    if core_end > 3:
        stroke(ink, role, P[:core_end + 2], w * wprof(dry * 0.5), t0=t0 / max(dry, 0.05), t1=0.08, seed=seed, wob=wob, smooth=False, minw=0.55, belly=0.08)
    for s in range(strands):
        o = (s / (strands - 1) - 0.5) * w * 0.82
        start = max(0, core_end - int(n * 0.06)) if core_end > 3 else 0
        end = int(n * rnd.uniform(dry + (1 - dry) * 0.45, 1.0))
        if end - start < 4:
            continue
        line = []
        for i in range(start, end):
            i0, i1 = max(0, i - 1), min(n - 1, i + 1)
            tx, ty = P[i1][0] - P[i0][0], P[i1][1] - P[i0][1]
            l = math.hypot(tx, ty) or 1e-9
            nx, ny = -ty / l, tx / l
            t = i / (n - 1)
            oo = o * wprof(t) * (1 + 0.15 * math.sin(9 * t + s))
            line.append((P[i][0] + nx * oo, P[i][1] + ny * oo))
        # break the streak once or twice
        pieces = [line]
        if gaps and len(line) > 20 and rnd.random() < 0.7:
            c = rnd.randint(int(len(line) * 0.4), int(len(line) * 0.85))
            g = rnd.randint(2, 5)
            pieces = [line[:c], line[c + g:]]
        for k, pc in enumerate(pieces):
            if len(pc) < 3:
                continue
            sw = w / strands * rnd.uniform(1.5, 2.3) * wprof(start / (n - 1))
            stroke(ink, role, pc, sw, t0=0.02 if k == 0 else 0.2, t1=0.5, seed=seed * 31 + s * 7 + k, wob=0.2, smooth=False, minw=0.15, belly=0.0, lat=0.02)

def smoke(ink, role, pts, w, seed=0, strands=2):
    stroke(ink, role, pts, w, t0=0.12, t1=0.75, seed=seed, wob=0.3, belly=0.35, pw=0.6, lat=0.15)
    for s in range(1, strands):
        q = offset_path(pts, lambda t, s=s: (0.6 + 1.4 * t) * (1 if s % 2 else -1) * w * 0.9)
        stroke(ink, role, q[: int(len(q) * 0.85)], w * 0.45, t0=0.25, t1=0.7, seed=seed + 50 * s, wob=0.35, smooth=False, lat=0.2)

# ---------------------------------------------------------------- sumi AMR ground (round 3)
SUMI_LEVELS = [
    # spacing, width, opacity, entry press, bow, break prob, jitter
    dict(sp=16, w=0.46, op=0.58, t0=0.05, t1=0.30, pw=0.35, wob=0.28, lat=0.10, brk=0.25, jit=0.15),
    dict(sp=8,  w=0.30, op=0.40, t0=0.08, t1=0.40, pw=0.45, wob=0.32, lat=0.16, brk=0.40, jit=0.20),
    dict(sp=4,  w=0.19, op=0.26, t0=0.12, t1=0.50, pw=0.55, wob=0.40, lat=0.22, brk=0.55, jit=0.22),
]

def sumi_seg(out, p, q, L, rnd):
    tmp = Ink()
    (x0, y0), (x1, y1) = p, q
    j = L["jit"]
    pts = [(x0 + rnd.uniform(-j, j), y0 + rnd.uniform(-j, j)), (x1 + rnd.uniform(-j, j), y1 + rnd.uniform(-j, j))]
    # bow: a gentle mid-point offset, as a hand-drawn ruled line
    mx, my = (pts[0][0] + pts[1][0]) / 2, (pts[0][1] + pts[1][1]) / 2
    dx, dy = pts[1][0] - pts[0][0], pts[1][1] - pts[0][1]
    l = math.hypot(dx, dy) or 1
    b = rnd.uniform(-1, 1) * L["lat"] * 1.4
    mid = (mx - dy / l * b, my + dx / l * b)
    C = catmull([pts[0], mid, pts[1]])
    P, tot = resample(C, 0.5)
    pieces = [P]
    if rnd.random() < L["brk"] and len(P) > 12:
        c = rnd.randint(int(len(P) * 0.45), int(len(P) * 0.8))
        g = rnd.randint(1, 3)
        pieces = [P[:c], P[c + g:]]
    for k, pc in enumerate(pieces):
        if len(pc) < 3:
            continue
        w = L["w"] * rnd.uniform(0.85, 1.2) * (1.0 if k == 0 else 0.75)
        stroke(tmp, "m", pc, w, t0=L["t0"] if k == 0 else 0.25, t1=L["t1"], seed=rnd.randint(0, 9999), wob=L["wob"],
               minw=0.18, belly=0.12, pw=L["pw"], smooth=False, track=False, lat=0.03, step=0.5)
    if ("m", False) in tmp.layers:
        out.extend(tmp.layers[("m", False)])

def amr_sumi(ink, seed=7):
    rnd = random.Random(seed)
    out = [[], [], []]
    L0, L1, L2 = SUMI_LEVELS
    # level 0: the 16-grid, one stroke per cell edge, joints overlapping
    for k in (16, 32, 48):
        for c in range(4):
            a, b = 16 * c - rnd.uniform(0.2, 0.9), 16 * c + 16 + rnd.uniform(-0.2, 0.5)
            sumi_seg(out[0], (a, k), (b, k), L0, rnd)
            sumi_seg(out[0], (k, a), (k, b), L0, rnd)
    for i in range(4):
        for j in range(4):
            x0, y0 = 16 * i, 16 * j
            if not hit(ink.pts, x0, y0, 16):
                continue
            o = rnd.uniform(-0.6, 0.2)
            sumi_seg(out[1], (x0 + 8, y0 + o), (x0 + 8, y0 + 16 + rnd.uniform(-0.3, 0.4)), L1, rnd)
            sumi_seg(out[1], (x0 + o, y0 + 8), (x0 + 16 + rnd.uniform(-0.3, 0.4), y0 + 8), L1, rnd)
            for a in range(2):
                for b in range(2):
                    u0, v0 = x0 + 8 * a, y0 + 8 * b
                    if hit(ink.pts, u0, v0, 8):
                        o = rnd.uniform(-0.4, 0.2)
                        sumi_seg(out[2], (u0 + 4, v0 + o), (u0 + 4, v0 + 8 + rnd.uniform(-0.3, 0.3)), L2, rnd)
                        sumi_seg(out[2], (u0 + o, v0 + 4), (u0 + 8 + rnd.uniform(-0.3, 0.3), v0 + 4), L2, rnd)
    return [(" ".join(o), SUMI_LEVELS[k]["op"]) for k, o in enumerate(out) if o]
