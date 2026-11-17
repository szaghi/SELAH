import math, random
from brush import Ink, stroke, blob, drop, arcpts, resample, catmull, kasure, smoke, offset_path, poly_d

def earc(cx, cy, rx, ry, a0, a1, n=None, rot=0.0):
    n = n or max(10, int(abs(a1 - a0) * max(rx, ry) / 0.8))
    c, s = math.cos(rot), math.sin(rot)
    out = []
    for i in range(n + 1):
        a = a0 + (a1 - a0) * i / n
        x, y = rx * math.cos(a), ry * math.sin(a)
        out.append((cx + x * c - y * s, cy + x * s + y * c))
    return out

def S(ink, role, pts, w, seed, **kw):
    kw.setdefault("t0", 0.12); kw.setdefault("t1", 0.3); kw.setdefault("wob", 0.1)
    stroke(ink, role, pts, w, seed=seed, **kw)

# ---------------------------------------------------------------- 4.1 Ma 間
def ma(ws=1.0):
    """The calligrapher: 間, the sun seen through the gap of a gate. SELAH means a pause."""
    ink = Ink()
    w = 3.1 * ws
    # 門, written: left half
    S(ink, "ink", [(12.5, 7.0), (12.0, 32.0), (12.4, 57.0)], w * 1.05, 1, t0=0.06, t1=0.35, pw=0.4, belly=0.05)
    S(ink, "ink", [(13.0, 8.0), (25.5, 7.6), (26.2, 9.0)], w * 0.8, 2, t0=0.1, t1=0.15)
    S(ink, "ink", [(26.0, 8.5), (25.6, 24.5)], w * 0.8, 3, t0=0.1, t1=0.4)
    S(ink, "ink", [(13.5, 15.8), (24.8, 15.5)], w * 0.6, 4, t0=0.15, t1=0.3)
    S(ink, "ink", [(13.5, 23.6), (24.8, 23.4)], w * 0.6, 5, t0=0.15, t1=0.3)
    # right half
    S(ink, "ink", [(35.4, 8.5), (35.6, 24.5)], w * 0.8, 6, t0=0.1, t1=0.4)
    S(ink, "ink", [(35.6, 7.8), (49.5, 7.4), (50.4, 9.5), (50.2, 34.0), (49.6, 54.0), (47.0, 55.0), (44.6, 53.0)], w * 1.0, 7, t0=0.08, t1=0.12, pw=0.4, belly=0.05)
    S(ink, "ink", [(36.5, 15.6), (49.0, 15.4)], w * 0.6, 8, t0=0.15, t1=0.3)
    S(ink, "ink", [(36.5, 23.4), (49.0, 23.2)], w * 0.6, 9, t0=0.15, t1=0.3)
    # the sun in the gap
    blob(ink, "acc", 31.0, 40.0, 8.0, 8.0, seed=10, rough=0.03)
    return ink

# ---------------------------------------------------------------- 4.2 Toki no kane 時の鐘
def kane(ws=1.0):
    """The bell of time: struck to tell the hour; the sound spreads in rings."""
    ink = Ink()
    kasure(ink, "ink", [(4.0, 6.0), (32.0, 5.0), (60.0, 6.2)], 3.4 * ws, seed=20, dry=0.6, strands=6)
    S(ink, "ink", [(32.0, 6.5), (32.0, 10.0)], 1.4 * ws, 21, t0=0.05, t1=0.05)
    S(ink, "ink", [(29.0, 11.6), (32.0, 9.4), (35.0, 11.6)], 1.6 * ws, 22, t0=0.2, t1=0.2)
    # the body: dome, sides, lip
    kasure(ink, "ink", earc(32.0, 19.5, 11.5, 7.5, math.pi * 1.0, math.pi * 2.0), 2.6 * ws, seed=23, dry=0.7, strands=5)
    S(ink, "ink", [(20.6, 19.5), (20.0, 33.0), (18.8, 46.0)], 2.6 * ws, 24, t0=0.05, t1=0.1, belly=0.05)
    S(ink, "ink", [(43.4, 19.5), (44.0, 33.0), (45.2, 46.0)], 2.6 * ws, 25, t0=0.05, t1=0.1, belly=0.05)
    kasure(ink, "ink", [(16.5, 47.0), (32.0, 46.4), (47.5, 47.2)], 3.8 * ws, seed=26, dry=0.6, strands=6)
    # bands and bosses (pale)
    S(ink, "ink55", [(20.6, 26.0), (43.6, 25.8)], 0.9 * ws, 27, t0=0.1, t1=0.1)
    S(ink, "ink55", [(20.2, 38.5), (44.4, 38.3)], 0.9 * ws, 28, t0=0.1, t1=0.1)
    S(ink, "ink55", [(32.0, 26.0), (32.0, 38.5)], 0.8 * ws, 29, t0=0.1, t1=0.1)
    for i in range(3):
        for j in range(2):
            blob(ink, "ink55", 24.5 + 3.4 * i, 29.0 + 3.6 * j, 0.9, 0.9, seed=30 + 3 * i + j)
            blob(ink, "ink55", 35.5 + 3.4 * i, 29.0 + 3.6 * j, 0.9, 0.9, seed=40 + 3 * i + j)
    # the striker log
    kasure(ink, "ink", [(63.5, 40.0), (54.0, 41.0), (47.6, 41.6)], 3.2 * ws, seed=50, dry=0.5, strands=5)
    S(ink, "ink55", [(58.0, 40.5), (58.6, 6.5)], 0.6 * ws, 51, t0=0.05, t1=0.05)
    # sound: four rings to the left of the bell, the last in colour
    for k, r in enumerate([4.0, 7.2, 10.4, 13.6]):
        role = "acc" if k == 3 else ("ink" if k == 0 else "ink55")
        S(ink, role, arcpts(16.0, 36.0, r, math.radians(140), math.radians(220), seed=60 + k), (1.5 - 0.1 * k) * ws, 61 + k, t0=0.2, t1=0.4, belly=0.2)
    return ink

# ---------------------------------------------------------------- 4.3 Wadokei 和時計
def wadokei(ws=1.0):
    """The Edo clock: six day hours and six night hours of unequal length."""
    ink = Ink()
    cx, cy, R = 30.0, 30.0, 21.0
    S(ink, "ink", arcpts(cx, cy, R, math.radians(-80), math.radians(262), rj=0.3, seed=70), 2.8 * ws, 71, t0=0.06, t1=0.4, belly=0.2, pw=0.5)
    # night half pale (summer: short nights): a wash band on the inside
    night = (math.radians(30), math.radians(150))
    S(ink, "ink55", arcpts(cx, cy, R - 4.0, night[0], night[1], seed=72), 3.2 * ws, 73, t0=0.2, t1=0.3)
    # hour plates: day hours long, night hours short
    angs = []
    a = night[1]
    day = (2 * math.pi - (night[1] - night[0])) / 6
    nig = (night[1] - night[0]) / 6
    for k in range(6):
        angs.append(a); a += day
    a = night[0]
    for k in range(6):
        angs.append(a); a += nig
    for k, a in enumerate(angs):
        c, s = math.cos(a), math.sin(a)
        S(ink, "ink", [(cx + c * (R - 7.0), cy + s * (R - 7.0)), (cx + c * (R - 2.2), cy + s * (R - 2.2))], 1.9 * ws, 80 + k, t0=0.05, t1=0.1, belly=0.0)
    # the hand
    a = angs[3] + day * 0.4
    S(ink, "acc", [(cx - math.cos(a) * 3.0, cy - math.sin(a) * 3.0), (cx + math.cos(a) * (R - 8.0), cy + math.sin(a) * (R - 8.0))], 2.4 * ws, 95, t0=0.15, t1=0.7, pw=0.5)
    blob(ink, "ink", cx, cy, 2.0, 2.0, seed=96)
    return ink

# ---------------------------------------------------------------- 4.4 Tsuzumi 鼓
def tsuzumi(ws=1.0):
    """The hand drum: an hourglass laced head to head; struck, then the pause."""
    ink = Ink()
    # two heads, as brushed ellipses
    for k, y in enumerate([11.0, 49.0]):
        S(ink, "ink", earc(30.0, y, 18.5, 4.0, math.radians(-170), math.radians(185)), 2.4 * ws, 100 + k, t0=0.06, t1=0.3, belly=0.15)
    # the hourglass body between them
        body_l = [(20.0, 14.5), (26.4, 24.0), (27.8, 30.0), (26.4, 36.0), (20.0, 45.5)]
    body_r = [(40.0, 14.5), (33.6, 24.0), (32.2, 30.0), (33.6, 36.0), (40.0, 45.5)]
    ink.add("ink", poly_d(catmull(body_l) + catmull(body_r)[::-1]))
    for x, y in catmull(body_l)[::5] + catmull(body_r)[::5]:
        ink.pts.append((x, y, 1.0))
    # the lacing (shirabe), traditionally vermilion: head to head, crossing at the waist
    xs = [13.5, 25.0, 36.5]
    for k, x in enumerate(xs):
        S(ink, "acc", [(x, 13.5), (x + 11.0, 46.5)], 0.85 * ws, 110 + k, t0=0.05, t1=0.05, belly=0.0, wob=0.15)
        S(ink, "acc", [(x + 11.0, 13.5), (x, 46.5)], 0.85 * ws, 120 + k, t0=0.05, t1=0.05, belly=0.0, wob=0.15)
    # the binding cord round the waist
    S(ink, "acc", [(24.5, 30.6), (30.0, 32.2), (35.5, 30.5)], 2.0 * ws, 130, t0=0.15, t1=0.2)
    S(ink, "acc", [(35.0, 31.2), (42.0, 38.0), (44.0, 46.0), (42.5, 55.0)], 1.1 * ws, 131, t0=0.1, t1=0.6)
    return ink

# ---------------------------------------------------------------- 4.5 Hamon 波紋
def hamon(ws=1.0):
    """The ripple: one drop, four rings spreading on still water."""
    ink = Ink()
    drop(ink, "acc", 30.0, 10.0, 3.3 * ws)
    cy = 40.0
    for k, (rx, w, role) in enumerate([(4.5, 1.3, "ink"), (10.5, 1.15, "ink"), (17.0, 1.0, "ink55"), (24.0, 0.9, "ink55")]):
        ry = rx * 0.32
        a0 = math.radians(-160 + 15 * k)
        stroke(ink, role, earc(30.0, cy, rx, ry, a0, a0 + math.radians(330)), w * ws, t0=0.06, t1=0.45, seed=140 + k, wob=0.15, belly=0.2, pw=0.5)
    # splash: two short strokes up from the centre
    blob(ink, "ink", 30.0, 39.6, 1.2, 0.7, seed=150)
    # the drop's fall, as a faint dry trace
    S(ink, "ink55", [(30.0, 15.5), (30.1, 30.0)], 0.5 * ws, 151, t0=0.3, t1=0.3, wob=0.3)
    return ink

# ---------------------------------------------------------------- 4.6 Soroban 算盤
def soroban(ws=1.0):
    """The abacus of the wasan mathematicians: RK4's weights set as 1 2 2 1 over 6."""
    ink = Ink()
    x0, x1, y0, y1, beam = 5.0, 51.0, 11.0, 49.5, 20.5
    kasure(ink, "ink", [(x0, y0), (x1, y0 - 0.4)], 3.2 * ws, seed=160, dry=0.7, strands=6)
    kasure(ink, "ink", [(x0 - 0.2, y1), (x1, y1 + 0.3)], 3.2 * ws, seed=161, dry=0.7, strands=6)
    S(ink, "ink", [(x0, y0 - 1.0), (x0 + 0.2, y1 + 1.0)], 2.6 * ws, 162, t0=0.05, t1=0.05, belly=0.0)
    S(ink, "ink", [(x1, y0 - 1.0), (x1 - 0.2, y1 + 1.0)], 2.6 * ws, 163, t0=0.05, t1=0.05, belly=0.0)
    S(ink, "ink", [(x0, beam), (x1, beam - 0.2)], 1.8 * ws, 164, t0=0.05, t1=0.1, belly=0.0)
    xs = [11.5, 19.5, 27.5, 35.5, 43.8]
    counts = [1, 2, 2, 1, 6]
    for k, (x, v) in enumerate(zip(xs, counts)):
        S(ink, "ink55", [(x, y0 + 1.0), (x, y1 - 1.0)], 0.6 * ws, 170 + k, t0=0.02, t1=0.02)
        role = "acc" if k == 4 else "ink"
        five, ones = divmod(v, 5)
        def bead(y, sd):
            pts = [(x - 3.3, y), (x - 0.4, y - 1.75), (x + 0.4, y - 1.75), (x + 3.3, y), (x + 0.4, y + 1.75), (x - 0.4, y + 1.75)]
            ink.add(role, poly_d(pts))
            ink.pts.append((x, y, 3.0))
        bead(beam - 2.4 if five else y0 + 3.4, 180 + k)
        for i in range(ones):
            bead(beam + 2.4 + 3.6 * i, 190 + 10 * k + i)
        for j in range(4 - ones):
            bead(y1 - 3.0 - 3.6 * j, 195 + 10 * k + j)
    return ink

# ---------------------------------------------------------------- 4.7 Ume 梅
def blossom(ink, role, x, y, r, seed, petals=5, open_=1.0):
    rnd = random.Random(seed)
    a0 = rnd.uniform(0, 1.2)
    for k in range(petals):
        a = a0 + 2 * math.pi * k / 5
        blob(ink, role, x + math.cos(a) * r * 0.62 * open_, y + math.sin(a) * r * 0.62 * open_, r * 0.46, r * 0.46, seed=seed + k, rough=0.1)

def ume(ws=1.0):
    """The plum branch: bud, half-open, full, falling, four stages of one flower."""
    ink = Ink()
    kasure(ink, "ink", [(1.0, 60.0), (10.0, 46.0), (16.0, 38.0), (26.0, 30.5), (34.0, 21.0), (46.0, 14.0), (58.0, 5.0)], 4.2 * ws, seed=200, dry=0.55, strands=7,
           wprof=lambda t: 1.0 - 0.55 * t)
    S(ink, "ink", [(16.5, 38.0), (22.0, 44.0), (30.5, 46.0)], 1.5 * ws, 201, t0=0.05, t1=0.6)
    S(ink, "ink", [(34.0, 21.5), (34.5, 13.0), (30.0, 6.0)], 1.3 * ws, 202, t0=0.05, t1=0.6)
    # stage 1: a bud
    blob(ink, "ink", 30.6, 6.6, 1.4, 1.9, rot=0.4, seed=210)
    # stage 2: half open
    blossom(ink, "ink55", 31.0, 46.2, 3.2, 220, petals=3, open_=0.7)
    blob(ink, "ink", 31.0, 46.2, 0.8, 0.8, seed=225)
    # stage 3: full, in colour
    blossom(ink, "acc", 44.0, 18.0, 6.0, 230)
    for k in range(5):
        a = 0.3 + 2 * math.pi * k / 5
        S(ink, "ink", [(44.0, 18.0), (44.0 + math.cos(a) * 2.6, 18.0 + math.sin(a) * 2.6)], 0.5 * ws, 236 + k, t0=0.1, t1=0.2)
    # stage 4: a petal falling
    blob(ink, "ink55", 47.0, 40.0, 1.9, 1.3, rot=0.7, seed=245)
    blob(ink, "ink55", 41.0, 48.5, 1.5, 1.0, rot=-0.3, seed=246)
    return ink

# ---------------------------------------------------------------- 4.8 Hyōshigi 拍子木
def hyoshigi(ws=1.0):
    """Kabuki clappers: two blocks, a cord, and the pause between two claps."""
    ink = Ink()
    # two blocks, a narrow gap between them (ma)
    for k, y in enumerate([24.5, 35.5]):
        stroke(ink, "ink", [(9.0, y), (30.0, y - 0.3), (47.0, y + 0.2)], 7.6 * ws, t0=0.015, t1=0.015, seed=301 + k, wob=0.04, belly=0.0, minw=0.92, cap=0.92)
        kasure(ink, "ink", [(30.0, y - 0.3), (47.0, y + 0.2), (48.6, y + 0.2)], 7.6 * ws, seed=303 + k, dry=0.05, strands=8, gaps=False)
    # the cord through both blocks, a loop and a tassel
    S(ink, "ink", [(10.5, 22.0), (5.0, 26.5), (4.0, 30.0), (5.0, 33.5), (10.5, 38.0)], 1.1 * ws, 305, t0=0.05, t1=0.05)
    S(ink, "ink", [(4.4, 31.0), (3.6, 40.0), (5.0, 50.0)], 1.0 * ws, 306, t0=0.05, t1=0.7)
    # the clap at the far ends: four short strokes, the fourth in colour
    cx, cy = 51.0, 30.0
    for k, a in enumerate([-60, -25, 25, 60]):
        aa = math.radians(a)
        role = "acc" if k == 3 else "ink"
        S(ink, role, [(cx + math.cos(aa) * 3.0, cy + math.sin(aa) * 3.0), (cx + math.cos(aa) * 9.5, cy + math.sin(aa) * 9.5)], 1.5 * ws, 310 + k, t0=0.45, t1=0.3)
    return ink

MARKS4 = {"ma": ma, "kane": kane, "wadokei": wadokei, "tsuzumi": tsuzumi, "hamon": hamon, "soroban": soroban, "ume": ume, "hyoshigi": hyoshigi}
