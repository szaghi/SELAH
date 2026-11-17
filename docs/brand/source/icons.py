import math
from brush import Ink, stroke, blob, drop, arcpts

def element(ink, style, x, y, s, role, seed, ws=1.0):
    """The stage glyph of each hand, centred on (x, y), size s."""
    if style == "arc":
        stroke(ink, role, arcpts(x, y + 0.6 * s, 1.15 * s, math.radians(200), math.radians(330), seed=seed), 0.85 * s * ws, t0=0.15, t1=0.5, seed=seed, wob=0.1, pw=0.5)
    elif style == "ember":
        blob(ink, role, x, y, 0.62 * s, 0.62 * s, seed=seed, rough=0.08)
    elif style == "drop":
        drop(ink, role, x, y, 0.85 * s)
    elif style == "leaf":
        stroke(ink, role, [(x - 1.1 * s, y + 0.55 * s), (x, y - 0.15 * s), (x + 1.2 * s, y - 0.5 * s)], 0.95 * s * ws, t0=0.2, t1=0.6, seed=seed, wob=0.08, pw=0.5, belly=0.25)
    elif style == "roof":
        stroke(ink, role, [(x - 1.15 * s, y + 0.35 * s), (x, y - 0.55 * s), (x + 1.15 * s, y + 0.35 * s)], 0.62 * s * ws, t0=0.2, t1=0.3, seed=seed, wob=0.1, smooth=False)
        stroke(ink, role, [(x - 0.7 * s, y + 0.9 * s), (x + 0.7 * s, y + 0.85 * s)], 0.45 * s * ws, t0=0.2, t1=0.2, seed=seed + 1, wob=0.1)
    elif style == "thread":
        stroke(ink, role, [(x - 1.2 * s, y + 0.15 * s), (x, y - 0.2 * s), (x + 1.2 * s, y + 0.2 * s)], 0.5 * s * ws, t0=0.25, t1=0.3, seed=seed, wob=0.12)
        stroke(ink, role, [(x + 1.1 * s, y + 0.2 * s), (x + 1.4 * s, y + 1.1 * s)], 0.38 * s * ws, t0=0.2, t1=0.6, seed=seed + 1, wob=0.15)
        blob(ink, role, x, y, 0.32 * s, 0.32 * s, seed=seed + 2, rough=0.1)
    elif style == "coil":
        from marks2 import spiral_pts
        stroke(ink, role, spiral_pts(x, y, 0.15 * s, 1.0 * s, 1.6, step=0.15), 0.4 * s * ws, t0=0.05, t1=0.2, seed=seed, wob=0.08, smooth=False, step=0.15)
    elif style == "plate":
        stroke(ink, role, [(x, y - 1.0 * s), (x + 0.05 * s, y + 1.0 * s)], 0.65 * s * ws, t0=0.05, t1=0.1, seed=seed, wob=0.1, belly=0.0, minw=0.8, cap=0.8)
        stroke(ink, role, [(x - 0.7 * s, y - 1.05 * s), (x + 0.75 * s, y - 1.15 * s)], 0.45 * s * ws, t0=0.2, t1=0.3, seed=seed + 1, wob=0.1)
    elif style == "curl":
        stroke(ink, role, [(x - 0.3 * s, y + 1.0 * s), (x - 0.9 * s, y + 0.1 * s), (x, y - 0.6 * s), (x + 0.9 * s, y - 0.2 * s), (x + 0.4 * s, y + 0.4 * s)], 0.5 * s * ws, t0=0.2, t1=0.5, seed=seed, wob=0.2, belly=0.3)
    elif style == "cell":
        stroke(ink, role, [(x - 0.75 * s, y), (x + 0.75 * s, y)], 1.4 * s * ws, t0=0.02, t1=0.02, seed=seed, wob=0.08, belly=0.0, minw=0.85, cap=0.85)
    elif style == "pine":
        stroke(ink, role, [(x, y + 1.1 * s), (x + 0.05 * s, y - 0.9 * s)], 0.3 * s * ws, t0=0.1, t1=0.3, seed=seed, wob=0.1)
        for k, (yy, l) in enumerate([(-0.6, 0.9), (0.0, 1.4), (0.6, 1.9)]):
            stroke(ink, role, [(x - l / 2 * s, y + (yy + 0.15) * s), (x, y + (yy - 0.15) * s), (x + l / 2 * s, y + (yy + 0.15) * s)], 0.42 * s * ws, t0=0.25, t1=0.3, seed=seed + 2 + k, wob=0.12)
    elif style == "bead":
        from brush import poly_d
        pts = [(x - 1.25 * s, y), (x - 0.15 * s, y - 0.7 * s), (x + 0.15 * s, y - 0.7 * s), (x + 1.25 * s, y), (x + 0.15 * s, y + 0.7 * s), (x - 0.15 * s, y + 0.7 * s)]
        ink.add(role, poly_d(pts))
    elif style == "blossom":
        for k in range(5):
            a = 0.3 + 2 * math.pi * k / 5
            blob(ink, role, x + math.cos(a) * 0.62 * s, y + math.sin(a) * 0.62 * s, 0.45 * s, 0.45 * s, seed=seed + k, rough=0.1)
    elif style == "ring":
        from marks4 import earc
        stroke(ink, role, earc(x, y, 1.25 * s, 0.5 * s, math.radians(-150), math.radians(160)), 0.42 * s * ws, t0=0.05, t1=0.5, seed=seed, wob=0.12, belly=0.08, minw=0.25)
        blob(ink, role, x, y, 0.3 * s, 0.18 * s, seed=seed + 1)
    elif style == "stone":
        blob(ink, role, x, y, 0.85 * s, 0.6 * s, rot=-0.3, seed=seed, rough=0.1)

def line(ink, role, pts, w, seed, t0=0.15, t1=0.35, **kw):
    stroke(ink, role, pts, w, t0=t0, t1=t1, seed=seed, wob=kw.pop("wob", 0.12), **kw)

def pack(style):
    """Eight SELAH feature icons in a 24-unit box, drawn in the hand of one artist."""
    icons = []
    # 1 Butcher tableaux: strictly lower-triangular A, then the weights b
    k = Ink()
    rows = [(4.5, 1), (9.0, 2), (13.5, 3)]
    for i, (y, m) in enumerate(rows):
        for j in range(m):
            element(k, style, 6.0 + 4.6 * j, y, 1.55, "ink", 10 * i + j)
    line(k, "ink", [(3.0, 16.9), (21.0, 16.7)], 1.0, 50, t0=0.1, t1=0.2)
    for j in range(4):
        element(k, style, 5.3 + 4.5 * j, 20.6, 1.55, "acc" if j == 3 else "ink", 60 + j)
    icons.append(("Butcher tableaux", k))
    # 2 SSP: a bowl that every stage stays inside
    k = Ink()
    line(k, "ink", [(3.0, 6.5), (5.5, 15.5), (12.0, 20.3), (18.5, 15.5), (21.0, 6.5)], 1.7, 70, t0=0.1, t1=0.3)
    line(k, "acc", [(7.0, 10.5), (12.0, 10.8), (17.0, 10.4)], 1.3, 71, t0=0.2, t1=0.3)
    for j, x in enumerate([8.4, 12.0, 15.6]):
        element(k, style, x, 15.0, 1.25, "ink", 72 + j)
    icons.append(("Strong stability (SSP)", k))
    # 3 Low storage: two registers, one stroke looping between them
    k = Ink()
    line(k, "ink", [(3.0, 6.0), (21.0, 5.8)], 1.8, 80, t0=0.08, t1=0.25)
    line(k, "ink", [(3.0, 18.2), (21.0, 18.0)], 1.8, 81, t0=0.08, t1=0.25)
    line(k, "acc", [(17.5, 8.3), (17.0, 12.0), (12.0, 12.6), (7.2, 12.4), (6.6, 15.8)], 1.2, 82, t0=0.1, t1=0.5)
    element(k, style, 17.6, 9.3, 1.2, "ink", 83)
    icons.append(("Low storage", k))
    # 4 Operator splitting: one path, two operators, rejoined
    k = Ink()
    line(k, "ink", [(2.5, 12.2), (7.5, 12.0)], 1.6, 90, t0=0.05, t1=0.05)
    line(k, "ink", [(7.5, 12.0), (12.0, 5.8), (16.5, 12.0)], 1.4, 91, t0=0.1, t1=0.1)
    line(k, "acc", [(7.5, 12.0), (12.0, 18.2), (16.5, 12.0)], 1.4, 92, t0=0.1, t1=0.1)
    line(k, "ink", [(16.5, 12.0), (21.5, 11.8)], 1.6, 93, t0=0.05, t1=0.5)
    element(k, style, 12.0, 12.0, 1.3, "ink", 94)
    icons.append(("Operator splitting", k))
    # 5 Interior-only kernels: a block, its halo left untouched
    k = Ink()
    sq = [((6.0, 6.0), (18.0, 6.1)), ((18.1, 6.0), (18.0, 18.0)), ((18.0, 18.1), (6.0, 18.0)), ((5.9, 18.0), (6.0, 6.0))]
    for j, (p, q) in enumerate(sq):
        line(k, "ink", [p, q], 1.4, 100 + j, t0=0.08, t1=0.12)
    for j in range(4):
        a = j * math.pi / 2
        for m in (-1, 1):
            cx, cy = 12 + 9.2 * math.cos(a), 12 + 9.2 * math.sin(a)
            tx, ty = -math.sin(a), math.cos(a)
            line(k, "ink55", [(cx + tx * (m * 2.5 - 1.6), cy + ty * (m * 2.5 - 1.6)), (cx + tx * (m * 2.5 + 1.6), cy + ty * (m * 2.5 + 1.6))], 0.8, 110 + 2 * j + (m > 0), t0=0.2, t1=0.2)
    for j, (x, y) in enumerate([(9.6, 9.6), (14.4, 9.6), (9.6, 14.4)]):
        element(k, style, x, y, 1.15, "ink", 120 + j)
    element(k, style, 14.4, 14.4, 1.15, "acc", 124)
    icons.append(("Interior-only kernels", k))
    # 6 CPU and GPU: one source, many lanes
    k = Ink()
    line(k, "ink", [(3.0, 5.5), (21.0, 5.3)], 1.9, 130, t0=0.06, t1=0.3)
    for j in range(6):
        x = 4.6 + 3.0 * j
        line(k, "acc" if j == 4 else "ink", [(x, 9.0), (x + 0.1, 19.5)], 1.15, 131 + j, t0=0.15, t1=0.45)
    icons.append(("CPU and GPU", k))
    # 7 Stage contract: begin_stage, end_stage, finish_step
    k = Ink()
    line(k, "ink", [(2.5, 15.0), (21.5, 14.8)], 1.3, 140, t0=0.05, t1=0.3)
    for j, x in enumerate([5.5, 11.5, 18.0]):
        line(k, "acc" if j == 2 else "ink", [(x, 10.5 if j < 2 else 8.5), (x + 0.1, 19.5)], 1.4 if j < 2 else 1.7, 141 + j, t0=0.12, t1=0.3)
    element(k, style, 8.5, 7.0, 1.3, "ink", 145)
    icons.append(("Stage contract", k))
    # 8 Order conditions: a rooted tree of Butcher's theory
    k = Ink()
    line(k, "ink", [(12.0, 20.5), (7.0, 13.0)], 1.25, 150, t0=0.05, t1=0.1)
    line(k, "ink", [(7.0, 13.0), (4.6, 6.0)], 1.1, 151, t0=0.05, t1=0.1)
    line(k, "ink", [(7.0, 13.0), (10.4, 6.0)], 1.1, 152, t0=0.05, t1=0.1)
    line(k, "acc", [(12.0, 20.5), (17.2, 12.6), (17.6, 6.4)], 1.25, 153, t0=0.05, t1=0.1)
    for j, (x, y) in enumerate([(12.0, 20.5), (7.0, 13.0), (4.6, 5.6), (10.4, 5.6), (17.6, 5.8)]):
        element(k, style, x, y, 1.35, "ink", 160 + j)
    icons.append(("Order conditions", k))
    return icons
