import math
from brush import Ink, stroke, blob, drop
from marks4 import earc
from marks5 import ring
def fav_mark():
    ink = Ink()
    drop(ink, "acc", 32.0, 15.0, 8.0)
    ring(ink, "ink", 32.0, 44.5, 9.0, 3.4, 4.2, 900, gap0=-130, span=310, t1=0.35)
    ring(ink, "ink55", 32.0, 44.5, 22.0, 8.2, 3.6, 901, gap0=-20, span=300, t1=0.4)
    blob(ink, "ink", 32.0, 44.5, 2.2, 1.2, seed=902)
    return ink
