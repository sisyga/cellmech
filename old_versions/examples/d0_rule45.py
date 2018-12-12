import cell
import math
import random
import numpy as np


def check_d0(d0):
    if (d0 < 0.1):
        d0 = 0.1
    if (d0 > 2.0):
        d0 = 2.0
    return d0


def update_d0(g):
    adjust_rate = 0.01
    for l in cell.links:
        F = -np.dot(l.F1, l.e)  # F<0 for tensile link
        Ftarget = float(g['defaults']['Ftgt']) * (2.0 - l.d0)
        l.d0 = (l.d0 +
                adjust_rate * (F - Ftarget) * g['dt'] +
                0.05 * random.uniform(
                    -math.sqrt(g['dt']), math.sqrt(g['dt']))
                )
        l.d0 = check_d0(l.d0)
    for n in cell.nodes:
        if n.kr < 99:  # not a fixed boundary
            dr0 = n.kr * (cell.x[n.r] - n.r0) * g['dt']
            n.r0 = n.r0 + float(g['defaults']['adhesionAdjust']) * dr0
