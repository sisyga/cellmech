#!/usr/bin/python  -u

import numpy as np
import mayavi.mlab as mlab
from cell import *
from voronoi_neighbors import *

#######################################################


def generatePoint(L):
    X0 = (np.random.rand() - .5) * L
    Y0 = (np.random.rand() - .5) * L
    Z0 = 0.
    return np.array([X0, Y0, Z0])

def rand3d(L):
    return (np.random.random(3) - 0.5) * L

def generate_initial_config(L=10, d2=True, N=None):
    if d2:
        if N is None:
            N = int(L**2)
        c = Configuration(N, d2=d2)

        for i in range(N):
            while True:
                R1 = generatePoint(L)
                OK = True
                for n in c.nodes:
                    d = norm(n.getR() - R1, "len")
                    if d < d0min:
                        OK = False
                        break
                if OK: break
            n = node(c, R1)
            n.twist = twist
            n.bend = bend

        return c

    else:
        if N is None:
            N = int(L ** 2)
        c = Configuration(N, d2=d2)

        for i in range(N):
            while True:
                R1 = rand3d(L)
                OK = True
                for n in c.nodes:
                    d = norm(n.getR() - R1, "len")
                    if d < d0min:
                        OK = False
                        break
                if OK: break
            n = node(c, R1)
            n.twist = twist
            n.bend = bend

        return c



if __name__ == '__main__':

    bend = 10.0
    twist = 1.0
    L = 10.
    dt = 0.01
    nmax = 3000
    qmin = 0.001
    d2 = True

    d0min = 0.8  # min distance between cells
    d0max = 2.  # max distance connected by links
    d0_0 = 1.  # equilibrium distance of links
    p_add = 1.0  # rate to add links
    p_del = 0.1  # base rate to delete links
    anis = 1.0  # anisotropy of building links (we don't need this so far)
    chkx = True  # check if links overlap?

    N = int(L ** 2)  # 1 cell / unit area
    c = generate_initial_config(L, N=N, d2=d2)

    for i, j in VoronoiNeighbors([n.getR() for n in c.nodes], d0max, is3D=c.is3d):
        if norm(c.x[c.nodes[i].r] - c.x[c.nodes[j].r], 'mag') <= d0max:
            c.nodes[i].addLinkTo(c.nodes[j])

    configs, ts = c.timeevo(25, record=True)
    animateconfigs(configs, ts=ts)
    mlab.show()
