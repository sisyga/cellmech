#!/usr/bin/python  -u

from cell2 import *
from scipy.spatial import Delaunay
import cProfile

#######################################################


def VoronoiNeighbors(positions, is3D=False):
    p = [n[:2] for n in positions]
    # tri: list of interconnectedparticles: [ (a, b, c), (b, c, d), ... ]
    tri = Delaunay(p, qhull_options='QJ')
    # neighbors contain pairs of adjacent particles: [ (a,b), (c,d), ... ]
    neighbors = [list(itertools.combinations(v, 2)) for v in tri.simplices]
    # neighbors = set([tuple(itertools.combinations(v, 2)) for v in tri.simplices])
    # neighbors = tri.neighbors
    n = []
    for (i, j) in itertools.chain.from_iterable(neighbors):
        if i < j:
            n.append((i, j))
        else:
            n.append((j, i))
    neighbors = set(n)
    return neighbors


def generatePoint(L):
    X0 = (npr.rand() - .5) * L
    Y0 = (npr.rand() - .5) * L
    Z0 = 0.
    return np.array([X0, Y0, Z0])


def rand3d(L):
    return (np.random.random(3) - 0.5) * L


def generate_initial_config(L=10, N=None):
    if N is None:
        N = int(L ** 2)
    c = Configuration(N)

    for ni in range(N):
        while True:
            R1 = rand3d(L)
            OK = True
            for nj in range(ni):
                d = np.linalg.norm(c.nodes[nj][0] - R1)
                if d < d0min:
                    OK = False
                    break
            if OK:
                break
        c.nodes[ni][0] = R1

    return c


if __name__ == '__main__':

    bend = 10.0
    twist = 1.0
    Lmax = 10
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

    config = generate_initial_config(Lmax)

    for i, j in VoronoiNeighbors([config.nodes[:, 0, :]]):
        if np.linalg.norm(i - j) <= d0max:
            n1 = np.where(config.nodes == i)
            n2 = np.where(config.nodes == j)
            config.addlink(n1, n2)

    cProfile.run('config.timeevo(2, record=True)', sort='cumtime')
    # configs, ts = c.timeevo(200, record=True)
    # animateconfigs(configs, ts=ts)
    # mlab.show()
