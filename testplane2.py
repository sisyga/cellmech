#!/usr/bin/python  -u

from cell2 import *
import cProfile

npr.seed(seed=0)

#######################################################


def generatePoint(L):
    X0 = (npr.rand() - .5) * L
    Y0 = (npr.rand() - .5) * L
    Z0 = 0.
    return np.array([X0, Y0, Z0])


def rand3d(L):
    return (np.random.random(3) - 0.5) * L


def generate_initial_config(L=10, N=None):
    if N is None:
        # N = int(L ** 2)
        N = 70
    c = Configuration(N, dims=2)

    for ni in range(N):
        while True:
            R1 = generatePoint(L)
            OK = True
            for nj in range(ni):
                d = np.linalg.norm(c.nodesX[nj] - R1)
                if d < d0min:
                    OK = False
                    break
            if OK:
                break
        c.nodesX[ni] = R1

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
    config.updateDists(config.nodesX)

    for i, j in voronoi_neighbors2.VoronoiNeighbors(config.nodesX, d0max=config.d0max, vodims=2):
        if np.linalg.norm(config.nodesX[i] - config.nodesX[j]) <= d0max:
            config.addlink(i, j)

    # cProfile.run('config.timeevo(2, record=True)', sort='cumtime')
    configs, links, nodeforces, linkforces, ts = config.timeevo(2, record=True)
    animateconfigs(configs, links, nodeforces, linkforces, ts)
    mlab.show()