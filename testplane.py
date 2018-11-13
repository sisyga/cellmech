#!/usr/bin/python  -u

from cell import *
import cProfile
import animate
import matplotlib.pyplot as plt
npr.seed(seed=0)

#######################################################


def generatePoint(L):
    X0 = (npr.rand() - .5) * L
    Y0 = (npr.rand() - .5) * L
    Z0 = 0.
    return np.array([X0, Y0, Z0])


def rand3d(L):
    return (np.random.random(3) - 0.5) * L


def generate_initial_config(L, N, dt, nmax, qmin, d0_0, p_add, p_del, chkx, d0max, dims):
    if N is None:
        N = int(L ** 2)

    c = Configuration(N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max,
                      dims=dims)

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


def generate_default_initial(L=10, N=None):
    if N is None:
        N = int(L ** 2)
    R = []
    for ni in range(N):
        while True:
            R1 = generatePoint(L)
            OK = True
            for r in R:
                d = np.linalg.norm(R1 - r)
                if d < d0min:
                    OK = False
                    break
            if OK:
                break
        R.append(R1)
    R = np.array(R)
    np.save("Rinit", R)

    return R


def generate_config_from_default(R):
    N = len(R)
    c = Configuration(N, dims=2)
    for ni in range(N):
        c.nodesX[ni] = R[ni]

    return c


def generate_cube(L, N, dt, nmax, qmin, d0_0, p_add, p_del, chkx, d0max, dims, stretch=1.):
    if N is not None:
        print "N decission overriden because of possible conflict with Lmax"
    N = int(L ** 3)

    c = Configuration(N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max,
                      dims=dims)

    x0 = - stretch * L / 2.
    for ni in range(L):
        for nj in range(L):
            for nk in range(L):
                c.nodesX[L * L * ni + L * nj + nk] = np.array([x0 + stretch * ni, x0 + stretch * nj, x0 + stretch * nk])
    return c


if __name__ == '__main__':

    Lmax = 20
    N = None

    bend = 10.0
    twist = 1.0
    dt = 0.01
    nmax = 3000
    qmin = 0.001
    dims = 2

    d0min = 0.8  # min distance between cells
    d0max = 2.  # max distance connected by links
    d0_0 = 1.  # equilibrium distance of links
    p_add = 1.0  # rate to add links
    p_del = 0.1  # base rate to delete links
    chkx = True  # check if links overlap?

    config = generate_initial_config(L=Lmax, N=N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del,
                                     chkx=chkx, d0max=d0max, dims=dims)

    config.updateDists(config.nodesX)

    for i, j in voronoi_neighbors.VoronoiNeighbors(config.nodesX, d0max=config.d0max, vodims=2):
        if np.linalg.norm(config.nodesX[i] - config.nodesX[j]) <= d0max:
            config.addlink(i, j)

    cProfile.run('config.timeevo(2, record=True)', sort='cumtime')
    # configs, links, nodeforces, linkforces, ts = config.timeevo(5., record=True)
    # animateconfigs(configs, links, nodeforces, linkforces, ts)
    # mlab.show()

