#!/usr/bin/python  -u

from cell import *
from animate import *
import cProfile
import matplotlib.pyplot as plt

npr.seed(seed=0)

#######################################################


def generatePoint(L):
    X0 = (npr.rand() - .5) * L
    Y0 = (npr.rand() - .5) * L
    Z0 = 0.
    return np.array([X0, Y0, Z0])


def generate_config_from_default(R, L, N, dt, nmax, qmin, d0_0, p_add, p_del, chkx, d0max, dims):
    N = len(R)
    c = Configuration(N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max,
                      dims=dims)
    for ni in range(N):
        c.nodesX[ni] = R[ni]

    return c


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

def generate_initial_bilayer(L, N, dt, nmax, qmin, d0_0, p_add, p_del, chkx, d0max, dims):
    if N is None:
        N = int(2 * (L ** 2))

    c = Configuration(N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max,
                      dims=dims)

    for ni in range(N):
        while True:
            R1 = generatePoint(L)
            if ni >= N / 2:
                R1[2] = d0_0
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

    Lmax = 5
    N = None

    bend = 10.0
    twist = 1.0
    dt = 0.01
    nmax = 3000
    qmin = 0.001
    dims = 3

    d0min = 0.8  # min distance between cells
    d0max = 2.  # max distance connected by links
    d0_0 = 1.  # equilibrium distance of links
    p_add = 1.0  # rate to add links
    p_del = 0.1  # base rate to delete links
    chkx = True  # check if links overlap?

    # config = generate_initial_config(Lmax)

    xpos = np.sqrt(2 * 1.3 * 1.3) / 2.
    # R = np.array([[0, 0, 0], [1.3, 0, 0]])
    # R = np.array([[-xpos, xpos, 0], [0, 0., 0], [xpos, xpos, 0], [0, 2 * xpos, 0]])
    # R = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    # R = np.array([[0, 0, 0], [0, 1, 0], [2, 2, 0], [0.5, 1, 0]])
    # R = np.array([[0.09762701, 0.43037873, 0.], [0.92732552, -0.23311696,  0.], [-0.85792788, -0.8257414,  0.], [-0.95956321,  0.66523969,  0.]])

    R = np.array([[0, 0, 0], [1, 0, 0], [1, 1.5, 0], [0.8, 1.2, 0], [1.6, 0.8, 1.5]])

    config = generate_config_from_default(R, L=Lmax, N=N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max, dims=dims)

    # config = generate_initial_config(L=Lmax, N=N, dt=dt, nmax=nmax, qmin=qmin, d0_0=d0_0, p_add=p_add, p_del=p_del, chkx=chkx, d0max=d0max, dims=dims)
    """
    for i, j in VoronoiNeighbors(config.nodesX, d0max=config.d0max, vodims=2):
        if np.linalg.norm(config.nodesX[i] - config.nodesX[j]) <= d0max:
            config.addlink(i, j, d0=1.)

    """

    config.updateDists(config.nodesX)

    config.addlink(0, 1, d0 = 1.)
    config.addlink(1, 2, d0 = 1.)
    config.addlink(2, 3, d0 = 1.)
    config.addlink(3, 0, d0 = 1.)
    config.addlink(1, 3, d0 = 1.)

    config.addlink(0, 4, d0 = 1.)
    config.addlink(1, 4, d0 = 1.)
    config.addlink(2, 4, d0 = 1.)
    config.addlink(3, 4, d0 = 1.)


    # cProfile.run('config.oneequil()', sort='cumtime')

    configs, links, nodeforces, linkforces, ts = config.oneequil()
    animateconfigs(configs, links, nodeforces, linkforces, ts)
    mlab.show()

