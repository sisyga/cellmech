#!/usr/bin/python  -u

from cell2opt import *
import cProfile
import matplotlib.pyplot as plt

npr.seed(seed=0)

#######################################################

def generate_config_from_default(R):
    N = len(R)
    c = Configuration(N, dims=2)
    for ni in range(N):
        c.nodesX[ni] = R[ni]

    return c

if __name__ == '__main__':

    bend = 10.0
    twist = 1.0
    Lmax = 3
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

    N = 3

    # config = generate_initial_config(Lmax)

    xpos = np.sqrt(2 * 1.3 * 1.3) / 2.

    R = np.array([[-xpos, xpos, 0], [0, 0., 0], [xpos, xpos, 0], [0, 2 * xpos, 0]])
    config = generate_config_from_default(R)

    config.updateDists(config.nodesX)

    config.addlink(0, 1)
    config.addlink(1, 2)
    config.addlink(2, 3)
    config.addlink(3, 0)

    config.addlink(0, 2)

    configs, links, nodeforces, linkforces, ts = config.minitimeevo(4., record=True, now=False)
    """
    t = 0.
    config.makesnap(0)

    thisdt = config.mechEquilibrium(now=False)
    print "equil: dt = ", thisdt
    t += thisdt
    print "t = ", t
    config.makesnap(t)


    myd0 = scipy.linalg.norm(config.nodesX[0] - config.nodesX[2])



    config.addlink(0, 2, d0 = 2.)
    config.updateDists(config.nodesX)

    thisdt = config.mechEquilibrium(now=False)
    print "equil: dt = ", thisdt
    t += thisdt
    print "t = ", t
    config.makesnap(t)

    Qtie = np.linspace(1, len(config.Qtrack), len(config.Qtrack))
    plt.plot(Qtie, config.Qtrack, color="blue")
    Ferrtie = np.linspace(1, len(config.Ferrtrack), len(config.Ferrtrack))
    plt.plot(Ferrtie, config.Ferrtrack, color="red")
    # plt.show()

    config.nodesnap = np.array(config.nodesnap)
    config.fnodesnap = np.array(config.fnodesnap)
    """
    animateconfigs(config.nodesnap, config.linksnap, config.fnodesnap, config.flinksnap, config.snaptimes)
    mlab.show()
