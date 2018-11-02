#!/usr/bin/python  -u

from cell import *
import cProfile
import matplotlib.pyplot as plt


#######################################################

def generate_config_from_default(R):
    N = len(R)
    c = Configuration(num=N, d2=True)
    for ni in range(N):
        n = node(c, R[ni])
        n.twist = twist
        n.bend = bend

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

    N = 4

    # config = generate_initial_config(Lmax)

    xpos = np.sqrt(2 * 1.3 * 1.3) / 2.

    R = np.array([[-xpos, xpos, 0], [0, 0., 0], [xpos, xpos, 0], [0, 2 * xpos, 0]])
    config = generate_config_from_default(R)

    for i in range(N):
        config.nodes[i].addLinkTo(config.nodes[(i + 1)%N], d0 = 1.)

    configs, ts = [copy.deepcopy(config)], [0.]
    t = 0.
    dt = config.mechEquilibrium()
    t += dt
    configs.append(copy.deepcopy(config))
    ts.append(t)

    config.nodes[0].addLinkTo(config.nodes[2], d0 = 2.)

    dt = config.mechEquilibrium()
    t += dt
    configs.append(copy.deepcopy(config))
    ts.append(t)

    # Qtie = np.linspace(1, len(config.Qtrack), len(config.Qtrack))
    # plt.plot(Qtie, config.Qtrack)
    # plt.show()


    animateconfigs(configs, ts=ts)
    mlab.show()