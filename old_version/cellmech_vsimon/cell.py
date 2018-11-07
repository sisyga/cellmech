#!/usr/bin/python  -u

PYTHONUNBUFFERED = 1

import argparse
import collections
import copy
import gzip
import itertools
import os
import pickle
import sys
from math import exp, log, sqrt

import numpy as np
from mayavi import mlab
from numpy.random import random

import voronoi_neighbors

null = np.array([0.0, 0.0, 0.0])
ex = np.array([1.0, 0.0, 0.0])
ey = np.array([0.0, 1.0, 0.0])
ez = np.array([0.0, 0.0, 1.0])


def norm(v, what=None):
    n = np.sqrt(np.dot(v, v))
    if n > 1e-6:
        u = v / n
    else:
        u, n = null, 0.
    if what == 'vec':
        return u
    if what == 'mag' or what == 'len':
        return n
    return [u, n]


def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


def intersect(A, B, C, D):
    if norm(A - C, "mag") < 0.01: return False
    if norm(A - D, "mag") < 0.01: return False
    if norm(B - C, "mag") < 0.01: return False
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def R(dphi, v):
    axis, theta = norm(dphi)
    if theta < 1e-6:
        return v
    a = np.cos(theta / 2)
    b, c, d = axis * np.sin(theta / 2)
    return np.dot(np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                            [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]]), v)


def update_progress(progress):
    """
    Simple progress bar update.
    :param progress: float. Fraction of the work done, to update bar.
    :return:
    """
    barLength = 20  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\rProgress: [{0}] {1} % {2}".format("#" * block + "-" * (barLength - block), round(progress * 100, 1),
                                               status)
    sys.stdout.write(text)
    sys.stdout.flush()


class Configuration():
    def __init__(self, num=None, dumpdir=None, dt=0.01, nmax=3000, qmin=0.001, d0_0=1, force_limit=15., p_add=1.,
                 p_del=0.2, chkx=True, anis=0.0, d0max=2., d2=False, unrestricted=True, distlimit=3):
        self.nodes = []
        self.links = []
        self.nNodes = 0  # bookkeeping for state vector position (I don't get that)
        self.dumpdir = dumpdir  # string for directory to save configs
        # parameters for mechanical equilibration
        self.dt = dt
        self.nmax = nmax
        self.qmin = qmin
        # 2d or 3d?
        self.d2 = d2
        self.is3d = 1 - self.d2
        # parameters to add/remove links
        self.d0_0 = d0_0
        self.force_limit = force_limit
        self.p_add = p_add
        self.p_del = p_del
        self.chkx = chkx
        self.anis = anis
        self.d0max = d0max
        # parameters to restrict mech. equilibration to surrounding of changed link
        self.unrestricted = unrestricted
        self.distlimit = distlimit
        # variables to store cell number and cell positions and angles
        if num is not None:
            self.N = num
            self.x = np.zeros(6 * num)
        else:
            self.N = 0
            self.x = None

    def updateLinkForces(self, l, x):
        l.e, l.d = norm(x[l.n2.r] - x[l.n1.r])  # actual link direction and length
        l.n = norm(l.norm1 + l.norm2, "vec")  # Eq. 7
        norm2 = np.cross(l.e, l.n)  # Eq. 9
        if self.d2:
            l.M1 = l.n1.bend * np.cross(l.e, R(x[l.n1.phi], l.t1))
            l.M2 = l.n2.bend * np.cross(-l.e, R(x[l.n2.phi], l.t2))
        else:
            l.M1 = l.n1.bend * np.cross(l.e, R(x[l.n1.phi], l.t1)) + \
                   l.n1.twist * np.cross(l.n, R(x[l.n1.phi], l.norm1))  # Eqs 3, 18, 19
            l.M2 = l.n2.bend * np.cross(-l.e, R(x[l.n2.phi], l.t2)) + \
                   l.n2.twist * np.cross(l.n, R(x[l.n2.phi], l.norm2))  # Eqs 5, 18, 19
        M = l.M1 + l.M2
        l.F2 = np.dot(M, l.n) / l.d * norm2 - np.dot(M, norm2) / l.d * l.n + l.k * (
                l.d - l.d0) * l.e  # Eqs (11), (12) and (13) for F2=-F1!
        l.F1 = -l.F2

    def F(self, x):
        ret = np.zeros(len(x))  # change in state vector, 6 dof for each unit
        for l in self.links:
            self.updateLinkForces(l, x)
        for n in self.nodes:
            n.F = -n.kr * (x[n.r] - n.r0) + n.F0
            n.M = -n.kphi * (x[n.phi] - n.phi0)
        for l in self.links:
            l.n1.F = l.n1.F - l.F1
            l.n2.F = l.n2.F - l.F2
            l.n1.M = l.n1.M - l.M1
            l.n2.M = l.n2.M - l.M2
        for n in self.nodes:
            ret[n.r] = n.F
            ret[n.phi] = n.M
        return ret

    def mechEquilibrium(self):
        x = self.x.copy()
        h = self.dt
        for i in range(self.nmax):
            k1 = self.F(x)
            Q = np.dot(k1, k1) / len(self.nodes)
            # print i, Q, max([x[n.r][0] for n in self.nodes])
            if i % 100 == 0 and self.dumpdir is not None:
                dumpData('s/%05d.pickle' % (self.dumpdir, i))
            if Q < self.qmin:
                break
            k1 *= h
            k2 = h * self.F(x + k1 / 2)
            k3 = h * self.F(x + k2 / 2)
            k4 = h * self.F(x + k3)
            x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.
        self.x = x
        return (i + 1) * h

    def delAttrDist(self):
        for n in self.nodes:
            if hasattr(n, 'dist'):
                delattr(n, "dist")

    def checkLinkX(self):
        Xs = []
        delete_list = []
        for l1, l2 in itertools.combinations(self.links, 2):
            if intersect(l1.n1.getR(), l1.n2.getR(), l2.n1.getR(), l2.n2.getR()):
                Xs.append([l1, l2])
        # print "Xs:", Xs
        while len(Xs) > 0:
            # Counter: count occurence
            # from_iterable: unpack list of lists
            counts = collections.Counter(itertools.chain.from_iterable(Xs))
            # get item most in conflict (highest occurence in list)
            l = max(counts, key=counts.get)
            delete_list.append(l)
            newXs = [x for x in Xs if l not in x]
            Xs = newXs
        for l in delete_list:
            l.remove()

    def delLinkList(self):
        to_del = []
        for l in self.links:
            if l.d < l.d0:
                continue  # compressed links are stable
            f = norm(l.F1, "mag")  # magnitude of force transmitted
            if hasattr(l, "limit"):
                p = exp(f / l.limit)
            else:
                p = exp(f / self.force_limit)
            if self.anis > 1e-5:  # horizontal links are more stable
                p = p * (1 - self.anis + self.anis * pow((l.n1.getR() - l.n2.getR())[1] / l.d, 2))  # what?
            to_del.append((l, p))
        return to_del

    def tryLink(self, n, m):
        if n.findLinkTo(m) != None:
            return -1, null
        for l in self.links:
            if intersect(self.x[n.r], self.x[m.r], self.x[l.n1.r], self.x[l.n2.r]):
                return -1, null  # false
        e, d = norm(self.x[n.r] - self.x[m.r])
        if d > self.d0max:
            return -1, null  # false
        return d, e  # true: d>0

    def addLinkList(self):
        to_add = []
        for i, j in voronoi_neighbors.VoronoiNeighbors(
                [n.getR() for n in self.nodes], self.d0max, is3D=self.is3d):
            d, e = self.tryLink(self.nodes[i], self.nodes[j])
            if d > 1e-5:
                p = (1 - (d / 2.0)) * (1 - self.anis + self.anis * e[0] * e[0])
                to_add.append(((self.nodes[i], self.nodes[j]), p))
        return to_add

    def pickEvent(self, to_del, to_add):
        s1 = 0.
        for (l, p) in to_del:
            s1 += p * self.p_del
        s2 = 0.
        for (q, p) in to_add:
            s2 += p * self.p_add

        S = s1 + s2
        if S < 1e-7:
            print "nothing to do!"
            return 1.
        dt = -log(random()) / S
        if dt > 1:
            print 'Must adjust d0 variables before the next event!'
            return 1.

        r = S * random()
        if r < s1:  # we will remove a link
            for (l, p) in to_del:
                r = r - p * self.p_del
                if r < 0:
                    l.n1.dist = 0
                    l.n2.dist = 0
                    l.remove()
                    return dt
        r = r - s1
        if r < s2:  # we will add a link
            for ((n1, n2), p) in to_add:
                r = r - p * self.p_add
                if r < 0:
                    n1.dist = 0
                    n2.dist = 0
                    n1.addLinkTo(n2)
                    return dt

    def default_update_d0(self, dt):
        for l in self.links:
            l.d0 += 0.2 * (self.d0_0 - l.d0) * dt + 0.05 * (
                        2 * sqrt(dt) * random() - sqrt(dt))  # magic number 0.2 and 0.05??

    def modlink(self):
        self.delAttrDist()
        if self.chkx:
            self.checkLinkX()

        to_del = self.delLinkList()
        to_add = self.addLinkList()
        dt = self.pickEvent(to_del, to_add)
        self.default_update_d0(dt)
        return dt

    def timeevo(self, tmax, record=False):
        # configs, ts = [copy.deepcopy(self)], [0.]
        t = 0.
        while t < tmax:
            dt = self.mechEquilibrium()
            t += dt
            if record:
                configs.append(copy.deepcopy(self))
                ts.append(t)
            dt = self.modlink()
            update_progress(t / tmax)
            t += dt
            if record:
                configs.append(copy.deepcopy(self))
                ts.append(t)
            update_progress(t / tmax)
        if record:
            return configs, ts


def showconfig(config, figure=None, figureindex=0, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), figsize=(1000, 1000),
               cmap='viridis', vmaxlinks=5, vmaxcells=5,
               cbar=False):
    if figure is None:
        fig = mlab.figure(figureindex, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    else:
        fig = figure
    x, y, z, phi1, phi2, phi3 = config.x.reshape(6, config.N, order='F')
    posl, rl, dl, fl, fc = [], [], [], [], []
    for l in config.links:
        posl.append(config.x[l.n1.r])
        rl.append(l.e)
        dl.append(l.d)
        fl.append(norm(l.F1, 'mag'))

    for c in config.nodes:
        fc.append(norm(c.F, 'mag'))

    dl = np.asarray(dl)
    fl = np.asarray(fl)
    xl, yl, zl = np.asarray(posl).T
    rxl, ryl, rzl = np.asarray(rl).T
    rxl *= dl
    ryl *= dl
    rzl *= dl

    cells = mlab.points3d(x, y, z, fc, scale_factor=1, opacity=0.5, resolution=16, scale_mode='none', vmin=0.,
                          colormap=cmap, vmax=vmaxcells)
    links = mlab.quiver3d(xl, yl, zl, rxl, ryl, rzl, scalars=fl, mode='2ddash', line_width=4., scale_mode='vector',
                          scale_factor=1, colormap=cmap, vmin=0., vmax=vmaxlinks)
    links.glyph.color_mode = "color_by_scalar"
    if cbar:
        mlab.scalarbar(links, nb_labels=2, title='Force on link')
    # mlab.draw()
    return cells, links


@mlab.animate(delay=500)
def animateconfigs(configs, ts=None, figureindex=0, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), figsize=(1000, 1000),
                   cmap='viridis', cbar=False):
    fig = mlab.figure(figureindex, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    c = configs[0]
    vmaxcells = max([max([norm(node.F, 'mag') for node in config.nodes]) for config in configs])
    vmaxlinks = max([max([norm(link.F1, 'mag') for link in config.links]) for config in configs])
    cells, links = showconfig(c, figure=fig, vmaxcells=vmaxcells, vmaxlinks=vmaxlinks)
    text = mlab.title('0.0', height=.9)  # (cells.mlab_source.x.max(), cells.mlab_source.y.max(), cells.mlab_source.z.max(), '0.0')
    if ts is None:
        ts = range(len(configs))
    while True:
        for (c, t) in zip(configs, ts):
            x, y, z, phi1, phi2, phi3 = c.x.reshape(6, c.N, order='F')
            posl, rl, dl, fl, fc = [], [], [], [], []
            for l in c.links:
                posl.append(c.x[l.n1.r])
                rl.append(l.e)
                dl.append(l.d)
                fl.append(norm(l.F1, 'mag'))

            for n in c.nodes:
                fc.append(norm(n.F, 'mag'))

            dl = np.asarray(dl)
            fl = np.asarray(fl)
            xl, yl, zl = np.asarray(posl).T
            rxl, ryl, rzl = np.asarray(rl).T
            rxl *= dl
            ryl *= dl
            rzl *= dl
            cells.mlab_source.set(x=x, y=y, z=z, scalars=fc)
            links.mlab_source.reset(x=xl, y=yl, z=zl, u=rxl, v=ryl, w=rzl, scalars=fl)
            text.set(text='{}'.format(round(t, 2)))
            print 'Updating... '
            yield

class node(object):
    def __init__(self, config, r, r0=None, state=-1):

        self.config = config  # reference to a specific configuration
        # state vector indices
        if config.nNodes == len(config.x):
            config.N += 1
            config.x.append(np.zeros(6))

        self.r = [config.nNodes, config.nNodes + 1, config.nNodes + 2]
        config.nNodes += 3
        config.x[self.r] = r
        if r0 is None:
            r0 = r
        self.r0 = r0
        # state vector indices
        self.phi = [config.nNodes, config.nNodes + 1, config.nNodes + 2]
        config.nNodes += 3
        config.x[self.phi] = null
        self.phi0 = null

        self.links = []
        self.F = null  # net force
        self.F0 = null  # external force
        self.M = null  # net torque

        self.bend = 10.  # bending rigidity
        self.twist = 1.  # torsion spring constant
        self.kr = 0.  # spring to anchor point
        self.kphi = 0.  # torsion spring to anchor point

        self.state = state  # discrete state variable, like 0, 1, 2, ..
        self.normal = None  # normal vector, if surface node
        config.nodes.append(self)

    def findLinkTo(self, n):
        for l in self.links:
            if ((l.n1 == self) & (l.n2 == n)) | ((l.n1 == n) & (l.n2 == self)):
                return l
        return None

    def addLinkTo(self, n, t1=None, t2=None, d0=None, k=None):
        if self.findLinkTo(n) is None:
            l = link(self.config, self, n, t1, t2, d0, k)
            self.links.append(l)
            n.links.append(l)
            return l

    def rmLinkTo(self, n):
        l = self.findLinkTo(n)
        if l is not None:
            self.links.remove(l)
            n.links.remove(l)
            self.config.links.remove(l)

    def getR(self):
        return self.config.x[self.r]


class link(object):
    def __init__(self, config, n1, n2, t1=None, t2=None, d0=None, k=None, n=None, norm1=None, norm2=None):

        self.config = config  # reference to a specific configuration
        self.n1, self.n2 = n1, n2  # nodes at endpoints
        # actual direction e and length d
        self.e, self.d = norm(config.x[n2.r] - config.x[n1.r])  # vector and its length

        if k is None:
            k = 15.
        self.k = k  # spring parameter

        if d0 is None:
            d0 = self.d
        self.d0 = d0  # equilibrium distance
        # preferred directions
        if t1 is None:
            self.t1 = R(-config.x[self.n1.phi], self.e)
        else:
            self.t1 = norm(t1, "vec")
        if t2 is None:
            self.t2 = R(-config.x[self.n2.phi], -self.e)
        else:
            self.t2 = norm(t2, "vec")
        if n is None:
            n, q = norm(np.cross(self.e, ez))  # n is perpendicular to e
            # n is perpendicular to z (l is in the x-y plane)
            if q < 1e-5:
                n = norm(np.cross(self.e, ex), "vec")  # e || ez   =>	n is perpendicular to x
        self.n = n  # perpendicular vector to link
        if norm1 is None:
            norm1 = n
        self.norm1 = norm1
        if norm2 is None:
            norm2 = n
        self.norm2 = norm2
        self.M1, self.M2, self.F1, self.F2 = null, null, null, null
        config.links.append(self)

    def remove(self):
        self.n1.links.remove(self)
        self.n2.links.remove(self)
        self.config.links.remove(self)


def checkDir(s):
    if os.path.exists(s):
        if not os.path.isdir(s):
            raise argparse.ArgumentTypeError("%s is not a dir." % s)
    else:
        os.makedirs(s)
    return s


def checkFile(s):
    if os.path.exists(s):
        if s.strip() == '/dev/stdin':
            return s
        if not os.path.isfile(s):
            raise argparse.ArgumentTypeError("%s is not a file." % s)
        return s
    else:
        argparse.ArgumentTypeError("%s does not exists." % s)


def dumpData(configuration, fname):
    sys.setrecursionlimit(50000)
    with gzip.GzipFile(fname, 'wb') as out:
        pickle.dump(configuration, out)


def loadData(fname):
    with gzip.GzipFile(fname, 'rb') as inf:
        config = pickle.load(inf)

    print >> sys.stderr, "N=", config.N
    print >> sys.stderr, len(config.x), "dof."
    print >> sys.stderr, len(config.nodes), "nodes."
    print >> sys.stderr, len(config.links), "links."
    return config
