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
import warnings

import numpy as np
from mayavi import mlab
from numpy.random import random

import voronoi_neighbors

warnings.filterwarnings("error")

null = np.array([0.0, 0.0, 0.0])
ex = np.array([1.0, 0.0, 0.0])
ey = np.array([0.0, 1.0, 0.0])
ez = np.array([0.0, 0.0, 1.0])

def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


def getNormvec(v):
    # returns normalized v
    d = np.linalg.norm(v)
    try:
        v = v/d
    except RuntimeWarning:
        pass
    return v

def getNormtoo(v):
    # returns norm of v and normalized v
    d = np.linalg.norm(v)
    try:
        v = v/d
    except RuntimeWarning:
        d = 0
    return v, d

def intersect(A, B, C, D):
    if np.linalg.norm(A - C) < 0.01:
        return False
    if np.linalg.norm(A - D) < 0.01:
        return False
    if np.linalg.norm(B - C) < 0.01:
        return False
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def R(dphi, v):
    axis, theta = getNormtoo(dphi)
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
        fl.append(np.linalg.norm(l.F1))

    for c in config.nodes:
        fc.append(np.linalg.norm(c.F))

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


@mlab.animate(delay=100)
def animateconfigs(configs, ts=None, figureindex=0, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), figsize=(1000, 1000),
                   cmap='viridis', cbar=False):
    fig = mlab.figure(figureindex, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    c = configs[0]
    vmaxcells = max([max([np.linalg.norm(node.F) for node in config.nodes]) for config in configs])
    vmaxlinks = max([max([np.linalg.norm(link.F1) for link in config.links]) for config in configs])
    cells, links = showconfig(c, figure=fig, vmaxcells=vmaxcells, vmaxlinks=vmaxlinks)
    # (cells.mlab_source.x.max(), cells.mlab_source.y.max(), cells.mlab_source.z.max(), '0.0')
    text = mlab.title('0.0', height=.9)
    if ts is None:
        ts = range(len(configs))
    while True:
        for (c, t) in zip(configs[1:], ts[1:]):
            x, y, z, phi1, phi2, phi3 = c.x.reshape(6, c.N, order='F')
            posl, rl, dl, fl, fc = [], [], [], [], []
            for l in c.links:
                posl.append(c.x[l.n1.r])
                rl.append(l.e)
                dl.append(l.d)
                fl.append(np.linalg.norm(l.F1))

            for n in c.nodes:
                fc.append(np.linalg.norm(n.F))

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


class Configuration:
    def __init__(self, num, dumpdir=None, dt=0.01, nmax=3000, qmin=0.001, d0_0=1, force_limit=15., p_add=1.,
                 p_del=0.2, chkx=True, anis=0.0, d0max=2., unrestricted=True, distlimit=3):
        self.nodes = []
        self.links = []
        self.nNodes = 0  # bookkeeping for state vector position (I don't get that)
        self.dumpdir = dumpdir  # string for directory to save configs
        # parameters for mechanical equilibration
        self.dt = dt
        self.nmax = nmax
        self.qmin = qmin
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
        self.N = num
        self.ones = np.ones(self.N)

        """description of nodes"""
        self.nodes = np.zeros((self.N, 2, 3))           # r and phi of nodes
        self.Fnode = np.zeros((self.N, 3))              # total force on node
        self.Mnode = np.zeros((self.N, 3))              # total torsion on node

        """description of links"""
        # islink[i, j] is True if nodes i and j are connected via link
        self.islink = np.full((self.N, self.N), False)

        self.e = np.zeros((self.N, self.N, 3))          # direction connecting nodes (a.k.a. "actual direction")
        self.d = np.zeros((self.N, self.N))             # distance between nodes (a.k.a. "actual distance")
        self.k = np.zeros((self.N, self.N))             # spring constant between nodes
        self.bend = 10. * np.ones(self.N, self.N)       # bending rigidity
        self.twist = 1. * np.ones(self.N, self.N)       # torsion spring constant
        self.d0 = np.zeros((self.N, self.N))            # equilibrium distance between nodes,
        self.t = np.zeros((self.N, self.N, 3))          # tangent vector of link at node (a.k.a. "preferred direction")
        self.norm = np.zeros((self.N, self.N, 3))       # normal vector of link at node
        self.Mlink = np.zeros((self.N, self.N, 3))      # Torsion from link on node
        self.Flink = np.zeros((self.N, self.N, 3))      # Force from link on node

    def addlink(self, n1, n2, t1=None, t2=None, d0=None, k=None, bend=None, twist=None, n=None, norm1=None, norm2=None):
        if n1 > n2:
            ni = n1
            mi = n2
        else:
            ni = n2
            mi = n1

        self.islink[ni, mi], self.islink[mi, ni] = True, True

        if k is None:
            k = 15.
        self.k[ni, mi], self.k[mi, ni] = k, k  # spring parameter
        if bend is not None:
            self.bend[mi, ni], self.bend[ni, mi] = bend, bend
        if twist is not None:
            self.twist[mi, ni], self.bend[ni, mi] = twist, twist

        if d0 is None:
            d0 = self.d[ni, mi]
        self.d0[ni, mi], self.d0[mi, ni] = d0, d0  # equilibrium distance
        # preferred directions
        if t1 is None:
            self.t[ni, mi] = R(-self.nodes[ni, 1], self.e[ni, mi])
        else:
            self.t[ni, mi] = getNormvec(t1)
        if t2 is None:
            self.t[mi, ni] = R(-self.nodes[mi, 1], -self.e[ni, mi])
        else:
            self.t[mi, ni] = getNormvec(t2)
        if n is None:
            n, q = getNormtoo(np.cross(self.e[ni, mi], ez))  # n is perpendicular to e
            # n is perpendicular to z (l is in the x-y plane)
            if q < 1e-5:
                n = getNormvec(np.cross(self.e[ni, mi], ex))  # e || ez   =>	n is perpendicular to x
        if norm1 is None:
            norm1 = n
        self.norm[ni, mi] = norm1
        if norm2 is None:
            norm2 = n
        self.norm[mi, ni] = norm2

    def removelink(self, n1, n2):
        if n1 > n2:
            ni = n1
            mi = n2
        else:
            ni = n2
            mi = n1
        self.islink[ni, mi], self.islink[mi, ni] = False, False
        self.Flink[ni, mi], self.Flink[mi, ni], self.Mlink[ni, mi], self.Mlink[mi, ni] = null, null, null, null

    def updateDists(self, X):
        dX = X - np.tile(np.expand_dims(X, 1), (1, X.shape[0], 1))
        self.d = np.linalg.norm(dX, axis=2)
        inds = np.where(self.d == 0)
        self.d[inds] = np.infty
        self.e = dX / self.d[..., None]

    def compactStuffINeed(self):
        #get only those parts of the big arrays that are actually needed

        t = self.norm[self.islink]
        norm  = self.norm[self.islink]
        normT = norm.transpose(norm, axes = (1, 0, 2))
        bend = self.bend[self.islink]
        twist = self.twist[self.islink]
        k = self.k[self.islink]
        d0 = self.d0[self.islink]
        nodeinds = np.where(self.islink is True)

        return t, norm, normT, bend, twist, k, d0, nodeinds

    def updateLinkForces(self, X, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):

        Nodes = X[Nodeinds[0]]
        NodesT = X[Nodeinds[1]]
        E = self.e[self.islink]
        D = self.d[self.islink]

        Norm2 = np.cross(R(Nodes[1], Norm), E)

        self.Mlink[self.islink] = Bend * np.cross(R(Nodes[1], T), E) + \
                              Twist * np.cross(R(Nodes[1], Norm), R(NodesT[1], NormT))  # Eqs. 5, 6

        M = self.Mlink + np.transpose(self.Mlink, axes = (1, 0, 2))

        self.Flink[self.islink] = np.dot(M[self.islink], Norm2 - Norm) / D[:, None] + \
                              K * (D  - D0)   # Eqs. 13, 14, 15

    def getForces(self, t, X, norm, normT, bend, twist, k, d0, nodeinds):
        self.updateDists(X)
        self.updateLinkForces(X, t, norm, normT, bend, twist, k, d0, nodeinds)
        self.Fnode = np.sum(self.Flink, axis = 0)
        self.Mnode = np.sum(self.Flink, axis = 0)

        return np.transpose(np.array([self.Fnode, self.Mnode]), axes = (1, 0, 2))

    def mechEquilibrium(self):
        x = self.nodes.copy()
        h = self.dt
        steps = self.nmax
        t, norm, normT, bend, twist, k, d0, nodeinds = self.compactStuffINeed()
        for i in range(self.nmax):
            k1 = self.getForces(x, t, norm, normT, bend, twist, k, d0, nodeinds)
            Q = np.dot(k1, k1) / len(self.nodes)
            # print i, Q, max([x[n.r][0] for n in self.nodes])
            if i % 100 == 0 and self.dumpdir is not None:
                dumpData('s/%05d.pickle' % (self.dumpdir, i))
            if Q < self.qmin:
                steps = i + 1
                break
            k1 *= h
            k2 = h * self.getForces(x + k1 / 2, t, norm, normT, bend, twist, k, d0, nodeinds)
            k3 = h * self.getForces(x + k2 / 2, t, norm, normT, bend, twist, k, d0, nodeinds)
            k4 = h * self.getForces(x + k3, t, norm, normT, bend, twist, k, d0, nodeinds)
            x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.
        self.nodes = x
        return steps * h

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
            f = np.linalg.norm(l.F1)  # magnitude of force transmitted
            if hasattr(l, "limit"):
                p = exp(f / l.limit)
            else:
                p = exp(f / self.force_limit)
            if self.anis > 1e-5:  # horizontal links are more stable
                p = p * (1 - self.anis + self.anis * pow((l.n1.getR() - l.n2.getR())[1] / l.d, 2))  # what?
            to_del.append((l, p))
        return to_del

    def tryLink(self, n, m):
        if n.findLinkTo(m) is not None:
            return -1, null
        for l in self.links:
            if intersect(self.x[n.r], self.x[m.r], self.x[l.n1.r], self.x[l.n2.r]):
                return -1, null  # false
        e, d = getNormtoo(self.x[n.r] - self.x[m.r])
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
        configs, ts = [copy.deepcopy(self)], [0.]
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
