#!/usr/bin/python  -u

import itertools
import sys
import warnings
from math import exp, log, sqrt

import numpy as np
import numpy.random as npr
import scipy.linalg
import voronoi_neighbors
from mayavi import mlab

warnings.filterwarnings("error")

null = np.array([0.0, 0.0, 0.0])
ex = np.array([1.0, 0.0, 0.0])
ey = np.array([0.0, 1.0, 0.0])
ez = np.array([0.0, 0.0, 1.0])


def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


def getNormvec(v):
    # returns normalized v
    d = scipy.linalg.norm(v)
    try:
        v = v/d
    except RuntimeWarning:
        pass
    return v


def getNormtoo(v):
    # returns norm of v and normalized v
    d = scipy.linalg.norm(v)
    try:
        v = v/d
    except RuntimeWarning:
        d = 0
    return v, d


def intersect(A, B, C, D):
    mynorm = scipy.linalg.norm(np.array([A - C, A - D, B - C]), axis=1)
    if (mynorm < 0.01).any():
        return False
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def getRotMatArray(Phis):
    Thetas = scipy.linalg.norm(Phis, axis=1)
    ThetaDiv = Thetas.copy()
    ThetaDiv[np.where(Thetas == 0)] = np.inf
    Axes = Phis / ThetaDiv[..., None]
    a = np.cos(Thetas / 2)
    b, c, d = np.transpose(Axes) * np.sin(Thetas / 2)
    RotMat = np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                      [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                      [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])
    return np.transpose(RotMat, axes=(2, 0, 1))


def getRotMat(Phis):
    Axis, Theta = getNormtoo(Phis)
    a = np.cos(Theta / 2)
    b, c, d = Axis * np.sin(Theta / 2)
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                    [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                    [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])


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


def showconfig(configs, links, nodeforces, fl, cmap='viridis', vmaxlinks=5, vmaxcells=5, cbar=False):
    # if figure is None:
    #     fig = mlab.figure(figureindex, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    # else:
    #     fig = figure
    x, y, z = configs.T
    xl, yl, zl = configs[links[..., 0]].T
    rxl, ryl, rzl = (configs[links[..., 1]] - configs[links[..., 0]]).T
    fc = scipy.linalg.norm(nodeforces, axis=1)

    cells = mlab.points3d(x, y, z, fc, scale_factor=1, opacity=0.5, resolution=16, scale_mode='none', vmin=0.,
                          colormap=cmap, vmax=vmaxcells)
    links = mlab.quiver3d(xl, yl, zl, rxl, ryl, rzl, scalars=fl, mode='2ddash', line_width=4., scale_mode='vector',
                          scale_factor=1, colormap=cmap, vmin=0., vmax=vmaxlinks)
    links.glyph.color_mode = "color_by_scalar"
    if cbar:
        mlab.scalarbar(links, nb_labels=2, title='Force on link')
    return cells, links


@mlab.animate(delay=100)
def animateconfigs(Configs, Links, nodeForces, linkForces, ts, figureindex=0, bgcolor=(1, 1, 1),
                     fgcolor=(0, 0, 0), figsize=(1000, 1000), cmap='viridis', cbar=False):

    mlab.figure(figureindex, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    vmaxcells = np.max(scipy.linalg.norm(nodeForces, axis=2))
    vmaxlinks = np.max(linkForces)
    cells, links = showconfig(Configs[0], Links[0], nodeForces[0], linkForces[0],
                              vmaxcells=vmaxcells, vmaxlinks=vmaxlinks)
    text = mlab.title('0.0', height=.9)

    while True:
        for (c, l, nF, fl, t) in zip(Configs[1:], Links[1:], nodeForces[1:], linkForces[1:], ts[1:]):
            x, y, z = c.T
            xl, yl, zl = c[l[..., 0]].T
            rxl, ryl, rzl = (c[l[..., 1]] - c[l[..., 0]]).T
            fc = scipy.linalg.norm(nF, axis=1)

            cells.mlab_source.set(x=x, y=y, z=z, scalars=fc)
            links.mlab_source.reset(x=xl, y=yl, z=zl, u=rxl, v=ryl, w=rzl, scalars=fl)
            text.set(text='{}'.format(round(t, 2)))
            print 'Updating... '
            yield



class Configuration:
    def __init__(self, num, dt=0.01, nmax=3000, qmin=0.001, d0_0=1, force_limit=15., p_add=1.,
                 p_del=0.2, chkx=True, anis=0.0, d0max=2.):
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
        # variables to store cell number and cell positions and angles
        self.N = num
        self.ones = np.ones(self.N)

        """stuff for documentation"""
        self.nodesnap = []
        self.linksnap = []
        self.fnodesnap = []
        self.flinksnap = []
        self.snaptimes = []

        """description of nodes"""
        self.nodes = np.zeros((self.N, 2, 3))  # r and phi of nodes in structure [whichnode][r/phi][x/y/z]
        self.Fnode = np.zeros((self.N, 3))              # total force on node
        self.Mnode = np.zeros((self.N, 3))              # total torsion on node

        """description of links"""
        # islink[i, j] is True if nodes i and j are connected via link
        self.islink = np.full((self.N, self.N), False)  # Describes link at node [0] leading to node [1]

        self.e = np.zeros((self.N, self.N, 3))          # direction connecting nodes (a.k.a. "actual direction")
        self.d = np.zeros((self.N, self.N))             # distance between nodes (a.k.a. "actual distance")
        self.k = np.zeros((self.N, self.N))             # spring constant between nodes
        self.bend = 10. * np.ones((self.N, self.N))       # bending rigidity
        self.twist = 1. * np.ones((self.N, self.N))       # torsion spring constant
        self.d0 = np.zeros((self.N, self.N))            # equilibrium distance between nodes,
        self.t = np.zeros((self.N, self.N, 3))          # tangent vector of link at node (a.k.a. "preferred direction")
        self.norm = np.zeros((self.N, self.N, 3))       # normal vector of link at node
        self.Mlink = np.zeros((self.N, self.N, 3))      # Torsion from link on node
        self.Flink = np.zeros((self.N, self.N, 3))      # Force from link on node

    def addlink(self, n1, n2, t1=None, t2=None, d0=None, k=None, bend=None, twist=None, n=None, norm1=None, norm2=None):
        ni = n1
        mi = n2

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
        RotMat1 = getRotMat(-self.nodes[ni, 1])
        RotMat2 = getRotMat(-self.nodes[mi, 1])
        if t1 is None:
            self.t[ni, mi] = np.dot(RotMat1, self.e[ni, mi])
        else:
            self.t[ni, mi] = getNormvec(t1)
        if t2 is None:
            self.t[mi, ni] = np.dot(RotMat2, self.e[mi, ni])
        else:
            self.t[mi, ni] = getNormvec(t2)
        if n is None:
            n, q = getNormtoo(np.cross(self.e[ni, mi], ez))  # n is perpendicular to e
            # n is perpendicular to z (l is in the x-y plane)
            if q < 1e-5:
                n = getNormvec(np.cross(self.e[ni, mi], ex))  # e || ez   =>	n is perpendicular to x
        if norm1 is None:
            norm1 = np.dot(RotMat1, n)
        self.norm[ni, mi] = norm1
        if norm2 is None:
            norm2 = np.dot(RotMat2, n)
        self.norm[mi, ni] = norm2

    def removelink(self, n1, n2):
        ni = n1
        mi = n2

        self.islink[ni, mi], self.islink[mi, ni] = False, False
        self.Flink[ni, mi], self.Flink[mi, ni], self.Mlink[ni, mi], self.Mlink[mi, ni] = null, null, null, null

    def updateDists(self, X):

        dX = X - np.tile(np.expand_dims(X, 1), (1, X.shape[0], 1))
        self.d = scipy.linalg.norm(dX, axis=2)
        inds = np.where(self.d == 0)
        self.d[inds] = np.infty
        self.e = dX / self.d[..., None]

    def compactStuffINeed(self):
        # get only those parts of the big arrays that are actually needed
        t = self.norm[self.islink]
        norm = self.norm[self.islink]
        normT = np.transpose(self.norm, axes=(1, 0, 2))[self.islink]
        bend = self.bend[self.islink]
        twist = self.twist[self.islink]
        k = self.k[self.islink]
        d0 = self.d0[self.islink]
        nodeinds = np.where(self.islink == True)

        return t, norm, normT, bend, twist, k, d0, nodeinds

    def updateLinkForces(self, X, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
        Nodes = X[Nodeinds[0]]
        NodesT = X[Nodeinds[1]]
        E = self.e[self.islink]
        D = self.d[self.islink]
        RotMat = getRotMatArray(Nodes[:, 1, :])
        RotMatT = getRotMatArray(NodesT[:, 1, :])

        NormRot = np.einsum("ijk, ik -> ij", RotMat, Norm)

        self.Mlink[self.islink] = Bend[:, None] * np.cross(np.einsum("ijk, ik -> ij", RotMat, T), E) + \
                                  Twist[:, None] * np.cross(NormRot, np.einsum("ijk, ik -> ij", RotMatT, NormT))  # Eqs. 5, 6

        M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))

        Norm2Rot = np.cross(NormRot, E)

        self.Flink[self.islink] = (np.einsum("ij, ij -> i", M[self.islink], NormRot) / D)[:, None] * Norm2Rot - \
                                  (np.einsum("ij, ij -> i", M[self.islink], Norm2Rot) / D)[:, None] * NormRot + \
                                  (K * (D - D0))[:, None] * E   # Eqs. 13, 14, 15

    def getForces(self, X, t, norm, normT, bend, twist, k, d0, nodeinds):
        self.updateDists(X[:, 0, :])
        self.updateLinkForces(X, t, norm, normT, bend, twist, k, d0, nodeinds)
        self.Fnode = np.sum(self.Flink, axis=1)
        self.Mnode = np.sum(self.Flink, axis=1)

        return np.transpose(np.array([self.Fnode, self.Mnode]), axes=(1, 0, 2))

    def mechEquilibrium(self):
        x = self.nodes.copy()
        h = self.dt
        steps = 0
        t, norm, normT, bend, twist, k, d0, nodeinds = self.compactStuffINeed()
        for i in range(self.nmax):
            k1 = self.getForces(x, t, norm, normT, bend, twist, k, d0, nodeinds)
            Q = np.einsum("ijk, ijk", k1, k1) / self.N
            if Q < self.qmin:
                print "broke"
                break
            k1 *= h
            k2 = h * self.getForces(x + k1 / 2, t, norm, normT, bend, twist, k, d0, nodeinds)
            k3 = h * self.getForces(x + k2 / 2, t, norm, normT, bend, twist, k, d0, nodeinds)
            k4 = h * self.getForces(x + k3, t, norm, normT, bend, twist, k, d0, nodeinds)
            x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.
            steps += 1
        print "Q = ", Q
        self.nodes = x
        return steps * h

    def getLinkList(self):
        allLinks0, allLinks1 = np.where(self.islink == True)
        return np.array([[allLinks0[i], allLinks1[i]] for i in range(len(allLinks0)) if allLinks0[i] > allLinks1[i]])

    def checkLinkX(self):
        Xs = []
        delete_list = []
        Links = self.getLinkList()
        for l1, l2 in itertools.combinations(Links, 2):
            if intersect(self.nodes[l1[0], 0], self.nodes[l1[1], 0], self.nodes[l2[0], 0], self.nodes[l2[1], 0]):
                Xs.append([l1, l2])
        while len(Xs) > 0:
            Xsflat = np.array(Xs).reshape(2 * len(Xs), 2)
            # u: unique elements in Xsflat (a.k.a. links), count: number of occurences in Xsflat
            u, count = np.unique(Xsflat, axis=0, return_counts=True)
            badlink = u[np.argmax(count)]
            delete_list.append(badlink)
            newXs = []
            for linkpair in Xs:
                if (badlink == linkpair[0]).all() or (badlink == linkpair[1]).all():
                    continue
                newXs.append(linkpair)
            Xs = newXs
        for badlink in delete_list:
            self.removelink(badlink[0], badlink[1])

    def delLinkList(self, linklist):
        to_del = []
        for link in linklist:
            if self.d[link[0], link[1]] < self.d0[link[0], link[1]]:
                continue            # compressed links are stable
            f = scipy.linalg.norm(self.Flink[link[0], link[1]])
            p = exp(f / self.force_limit)
            to_del.append((link, p))
        return to_del

    def tryLink(self, n1, n2, Links):
        if self.islink[n1, n2] is True:
            return -1, null
        for l in Links:
            if intersect(self.nodes[n1, 0], self.nodes[n2, 0], self.nodes[l[0], 0], self.nodes[l[1], 0]):
                return -1, null  # false
        e, d = getNormtoo(self.nodes[n1, 0] - self.nodes[n2, 0])
        if d > self.d0max:
            return -1, null  # false
        return d, e  # true: d>0

    def addLinkList(self, Links):
        to_add = []
        linkcands = np.transpose(np.where(self.d < self.d0max))
        for link in linkcands:
            d, e = self.tryLink(link[0], link[1], Links)
            if d > 1e-5:
                p = (1 - (d / 2.0))
                to_add.append(((link[0], link[1]), p))
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
        dt = -log(npr.random()) / S
        if dt > 1:
            print 'Must adjust d0 variables before the next event!'
            return 1.

        r = S * npr.random()
        if r < s1:  # we will remove a link
            for (l, p) in to_del:
                r = r - p * self.p_del
                if r < 0:
                    self.removelink(l[0], l[1])
                    return dt
        r = r - s1
        if r < s2:  # we will add a link
            for ((n1, n2), p) in to_add:
                r = r - p * self.p_add
                if r < 0:
                    self.addlink(n1, n2)
                    return dt

    def default_update_d0(self, dt):
        self.d0 += 0.2 * (self.d0_0 - self.d0) * dt + 0.05 * (
                   2 * sqrt(dt) * npr.random() - sqrt(dt))              # magic number 0.2 and 0.05??

    def modlink(self):
        if self.chkx:
            self.checkLinkX()
        Links = self.getLinkList()
        to_del = self.delLinkList(Links)
        to_add = self.addLinkList(Links)
        dt = self.pickEvent(to_del, to_add)
        self.default_update_d0(dt)
        return dt

    def makesnap(self, t):
        self.nodesnap.append(self.nodes[:, 0, :].copy())
        self.fnodesnap.append(self.Fnode.copy())
        linkList = self.getLinkList()
        self.linksnap.append(linkList)
        self.flinksnap.append(scipy.linalg.norm(self.Flink[linkList[..., 0], linkList[..., 1]], axis = 1))
        self.snaptimes.append(t)

    def timeevo(self, tmax, record=False):
        t = 0.
        if record:
            self.makesnap(0)
        while t < tmax:
            dt = self.mechEquilibrium()
            print "equil: dt = ", dt
            t += dt
            print "t = ", t
            if record:
                self.makesnap(t)
            dt = self.modlink()
            # update_progress(t / tmax)
            print "modify: dt = ", dt
            t += dt
            print "t = ", t
            if record:
                self.makesnap(t)
            # update_progress(t / tmax)
        if record:
            self.nodesnap = np.array(self.nodesnap)
            self.fnodesnap = np.array(self.fnodesnap)
            return self.nodesnap, self.linksnap, self.fnodesnap, self.flinksnap, self.snaptimes
