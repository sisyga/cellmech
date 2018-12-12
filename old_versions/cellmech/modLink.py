#!/usr/bin/python  -u

import os
import math
import random
import argparse
import ConfigParser
import itertools
import collections
import sys
import imp

#
# see http://stackoverflow.com/questions/1046628/importing-python-modules-from-different-working-directory
#
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

import cell
import voronoi_neighbors

dt = 0


def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


def intersect(A, B, C, D):
    if cell.norm(A - C, "mag") < 0.01: return False
    if cell.norm(A - D, "mag") < 0.01: return False
    if cell.norm(B - C, "mag") < 0.01: return False
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


# returns link length, direction vector
def tryLink(n, m, d0max):
    if n.findLinkTo(m) != None: return -1, cell.null  # false
    for l in c.links:
        if intersect(c.x[n.r], c.x[m.r], c.x[l.n1.r], c.x[l.n2.r]):
            return -1, cell.null  # false
    e, d = cell.norm(c.x[n.r] - c.x[m.r])
    if d > d0max: return -1, cell.null  # false
    return d, e  # true: d>0


def addLinkList(anis, is3d, d0max):
    for i, j in voronoi_neighbors.VoronoiNeighbors(
            [n.getR() for n in c.nodes], d0max, is3D=is3d):
        d, e = tryLink(c.nodes[i], c.nodes[j], d0max)
        if d > 1e-5:
            p = (1 - (d / 2.0)) * (1 - anis + anis * e[0] * e[0])
            to_add.append(((c.nodes[i], c.nodes[j]), p))


def checkLinkX():
    Xs = []
    delete_list = []
    for l1, l2 in itertools.combinations(c.links, 2):
        if intersect(l1.n1.getR(), l1.n2.getR(), l2.n1.getR(), l2.n2.getR()):
            Xs.append([l1, l2])
    print "Xs:", Xs
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


def delLinkList(limit, anis):
    for l in c.links:
        if l.d < l.d0: continue  # compressed links are stable
        f = cell.norm(l.F1, "mag")  # magnitude of force transmitted
        if hasattr(l, "limit"):
            p = math.exp(f / l.limit)
        else:
            p = math.exp(f / limit)
        if anis > 1e-5:  # horizontal links are more stable
            p = p * (1 - anis + anis * pow(
                (l.n1.getR() - l.n2.getR())[1] / l.d, 2))
        to_del.append((l, p))


def pickEvent(p_add, p_del):
    global dt

    s1 = 0
    for (l, p) in to_del:
        s1 = s1 + p * p_del
    s2 = 0
    for (q, p) in to_add:
        s2 = s2 + p * p_add

    S = s1 + s2
    if S < 1e-7:
        print "nothing to do!"
        dt = 1
        return
    dt = -math.log(random.uniform(0, 1.0)) / S
    if dt > 1:
        dt = 1  # must adjust d0 variables before the next event
        return

    r = random.uniform(0, S)
    if r < s1:  # we will remove a link
        for (l, p) in to_del:
            r = r - p * p_del
            if r < 0:
                l.n1.dist = 0
                l.n2.dist = 0
                l.remove()
                return
    r = r - s1
    if r < s2:  # we will add a link
        for ((n1, n2), p) in to_add:
            r = r - p * p_add
            if r < 0:
                n1.dist = 0
                n2.dist = 0
                n1.addLinkTo(n2)
                return


def default_update_d0(d0_0):
    for l in c.links:
        l.d0 += 0.2 * (d0_0 - l.d0) * dt + 0.05 * random.uniform(-math.sqrt(dt), math.sqrt(dt))


# log dt and step number
def log(dt):
    # logfile name: logFile
    # logfile path: same path modLink run is started
    try:
        with open(args.logFile, "a") as logfile:  # append
            logfile.write("output file = " + "\t" + args.outputFile + "\t")
            logfile.write("dt = " + "\t" + str(dt))
            logfile.write(os.linesep)
    except IOError as err:
        print('File error: ' + str(err))


def delAttrDist():
    for n in c.nodes:
        if hasattr(n, 'dist'):
            delattr(n, "dist")


#            MAIN

# see http://blog.vwelch.com/2011/04/combining-configparser-and-argparse.html
p1 = argparse.ArgumentParser(add_help=False)
p1.add_argument('-c', dest='conf_file', metavar='FILE',
                help='configuration FILE')
args, remaining_argv = p1.parse_known_args()

defaults = {
    "force_limit": 15.0,
    "p_add": 10.0,
    "p_del": 0.2,
    "chkx": True,
    "anis": 0.0,
    "d0max": 2.0,
    "d0_0": "1.17",
    "update_d0_func": ""
}

if args.conf_file:
    cp = ConfigParser.SafeConfigParser()
    cp.optionxform = str  # to keep upper case configs
    cp.read([args.conf_file])
    defaults.update(dict(cp.items("modLink params")))

p2 = argparse.ArgumentParser(
    parents=[p1],
    description='Change link configuration.',
    formatter_class=argparse.RawDescriptionHelpFormatter)
p2.set_defaults(**defaults)
p2.add_argument('-i', dest='inputFile', metavar='FILE', required=True,
                help='input pickle FILE')
p2.add_argument('-o', dest='outputFile', metavar='FILE', required=True,
                help='output pickle FILE')
p2.add_argument("--force_limit", dest='force_limit', metavar='0.4',
                type=float, help='force limit/link')
p2.add_argument("--p_add", dest='p_add', metavar='1.0',
                type=float, help='probability scaler of link addition')
p2.add_argument("--p_del", dest='p_del', metavar='0.1',
                type=float, help='probability scaler of link deletion')
p2.add_argument("--anis", dest='anis', metavar='1.0',
                type=float, help='link anisotropy ')
p2.add_argument("--chkx", dest='chkx',
                type=bool, help='check for overlapping links')
p2.add_argument("--is3d", dest='is3d', default=False,
                type=bool, help='3D simulation')
p2.add_argument("--d0max", default=2.0, dest='d0max', metavar='2.0',
                type=float, help='Maximal distance bridged. ')
p2.add_argument("--d0_0", default=1.17, dest='d0_0', metavar='1.17',
                type=float, help='Set d0 (equilibrium spring distance) value.')
p2.add_argument('--log', dest='logFile', metavar='FILE', default="log",
                help='LOG FILE')
p2.add_argument('--update_d0_func', dest='update_d0_func',
                help='update_d0 function`s path and name')

args = p2.parse_args(remaining_argv)

c = cell.loadData(args.inputFile)
delAttrDist()

if args.chkx > 0:
    checkLinkX()

to_del = []
delLinkList(args.force_limit, args.anis)
to_add = []
addLinkList(args.anis, args.is3d, args.d0max)
pickEvent(args.p_add, args.p_del)
print "dt:", dt

log(dt)

if args.update_d0_func:
    location = args.update_d0_func  # from conf file
    specName = "d0"
    specFun = "update_d0"
    mod = imp.load_source(specName, location)

    getattr(mod, specFun)(globals())
else:
    default_update_d0(args.d0_0)

print "---"
cell.dumpData(c, args.outputFile)
