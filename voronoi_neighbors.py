import numpy as np
import itertools
from scipy.spatial import Delaunay, delaunay_plot_2d

# describing a plane through point P with normal vector n
# P = (px, py, pz)
# n = (nx, ny, nz)
# for any point X of the plane, X-P will be perpendicular to the n normal vector
# -> (X-P)*n = 0 -> (x-px)*nx+(y-py)*ny+(z-pz)*nz=0

def getNorm(positions):
    cm = np.array([0.0, 0.0, 0.0])  # center of mass
    for p in positions:
        cm += p
    cm /= len(positions)
    M = np.zeros((3, 3))
    for i in positions:
        dr = i - cm
        X, Y, Z = dr[0], dr[1], dr[2]
        M = M + np.array([[X ** 2, X * Y, X * Z],
                          [X * Y, Y ** 2, Y * Z],
                          [X * Z, Y * Z, Z ** 2]])
    evalu, evect = np.linalg.eig(M)
    eind = np.argsort(evalu)
    n = evect[:, eind[0]]
    return n  # a normal vector


def distance(a, b):
    return np.sqrt(sum((a - b) ** 2))


def isNeighborNearby(actual, n, positions, d0max):
    for (p, f) in positions:
        if p is not n and distance(actual, p) < d0max: return True
    return False


def addGhosts(positions, d0max):
    ghosts = []
    for n, surface in positions:
        if not surface: continue  # work only with surface nodes
        normalVector = getNorm([m for (m, s) in positions if s and distance(n, m) < d0max])
        for pos in n + normalVector * 0.3, n - normalVector * 0.3:
            # if not isNeighborNearby(pos, n, positions, d0max):	# isn`t there a node beside pos
            ghosts.append(pos)  # top ghost particle
    return ghosts


#  call as VoronoiNeighbors( [ (cell.x[n.r], n.hasattr("surface") ) for n in cell.nodes ], d0max, is3D )
#  !!!! d0max is ignored in 2D...
#  positions: list of 2-tuples: ((x,y,z),flag)
#  boolean is3D
def VoronoiNeighbors(positions, d0max, vodims=2):
    if vodims == 3:
        # p = [n for (n, f) in positions]  # component [0] is a 3d vector
        # p += addGhosts(positions, d0max) # I don't understand ghost particles, for now, neglect them
        p = positions
    else:
        p = [n[:2] for n in positions]

    # tri: list of interconnectedparticles: [ (a, b, c), (b, c, d), ... ]
    # tri = DelaunayTri(p, joggle=True)
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
    return neighbors  # as there are no ghost particles yet, just return everything
    # if is3D:  # filter list for ghost particles
    #    return [(i, j) for (i, j) in neighbors if i < len(positions) and j < len(positions)]
    # else:
    #    return neighbors


# http://stackoverflow.com/questions/419163/what-does-if-name-main-do

if __name__ == "__main__":
    points = [np.array([10.0, 5.0, 0.1]),
              np.array([-5.0, 8.0, -0.1]),
              np.array([25.0, 13.0, 0.0]),
              np.array([15.0, 18.0, 0.15]),
              np.array([18.0, -17.0, -0.17]),
              np.array([20.0, -8.0, -0.1]),
              np.array([-12.0, -38.0, -0.22])]
    print getNorm(points)
    print distance(points[0], points[1])
