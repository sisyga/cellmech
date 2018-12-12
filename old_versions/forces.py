def updateLinkForces2D_newschick_filter(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]

    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(PHI[Nodeinds[0]]), T)

    self.Mlink[Nodeinds] = Bend[..., None] * np.cross(TNow, E)  # Eq 3

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    M = M[Nodeinds]

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E + np.cross(M, E) / D[:, None]  # Eqs. 10, 13, 14, 15


def updateLinkForces2D_unschick_nofilter(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e
    D = self.d
    # RotMat = np.transpose(np.tile(getRotMatArray(PHI), (self.N, 1, 1, 1)), axes=(1,0,2,3))
    RotMat = getRotMatArray(PHI)
    T = self.t
    Bend = self.bend
    K = self.k
    D0 = self.d0

    # vector perpendicular to E for projection of forces
    BaseNorm = getNormvec(np.cross(E, ez))
    BaseNorm2 = np.cross(BaseNorm, E)
    # TNow = np.einsum("hijk, hik -> hij", RotMat, T)  # rotated version of t to fit current setup
    TNow = np.einsum("hjk, hik -> hij", RotMat, T)  # rotated version of t to fit current setup

    self.Mlink = Bend[..., None] * np.cross(TNow, E)

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))

    torqueforce = - np.einsum("ij, ij -> i", M, BaseNorm2)[:, None] * BaseNorm  # Eq. 13, because Eq. 14 = 0 in 2D
    torqueforce /= D[:, None]

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E + torqueforce  # Eq. 10


def updateLinkForces3D_newschick_filter_twist1(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]
    NodesPhi = PHI[Nodeinds[0]]
    NodesPhiT = PHI[Nodeinds[1]]

    # Projection of the rotation vector onto the axis connecting cells A and B
    twistPhi = E * (np.einsum("ik, ik -> i", E, NodesPhi)[..., None])
    twistPhiT = E * (np.einsum("ik, ik -> i", E, NodesPhiT)[..., None])

    # rotated version of Norm and NormT to fit current setup
    NormNow = np.einsum("ijk, ik -> ij", getRotMatArray(twistPhi), Norm)
    NormTNow = np.einsum("ijk, ik -> ij", getRotMatArray(twistPhiT), NormT)

    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(NodesPhi), T)

    self.Mlink[Nodeinds] = Bend[..., None] * np.cross(TNow, E) + \
                           Twist[..., None] * np.cross(NormNow, NormTNow)  # Eq 5

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    M = M[Nodeinds]

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E + np.cross(M, E) / D[:, None]  # Eqs. 10, 13, 14, 15


def updateLinkForces3D_newschick_filter_simontwist(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]
    NodesPhi = PHI[Nodeinds[0]]
    NodesPhiT = PHI[Nodeinds[1]]

    # rotated version of Norm and NormT to fit current setup
    NormNow = np.einsum("ijk, ik -> ij", getRotMatArray(NodesPhi), Norm)
    NormTNow = np.einsum("ijk, ik -> ij", getRotMatArray(NodesPhiT), NormT)

    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(NodesPhi), T)

    # calculated new vector \bm{\tilde{n}}_{B, l}
    NormTilde = getNormvec(NormTNow - np.einsum("ij, ij -> i", NormTNow, TNow)[:, None] * TNow)

    self.Mlink = np.zeros((self.N, self.N, 3))

    self.Mlink[Nodeinds] += Twist[..., None] * np.cross(NormNow, NormTilde)  # Eq 5

    # print scipy.linalg.norm(self.Mlink, axis=-1) - np.transpose(scipy.linalg.norm(self.Mlink, axis=-1), axes=(1,0))

    self.Mlink[Nodeinds] += Bend[..., None] * np.cross(TNow, E)

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    M = M[Nodeinds]

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E + np.cross(M, E) / D[:, None]  # Eqs. 10, 13, 14, 15


def updateLinkForces3D_schick_filter_newtwist2(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]

    twistPHI = self.e * (np.einsum("ijk, ik -> ij", self.e, PHI)[..., None])

    NodesPhi = twistPHI[Nodeinds]
    NodesPhiT = np.transpose(twistPHI, axes=(1, 0, 2))[Nodeinds]

    RotMat = getRotMatArray(NodesPhi)
    RotMatT = getRotMatArray(NodesPhiT)

    NormNow = np.einsum("ijk, ik -> ij", RotMat, Norm)  # rotated version of norm to fit current setup
    NormTNow = np.einsum("ijk, ik -> ij", RotMatT, NormT)
    # vector perpendicular to E for projection of forces
    BaseNorm = getNormvec(np.cross(E, ez[None, :]))
    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(PHI[Nodeinds[0]]), T)

    self.Mlink[self.islink] = Twist[:, None] * np.cross(NormNow, NormTNow)
    # print self.Mlink
    # print self.Mlink + np.transpose(self.Mlink, axes=(1,0,2))

    self.Mlink[self.islink] += Bend[..., None] * np.cross(TNow, E)  # + \
    # Twist[:, None] * np.cross(NormNow, NormTNow)  # Eqs. 5, 6

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    omega = np.cross(self.e, M)
    omega = omega[Nodeinds]

    torqueforce = BaseNorm * (np.einsum("ij, ij -> i", BaseNorm, omega)[:, None])
    torqueforce /= ma.array(D[:, None])
    torqueforce = ma.getdata(ma.filled(torqueforce, 0))

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E - torqueforce  # Eqs. 13, 14, 15


def updateLinkForces3D_schick_filter_newtwist1(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]
    NodesPhi = PHI[Nodeinds[0]]
    NodesPhiT = PHI[Nodeinds[1]]

    twistPhi = E * (np.einsum("ik, ik -> i", E, NodesPhi)[..., None])
    twistPhiT = E * (np.einsum("ik, ik -> i", E, NodesPhiT)[..., None])

    RotMat = getRotMatArray(twistPhi)
    RotMatT = getRotMatArray(twistPhiT)

    NormNow = np.einsum("ijk, ik -> ij", RotMat, Norm)  # rotated version of norm to fit current setup
    NormTNow = np.einsum("ijk, ik -> ij", RotMatT, NormT)
    # vector perpendicular to E for projection of forces
    BaseNorm = getNormvec(np.cross(E, ez[None, :]))
    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(NodesPhi), T)

    self.Mlink[self.islink] = Bend[..., None] * np.cross(TNow, E) + \
                              Twist[:, None] * np.cross(NormNow, NormTNow)  # Eqs. 5, 6

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    omega = np.cross(self.e, M)
    omega = omega[Nodeinds]

    torqueforce = BaseNorm * (np.einsum("ij, ij -> i", BaseNorm, omega)[:, None])
    torqueforce /= ma.array(D[:, None])
    torqueforce = ma.getdata(ma.filled(torqueforce, 0))

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E - torqueforce  # Eqs. 13, 14, 15


def updateLinkForces3D_unschick_filter_oldtwist(self, PHI, T, Norm, NormT, Bend, Twist, K, D0, Nodeinds):
    NodesPhi = PHI[Nodeinds[0]]
    NodesPhiT = PHI[Nodeinds[1]]
    E = self.e[Nodeinds]
    D = self.d[Nodeinds]
    RotMat = getRotMatArray(NodesPhi)
    RotMatT = getRotMatArray(NodesPhiT)

    # vector perpendicular to E for projection of forces
    BaseNorm = getNormvec(np.cross(E, ez))
    BaseNorm2 = np.cross(BaseNorm, E)
    # rotated version of t to fit current setup
    TNow = np.einsum("ijk, ik -> ij", getRotMatArray(PHI[Nodeinds[0]]), T)
    # rotated version of norm to fit current setup
    NormNow = np.cross(E, np.einsum("ijk, ik -> ij", RotMat, Norm))
    NormTNow = np.cross(E, np.einsum("ijk, ik -> ij", RotMatT, NormT))

    self.Mlink[self.islink] = Bend[..., None] * np.cross(TNow, E) + \
                              Twist[..., None] * np.cross(NormNow, NormTNow)  # Eqs. 5, 6

    M = self.Mlink + np.transpose(self.Mlink, axes=(1, 0, 2))
    M = M[Nodeinds]

    torqueforce = - np.einsum("ij, ij -> i", M, BaseNorm2)[:, None] * BaseNorm + \
                  np.einsum("ij, ij -> i", M, BaseNorm)[:, None] * BaseNorm2  # Eq. 13, 14
    torqueforce /= D[:, None]

    self.Flink[Nodeinds] = (K * (D - D0))[..., None] * E + torqueforce  # Eq. 10


def oneequil(self):
    x = self.nodesX.copy()
    phi = self.nodesPhi.copy()
    h = self.dt
    steps = 0
    t, norm, normT, bend, twist, k, d0, nodeinds = self.compactStuffINeed()
    for i in range(self.nmax):
        self.nodesX = x
        self.makesnap(i)
        k1, j1 = self.getForces(x, phi, t, norm, normT, bend, twist, k, d0, nodeinds)
        Q = (np.einsum("ij, ij", k1, k1) + np.einsum("ij, ij", j1, j1)) * self.N_inv
        if Q < self.qmin:
            pass
            break
        if i == 800:
            pass
            # sys.exit()
            # break
        k1, j1 = h * k1, h * j1
        k2, j2 = self.getForces(x + k1 * 0.5, phi + j1 * 0.5, t, norm, normT, bend, twist, k, d0, nodeinds)
        k2, j2 = h * k2, h * j2
        k3, j3 = self.getForces(x + k2 * 0.5, phi + j2 * 0.5, t, norm, normT, bend, twist, k, d0, nodeinds)
        k3, j3 = h * k3, h * j3
        k4, j4 = self.getForces(x + k3, phi + j3, t, norm, normT, bend, twist, k, d0, nodeinds)
        k4, j4 = h * k4, h * j4
        x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.
        phi += (j1 + 2 * j2 + 2 * j3 + j4) / 6.
        steps += 1
    self.nodesPhi = phi
    self.nodesnap = np.array(self.nodesnap)
    self.fnodesnap = np.array(self.fnodesnap)
    return self.nodesnap, self.linksnap, self.fnodesnap, self.flinksnap, self.snaptimes
