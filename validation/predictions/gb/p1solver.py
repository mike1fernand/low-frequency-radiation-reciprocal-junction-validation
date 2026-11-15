"""Independent frequency-domain solver for the Paper One line-bulk junction.

Model (Paper One, Sec. II): scalar line (m_q, K_q) along x; simple-cubic vector bulk (m_b = 1,
axial K_axis, face-diagonal K_diag); segment s of the junction couples line sites q_{s-1}, q_s to the
four corners at (x_s, +-1/2, +-1/2) by paired links of stiffness g_s.

Bulk Green matrix: exact in k_x (generalized Bloch roots of the quadratic eigenproblem
M(w) = A1 w^2 + (A0 - z) w + A1^*, w = exp(i k_x)), midpoint quadrature in k_perp.
Only the centred A1 sector (one axial X and one radial R pattern per segment) is needed.

Units: a = m_b = K_diag = 1 unless stated.
"""
import numpy as np


# ----------------------------------------------------------------------------- bulk blocks
def bulk_blocks(ky, kz, Ka=1.0, Kd=1.0):
    """A0 (real) and A1 (complex) with D(k) = A0 + A1 e^{ik_x} + conj(A1) e^{-ik_x}."""
    cy, cz, sy, sz = np.cos(ky), np.cos(kz), np.sin(ky), np.sin(kz)
    sh = np.shape(ky)
    A1 = np.zeros(sh + (3, 3), complex)
    A1[..., 0, 0] = -(Ka + Kd * (cy + cz))
    A1[..., 0, 1] = A1[..., 1, 0] = -1j * Kd * sy
    A1[..., 0, 2] = A1[..., 2, 0] = -1j * Kd * sz
    A1[..., 1, 1] = -Kd * cy
    A1[..., 2, 2] = -Kd * cz
    A0 = np.zeros(sh + (3, 3))
    A0[..., 0, 0] = 2 * Ka + 4 * Kd
    A0[..., 1, 1] = 2 * Kd + 2 * Ka * (1 - cy)
    A0[..., 2, 2] = 2 * Kd + 2 * Ka * (1 - cz)
    fp = Kd * (1 - np.cos(ky + kz))      # (0,1,1)/sqrt2 : 2Kd(1-cos)(1/2)[[0,0,0],[0,1,1],[0,1,1]]
    fm = Kd * (1 - np.cos(ky - kz))      # (0,1,-1)/sqrt2
    A0[..., 1, 1] += fp + fm
    A0[..., 2, 2] += fp + fm
    A0[..., 1, 2] = A0[..., 2, 1] = fp - fm
    return A0, A1


def D_full(k, Ka=1.0, Kd=1.0):
    """Real symmetric dynamical matrix D(k) (m_b = 1) for k of shape (...,3)."""
    bonds = [((1, 0, 0), Ka), ((0, 1, 0), Ka), ((0, 0, 1), Ka),
             ((1, 1, 0), Kd), ((1, -1, 0), Kd), ((1, 0, 1), Kd), ((1, 0, -1), Kd),
             ((0, 1, 1), Kd), ((0, 1, -1), Kd)]
    k = np.asarray(k, float)
    M = np.zeros(k.shape[:-1] + (3, 3))
    for d, K in bonds:
        d = np.array(d, float); e = d / np.linalg.norm(d)
        M += 2 * K * (1 - np.cos(k @ d))[..., None, None] * np.outer(e, e)
    return M


def dD_dkx(k, Ka=1.0, Kd=1.0):
    bonds = [((1, 0, 0), Ka), ((1, 1, 0), Kd), ((1, -1, 0), Kd), ((1, 0, 1), Kd), ((1, 0, -1), Kd)]
    k = np.asarray(k, float)
    M = np.zeros(k.shape[:-1] + (3, 3))
    for d, K in bonds:
        d = np.array(d, float); e = d / np.linalg.norm(d)
        M += 2 * K * np.sin(k @ d)[..., None, None] * np.outer(e, e)   # d_x = 1
    return M


# ----------------------------------------------------------------------------- projections
def phi_vectors(ky, kz):
    """Transverse A1 form factors (X, R patterns), defined with exp(-i k_perp . rho_c)."""
    sh = np.shape(ky)
    PX = np.zeros(sh + (3,), complex)
    PR = np.zeros(sh + (3,), complex)
    PX[..., 0] = 2 * np.cos(ky / 2) * np.cos(kz / 2)
    PR[..., 1] = -np.sqrt(2) * 1j * np.sin(ky / 2) * np.cos(kz / 2)
    PR[..., 2] = -np.sqrt(2) * 1j * np.cos(ky / 2) * np.sin(kz / 2)
    return np.stack([PX, PR], -1)          # (..., 3, 2)


# ----------------------------------------------------------------------------- Green blocks
def green_blocks(omega, nmax, Nperp, shift=(0.5, 0.5), eps=0.0, Ka=1.0, Kd=1.0,
                 onshell=True, chunk=8192, tol_unit=1e-6, quadrant=False):
    """Projected bulk Green blocks for separations n = 0..nmax-1 (units of a) in the A1 basis.

    Returns dict with
      'G'   : (nmax, 2, 2) complex, g0 block for n >= 0 (block for -n is its transpose)
      'ImT','ImL' : (2*nmax-1, 2, 2) complex on-shell branch-resolved Hermitian-imaginary blocks for
              n = -(nmax-1)..(nmax-1) (only if eps == 0 and onshell)
    """
    z = omega**2 + 1j * eps
    h = 2 * np.pi / Nperp
    if quadrant:
        # integrand is even in k_y and in k_z (mirror symmetries of lattice and A1 patterns)
        assert Nperp % 2 == 0 and shift == (0.5, 0.5)
        ks = (np.arange(Nperp // 2) + 0.5) * h
        ky_all, kz_all = np.meshgrid(ks, ks, indexing="ij")
    else:
        ks = -np.pi + (np.arange(Nperp) + shift[0]) * h
        ky_all, kz_all = np.meshgrid(ks, -np.pi + (np.arange(Nperp) + shift[1]) * h, indexing="ij")
    ky_all = ky_all.ravel(); kz_all = kz_all.ravel()
    nvec = np.arange(nmax)
    nsym = np.arange(-(nmax - 1), nmax)
    G = np.zeros((nmax, 2, 2), complex)
    ImB = np.zeros((3, 2 * nmax - 1, 2, 2), complex)   # per branch index 0,1,2
    stats = dict(n_prop=0, n_fallback=0, min_absv=np.inf)
    for i0 in range(0, ky_all.size, chunk):
        ky = ky_all[i0:i0 + chunk]; kz = kz_all[i0:i0 + chunk]
        m = ky.size
        A0, A1 = bulk_blocks(ky, kz, Ka, Kd)
        A1c = np.conj(A1)
        Ph = phi_vectors(ky, kz)                       # (m,3,2)
        # companion matrices
        A1inv = np.linalg.inv(A1)
        Cm = np.zeros((m, 6, 6), complex)
        Cm[:, 0:3, 3:6] = np.eye(3)
        Cm[:, 3:6, 0:3] = -A1inv @ A1c
        Cm[:, 3:6, 3:6] = -A1inv @ (A0 - z * np.eye(3))
        w, X = np.linalg.eig(Cm)                       # (m,6), (m,6,6)
        phi = X[:, 0:3, :]                             # (m,3,6)
        Mw = A1[:, None] * (w**2)[:, :, None, None] + (A0 - z * np.eye(3))[:, None] * w[:, :, None, None] + A1c[:, None]
        resid = np.linalg.norm(np.einsum('mjcd,mdj->mjc', Mw, phi), axis=2) / \
            (np.linalg.norm(phi, axis=1) * (1 + np.abs(w)**2))
        stats['max_resid'] = max(stats.get('max_resid', 0.0), float(resid.max()))
        absw = np.abs(w)
        Mp = 2 * A1[:, None] * w[:, :, None, None] + (A0 - z * np.eye(3))[:, None]   # (m,6,3,3)
        den = np.einsum('mcj,mjcd,mdj->mj', phi, Mp, phi)
        # residue projected: Phi^dag phi phi^T Phi / den
        left = np.einsum('mca,mcj->maj', np.conj(Ph), phi)      # (m,2,6)
        right = np.einsum('mcj,mcb->mjb', phi, Ph)              # (m,6,2)
        inside = absw < 1 - tol_unit
        unit = np.abs(absw - 1) <= tol_unit
        incl = inside.copy()
        if np.any(unit):
            mi, ji = np.nonzero(unit)
            kx = np.angle(w[mi, ji])
            kk = np.stack([kx, ky[mi], kz[mi]], -1)
            Dk = D_full(kk, Ka, Kd)
            ev, evec = np.linalg.eigh(Dk)
            b = np.argmin(np.abs(ev - omega**2), axis=1)
            e = evec[np.arange(len(mi)), :, b]            # (p,3) real
            dD = dD_dkx(kk, Ka, Kd)
            vz = np.einsum('pc,pcd,pd->p', e, dD, e)
            stats['n_prop'] += len(mi)
            stats['min_absv'] = min(stats['min_absv'], np.min(np.abs(vz)))
            incl[mi, ji] = vz > 0
            if onshell and eps == 0.0:
                Pp = Ph[mi]                                  # (p,3,2)
                pe = np.einsum('pca,pc->pa', np.conj(Pp), e)  # Phi^dag e  (p,2)
                eP = np.einsum('pc,pcb->pb', e, Pp)           # e^T Phi   (p,2)
                wgt = 1.0 / (2 * np.abs(vz))
                ph = np.exp(1j * np.outer(kx, nsym))          # (p, 2nmax-1)
                contrib = wgt[:, None, None, None] * ph[:, :, None, None] * \
                    (pe[:, None, :, None] * eP[:, None, None, :])
                for bi in range(3):
                    sel = b == bi
                    if np.any(sel):
                        ImB[bi] += contrib[sel].sum(0)
        wn = np.where(incl[:, :, None], w[:, :, None] ** nvec[None, None, :], 0)   # (m,6,nmax)
        G += np.einsum('maj,mjn,mj,mjb->nab', left, wn, 1.0 / den, right)
    N2 = ky_all.size
    if quadrant:
        # the full-grid average equals the quadrant average for an even integrand;
        # the on-shell sums keep only even parts after the (n,-n) and ImH symmetrisation, which the
        # quadrant average reproduces for the real-symmetric assembled matrices
        pass
    out = dict(G=G / N2, stats=stats, nmax=nmax)
    if onshell and eps == 0.0:
        ImB /= N2
        out['ImT'] = ImB[0] + ImB[1]
        out['ImL'] = ImB[2]
    return out


def green_blocks_direct(omega, eps, nmax, N, Ka=1.0, Kd=1.0):
    """Brute-force check at complex frequency: full 3D grid sum (smooth for large eps)."""
    z = omega**2 + 1j * eps
    k1 = -np.pi + (np.arange(N) + 0.5) * 2 * np.pi / N
    G = np.zeros((nmax, 2, 2), complex)
    for kx in k1:
        ky, kz = np.meshgrid(k1, k1, indexing='ij')
        kk = np.stack([np.full_like(ky, kx), ky, kz], -1)
        Dk = D_full(kk, Ka, Kd)
        R = np.linalg.inv(Dk - z * np.eye(3))
        Ph = phi_vectors(ky, kz)
        blk = np.einsum('yzca,yzcd,yzdb->ab', np.conj(Ph), R[..., :, :], Ph) if False else \
            np.einsum('yzca,yzcd,yzdb->yzab', np.conj(Ph), R, Ph).sum((0, 1))
        for n in range(nmax):
            G[n] += np.exp(1j * kx * n) * blk
    return G / N**3


# ----------------------------------------------------------------------------- assembly
def assemble(blocks_G, Nc):
    """Full 2Nc x 2Nc complex symmetric g0 from n >= 0 blocks (segment s at x = s)."""
    g = np.zeros((2 * Nc, 2 * Nc), complex)
    for s in range(Nc):
        for t in range(Nc):
            n = s - t
            B = blocks_G[n] if n >= 0 else blocks_G[-n].T
            g[2 * s:2 * s + 2, 2 * t:2 * t + 2] = B
    return g


def assemble_im(blocks_im, Nc, nmax):
    g = np.zeros((2 * Nc, 2 * Nc), complex)
    for s in range(Nc):
        for t in range(Nc):
            g[2 * s:2 * s + 2, 2 * t:2 * t + 2] = blocks_im[(s - t) + nmax - 1]
    return g


def coupling_mats(gprof):
    """U (2Nc x 2Nc), V (2Nc x (Nc+1)), Kself ((Nc+1)x(Nc+1)) in the A1 basis (Paper One Eq. 45)."""
    Nc = len(gprof)
    U = np.zeros((2 * Nc, 2 * Nc)); V = np.zeros((2 * Nc, Nc + 1)); Ks = np.zeros((Nc + 1, Nc + 1))
    for s, g in enumerate(gprof):
        U[2 * s:2 * s + 2, 2 * s:2 * s + 2] = (2 * g / 3) * np.diag([1.0, 2.0])
        V[2 * s:2 * s + 2, s:s + 2] = (2 * g / 3) * np.array([[-1.0, -1.0], [-np.sqrt(2), np.sqrt(2)]])
        Ks[s:s + 2, s:s + 2] += (4 * g / 3) * np.eye(2)
    return U, V, Ks


def solve_junction(gprof, omega, g0, ImT, ImL, Kq=1.69, mq=1.0):
    """Caroli / Fisher-Lee solution. Returns dict of R, Tq, TqT, TqL, closure and the dressed source."""
    Nc = len(gprof)
    U, V, Ks = coupling_mats(gprof)
    kap = 2 * np.arcsin(omega * np.sqrt(mq / Kq) / 2)
    sig = -Kq * np.exp(1j * kap); gam = 2 * Kq * np.sin(kap)
    nd = Nc + 1
    KD = 2 * Kq * np.eye(nd) - Kq * (np.eye(nd, k=1) + np.eye(nd, k=-1))
    I2 = np.eye(2 * Nc)
    gc = np.linalg.solve(I2 + g0 @ U, g0)
    Sb = -V.T @ gc @ V
    SL = np.zeros((nd, nd), complex); SL[0, 0] = sig
    SR = np.zeros((nd, nd), complex); SR[-1, -1] = sig
    GD = np.linalg.inv(omega**2 * mq * np.eye(nd) - KD - Ks - SL - SR - Sb)
    chi = GD[:, 0] * np.sqrt(gam)
    Tq = gam**2 * abs(GD[-1, 0])**2
    r = 1 - 1j * gam * GD[0, 0]
    R = abs(r)**2
    eta = np.linalg.solve(I2 + U @ g0, V @ chi)
    TqT = 2 * np.real(np.conj(eta) @ ImT @ eta)
    TqL = 2 * np.real(np.conj(eta) @ ImL @ eta)
    Imtot = (g0 - g0.conj().T) / 2j
    Tb_tot = 2 * np.real(np.conj(eta) @ Imtot @ eta)
    return dict(R=R, Tq=Tq, TqT=TqT, TqL=TqL, Tbulk=TqT + TqL, Tbulk_fromImg0=Tb_tot,
                closure=R + Tq + TqT + TqL - 1, closure_g0=R + Tq + Tb_tot - 1,
                eta=eta, chi=chi, kappa=kap, gam=gam)


def reference_profile(Nc):
    if Nc == 2:
        return np.array([0.1, 0.1])
    g = np.full(Nc, 0.2); g[0] = g[-1] = 0.1
    return g


# ----------------------------------------------------------------------------- on-shell (spherical)
def gradD_dir(k, nh, Ka=1.0, Kd=1.0):
    """n_hat . grad_k D(k) for k, nh of shape (...,3)."""
    bonds = [((1, 0, 0), Ka), ((0, 1, 0), Ka), ((0, 0, 1), Ka),
             ((1, 1, 0), Kd), ((1, -1, 0), Kd), ((1, 0, 1), Kd), ((1, 0, -1), Kd),
             ((0, 1, 1), Kd), ((0, 1, -1), Kd)]
    M = np.zeros(k.shape[:-1] + (3, 3))
    for d, K in bonds:
        d = np.array(d, float); e = d / np.linalg.norm(d)
        M += (2 * K * np.sin(k @ d) * (nh @ d))[..., None, None] * np.outer(e, e)
    return M


def shell_radius(omega, branch, nh, Ka=1.0, Kd=1.0, kmax=np.pi, ncoarse=96, nbis=50):
    """Radius along unit directions nh (P,3) where the sorted branch eigenvalue equals omega^2
    (first crossing from k = 0)."""
    w2 = omega**2
    ks = np.linspace(1e-6, kmax, ncoarse)
    vals = np.stack([np.linalg.eigvalsh(D_full(k * nh, Ka, Kd))[:, branch] for k in ks], 1)  # (P,nc)
    above = vals >= w2
    first = np.argmax(above, axis=1)
    ok = above.any(axis=1) & (first > 0)
    lo = ks[np.maximum(first - 1, 0)]; hi = ks[first]
    for _ in range(nbis):
        mid = 0.5 * (lo + hi)
        f = np.linalg.eigvalsh(D_full(mid[:, None] * nh, Ka, Kd))[:, branch] - w2
        lo = np.where(f < 0, mid, lo); hi = np.where(f < 0, hi, mid)
    r = 0.5 * (lo + hi)
    r[~ok] = np.nan
    return r


def onshell_spherical(omega, nmax, nth=160, nph=320, Ka=1.0, Kd=1.0, branches=(0, 1, 2)):
    """Branch-resolved Hermitian-imaginary projected blocks via the constant-frequency surfaces.
    Returns dict branch -> (2nmax-1, 2, 2) complex blocks for n = -(nmax-1)..(nmax-1)."""
    x, wx = np.polynomial.legendre.leggauss(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    ct = np.repeat(x, nph); st = np.sqrt(1 - ct**2); pp = np.tile(ph, nth)
    W = np.repeat(wx, nph) * (2 * np.pi / nph)
    nh = np.stack([ct, st * np.cos(pp), st * np.sin(pp)], -1)     # x along the line
    nsym = np.arange(-(nmax - 1), nmax)
    out = {}
    for b in branches:
        r = shell_radius(omega, b, nh, Ka, Kd)
        assert np.all(np.isfinite(r)), 'shell not star-shaped / not found for some directions'
        k = r[:, None] * nh
        ev, evec = np.linalg.eigh(D_full(k, Ka, Kd))
        e = evec[:, :, b]
        dirder = np.einsum('pc,pcd,pd->p', e, gradD_dir(k, nh, Ka, Kd), e)   # d(omega_b^2)/dk along ray
        assert np.all(dirder > 0), 'non-monotonic branch along a ray'
        jac = r**2 / dirder
        Ph = phi_vectors(k[:, 1], k[:, 2])
        pe = np.einsum('pca,pc->pa', np.conj(Ph), e)
        eP = np.einsum('pc,pcb->pb', e, Ph)
        phase = np.exp(1j * np.outer(k[:, 0], nsym))
        blk = np.einsum('p,pn,pa,pb->nab', W * jac, phase, pe, eP) * (np.pi / (2 * np.pi)**3)
        out[b] = blk
    return out


def re_part_extrapolated(omega, nmax, eps_list=(0.025, 0.0375, 0.05, 0.075, 0.1, 0.125),
                         Nperp_small=512, Nperp_large=384, deg=3, Ka=1.0, Kd=1.0, verbose=False):
    """Retarded projected blocks by analytic (polynomial in i*eps) extrapolation of complex-frequency
    blocks. Returns c0 (nmax,2,2) complex and per-degree alternatives for an error estimate."""
    data = []
    for eps in eps_list:
        Np = Nperp_small if eps < 0.05 else Nperp_large
        data.append(green_blocks(omega, nmax, Np, eps=eps, onshell=False, Ka=Ka, Kd=Kd)['G'])
        if verbose:
            print('  eps', eps, 'done', flush=True)
    data = np.array(data)                       # (ne, nmax, 2, 2)
    t = 1j * np.array(eps_list)
    res = {}
    for d in (deg - 1, deg, deg + 1):
        V = np.vander(t, d + 1, increasing=True)            # (ne, d+1) complex
        flat = data.reshape(len(eps_list), -1)
        coef, *_ = np.linalg.lstsq(V, flat, rcond=None)
        res[d] = coef[0].reshape(data.shape[1:])
    return res, data


# ----------------------------------------------------------------------------- continuum comparators
def continuum_im_blocks(omega, nmax, cT=1.0, cL=np.sqrt(3), rho=1.0, nth=200, nph=400):
    """Hermitian-imaginary continuum Green blocks (shear S and compressional P) between the corner
    patterns (A1 basis), for n = -(nmax-1)..(nmax-1). Im G = pi/(2pi)^3 int d^3k delta(c^2k^2-w^2) P e^{ik.R}."""
    x, wx = np.polynomial.legendre.leggauss(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    ct = np.repeat(x, nph); st = np.sqrt(1 - ct**2); pp = np.tile(ph, nth)
    W = np.repeat(wx, nph) * (2 * np.pi / nph)
    nh = np.stack([ct, st * np.cos(pp), st * np.sin(pp)], -1)
    nsym = np.arange(-(nmax - 1), nmax)
    out = {}
    for lab, c in (('S', cT), ('P', cL)):
        k = omega / c
        kv = k * nh
        Ph = phi_vectors(kv[:, 1], kv[:, 2])                     # transverse form factors
        if lab == 'P':
            proj = nh[:, :, None] * nh[:, None, :]
        else:
            proj = np.eye(3)[None] - nh[:, :, None] * nh[:, None, :]
        left = np.einsum('pca,pcd->pad', np.conj(Ph), proj)       # Phi^dag P
        blk = np.einsum('p,pn,pad,pdb->nab', W, np.exp(1j * np.outer(kv[:, 0], nsym)), left, Ph)
        out[lab] = blk * (np.pi / (2 * np.pi)**3) * k**2 / (2 * c**2 * k) / rho
    return out


def bare_source(gprof, omega, Kq=1.69, mq=1.0):
    """Caroli-normalised bare (uncoupled) line response and bare core source eta = V chi0."""
    Nc = len(gprof)
    U, V, Ks = coupling_mats(gprof)
    kap = 2 * np.arcsin(omega * np.sqrt(mq / Kq) / 2)
    sig = -Kq * np.exp(1j * kap); gam = 2 * Kq * np.sin(kap)
    nd = Nc + 1
    KD = 2 * Kq * np.eye(nd) - Kq * (np.eye(nd, k=1) + np.eye(nd, k=-1))
    SL = np.zeros((nd, nd), complex); SL[0, 0] = sig
    SR = np.zeros((nd, nd), complex); SR[-1, -1] = sig
    GD0 = np.linalg.inv(omega**2 * mq * np.eye(nd) - KD - SL - SR)
    chi0 = GD0[:, 0] * np.sqrt(gam)
    return V @ chi0, chi0
