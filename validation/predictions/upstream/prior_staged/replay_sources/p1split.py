"""On-shell projected Green blocks split by branch (qT0,qT1,qL) and by direction of energy flow along x
(v_x > 0 forward, v_x < 0 backward). Mixed representation: exact real k_x roots per k_perp, midpoint
quadrature over the k_perp quadrant (integrand even in k_y and k_z)."""
import numpy as np
from p1solver import bulk_blocks, D_full, dD_dkx, phi_vectors

def onshell_split(omega, nmax, Nperp, Ka=1.0, Kd=1.0, chunk=8192, tol_unit=1e-6):
    z = omega**2
    h = 2*np.pi/Nperp
    ks = (np.arange(Nperp//2)+0.5)*h
    KY, KZ = np.meshgrid(ks, ks, indexing='ij'); KY = KY.ravel(); KZ = KZ.ravel()
    nsym = np.arange(-(nmax-1), nmax)
    Im = np.zeros((3, 2, 2*nmax-1, 2, 2), complex)
    minv = np.inf
    for i0 in range(0, KY.size, chunk):
        ky = KY[i0:i0+chunk]; kz = KZ[i0:i0+chunk]; m = ky.size
        A0, A1 = bulk_blocks(ky, kz, Ka, Kd); A1c = np.conj(A1)
        Cm = np.zeros((m, 6, 6), complex); A1inv = np.linalg.inv(A1)
        Cm[:, 0:3, 3:6] = np.eye(3); Cm[:, 3:6, 0:3] = -A1inv@A1c; Cm[:, 3:6, 3:6] = -A1inv@(A0-z*np.eye(3))
        w = np.linalg.eigvals(Cm)
        unit = np.abs(np.abs(w)-1) <= tol_unit
        if not unit.any(): continue
        mi, ji = np.nonzero(unit)
        kx = np.angle(w[mi, ji])
        kk = np.stack([kx, ky[mi], kz[mi]], -1)
        ev, evec = np.linalg.eigh(D_full(kk, Ka, Kd))
        b = np.argmin(np.abs(ev-z), axis=1)
        e = evec[np.arange(len(mi)), :, b]
        vz = np.einsum('pc,pcd,pd->p', e, dD_dkx(kk, Ka, Kd), e)
        minv = min(minv, np.abs(vz).min())
        Ph = phi_vectors(ky[mi], kz[mi])
        pe = np.einsum('pca,pc->pa', np.conj(Ph), e); eP = np.einsum('pc,pcb->pb', e, Ph)
        wgt = 1.0/(2*np.abs(vz))
        ph = np.exp(1j*np.outer(kx, nsym))
        contrib = wgt[:, None, None, None]*ph[:, :, None, None]*(pe[:, None, :, None]*eP[:, None, None, :])
        d = (vz < 0).astype(int)
        for bi in range(3):
            for di in range(2):
                sel = (b == bi) & (d == di)
                if sel.any(): Im[bi, di] += contrib[sel].sum(0)
    return Im/KY.size, minv
