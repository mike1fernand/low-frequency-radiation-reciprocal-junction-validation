"""Alpha-family junction solver (H_c = H_X + alpha*H_R), same-state X/R decomposition, and helpers to load
my independently built Green blocks. Corner energy (g/3)[A^2 + alpha*B^2], A = u_x - qbar,
B = s_y u_y + s_z u_z + dq/2 (Paper One Sec. II / eq. s6-alpha-family)."""
from pathlib import Path
CHECK_ROOT = Path(__file__).resolve().parents[1]
CHECK_REV = Path(__file__).resolve().parent
import sys
import numpy as np
sys.path.insert(0, str(CHECK_ROOT))
from p1solver import assemble, assemble_im, reference_profile


def coupling_mats_alpha(gprof, alpha=1.0, xfac=1.0):
    """U (2Nc x 2Nc), V (2Nc x (Nc+1)), Kself ((Nc+1)x(Nc+1)) for H_c = xfac*H_X + alpha*H_R."""
    Nc = len(gprof)
    U = np.zeros((2 * Nc, 2 * Nc)); V = np.zeros((2 * Nc, Nc + 1)); Ks = np.zeros((Nc + 1, Nc + 1))
    for s, g in enumerate(gprof):
        c = 2 * g / 3
        U[2 * s:2 * s + 2, 2 * s:2 * s + 2] = c * np.diag([xfac, 2.0 * alpha])
        V[2 * s, s:s + 2] = c * np.array([-xfac, -xfac])
        V[2 * s + 1, s:s + 2] = c * np.sqrt(2) * np.array([-alpha, alpha])
        Ks[s:s + 2, s:s + 2] += c * np.array([[xfac + alpha, xfac - alpha], [xfac - alpha, xfac + alpha]])
    return U, V, Ks


def solve_alpha(gprof, omega, g0, ImT, ImL, Kq, alpha=1.0, xfac=1.0, mq=1.0):
    Nc = len(gprof)
    U, V, Ks = coupling_mats_alpha(gprof, alpha, xfac)
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
    R = abs(1 - 1j * gam * GD[0, 0])**2
    M = np.linalg.inv(I2 + U @ g0)
    eta = M @ V @ chi
    VX = V.copy(); VX[1::2] = 0
    VR = V.copy(); VR[0::2] = 0
    etaX = M @ VX @ chi; etaR = M @ VR @ chi
    q = lambda a, B, b: np.real(np.conj(a) @ B @ b)
    TqT = 2 * q(eta, ImT, eta); TqL = 2 * q(eta, ImL, eta)
    return dict(R=R, Tq=Tq, TqT=TqT, TqL=TqL, Tbulk=TqT + TqL, closure=R + Tq + TqT + TqL - 1,
                XX=2 * q(etaX, ImT, etaX), RR=2 * q(etaR, ImT, etaR), XR=4 * q(etaX, ImT, etaR),
                eta=eta, chi=chi, kappa=kap)


def load_g0(omega, Nc=6, Np=2048, deg=2, nth=240, path=str(CHECK_REV)):
    """Assembled g0 (2Nc x 2Nc), ImT, ImL from g0_<omega>.npz (Re extrapolated, Im exact on-shell)."""
    d = np.load(f'{path}/g0_{omega}.npz')
    nmax = int(d['nmax'])
    ImTb = d[f'ImT_{nth}']; ImLb = d[f'ImL_{nth}']          # n = -(nmax-1)..(nmax-1)
    c0 = d[f'c0_{Np}_deg{deg}']
    Gb = c0.real + 1j * (ImTb + ImLb)[nmax - 1:].real
    g0 = assemble(Gb[:Nc], Nc)
    ImT = assemble_im(ImTb[nmax - Nc:nmax + Nc - 1].real, Nc, Nc)
    ImL = assemble_im(ImLb[nmax - Nc:nmax + Nc - 1].real, Nc, Nc)
    return g0, ImT, ImL


def load_g0_ref12(Nc=6):
    """omega = 1.2 blocks from Round 1 (eps scan N=2048 to n=47, spherical Im to n=47)."""
    sys.path.insert(0, str(CHECK_ROOT))
    import os
    cwd = os.getcwd(); os.chdir(str(CHECK_ROOT))
    from build_g0 import g0_blocks_from_scan
    Gb, ImTb, ImLb, nmax = g0_blocks_from_scan('eps_scan_n48.npy', 'onshell_w1.2_n48.npy')
    os.chdir(cwd)
    g0 = assemble(Gb[:Nc], Nc)
    ImT = assemble_im(ImTb[nmax - Nc:nmax + Nc - 1].real, Nc, Nc)
    ImL = assemble_im(ImLb[nmax - Nc:nmax + Nc - 1].real, Nc, Nc)
    return g0, ImT, ImL
