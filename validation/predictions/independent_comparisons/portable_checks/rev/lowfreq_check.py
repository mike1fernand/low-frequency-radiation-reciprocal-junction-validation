"""Independent checks of the report's low-frequency results with my own solver.
1. Static projected compliance G(0) (quadrant midpoint grid, exact in k_x), static device matrix A,
   force balance p^T D A^{-1} e_L = -1, and the line-transmission coefficient 4 Z^2 |[A^{-1}]_{RL}|^2.
2. Full solutions at omega = 0.4 (Green blocks from build_g0_freq.py) for s = 1, 0.1, 0.01 and the
   weak coefficient, against the report's rows.
"""
from pathlib import Path
CHECK_ROOT = Path(__file__).resolve().parents[1]
CHECK_REV = Path(__file__).resolve().parent
import sys
import numpy as np
sys.path.insert(0, str(CHECK_ROOT))
from p1solver import green_blocks, assemble, reference_profile, bare_source
from alpha_solver import coupling_mats_alpha, solve_alpha, load_g0

Nc = 6; Kq = 1.69; Z = 1.3
g = reference_profile(Nc)
# ---- 1. static compliance
Gs = {}
for Np in (512, 1024, 2048):
    Gs[Np] = green_blocks(0.0, Nc, Np, eps=0.0, onshell=False, quadrant=True)['G'].real
    print(f'G(0) N={Np}: n=0 block diag {Gs[Np][0,0,0]:.10f} {Gs[Np][0,1,1]:.10f}', flush=True)
# Richardson in 1/N (error ~ O(1/N) near the acoustic point) using the two finest grids
G0R = 2 * Gs[2048] - Gs[1024]
print('extrapolated n=0 diag', G0R[0, 0, 0], G0R[0, 1, 1], ' change vs N=2048:', np.abs(G0R - Gs[2048]).max())
for lab, Gb in (('N=2048', Gs[2048]), ('Richardson', G0R)):
    G = assemble(Gb, Nc).real
    U, V, Ks = coupling_mats_alpha(g, 1.0)
    nd = Nc + 1
    KD = 2 * Kq * np.eye(nd) - Kq * (np.eye(nd, k=1) + np.eye(nd, k=-1))
    L = KD.copy(); L[0, 0] -= Kq; L[-1, -1] -= Kq
    H = np.linalg.inv(G) + U
    A = L + Ks - V.T @ np.linalg.solve(H, V)
    D = np.linalg.solve(np.eye(2 * Nc) + U @ G, V)
    p = np.zeros(2 * Nc); p[0::2] = 2.0
    eL = np.zeros(nd); eL[0] = 1
    Ainv = np.linalg.inv(A)
    print(f'{lab}: min eig A = {np.linalg.eigvalsh(A).min():.10f}, min eig G = {np.linalg.eigvalsh(G).min():.10f}, '
          f'p^T D A^-1 e_L = {p @ D @ Ainv @ eL:.12f}, T_line/omega^2 -> {4 * Z**2 * abs(Ainv[-1, 0])**2:.10f}')
print('report: min eig static device 0.322886718105273, min eig G 0.105539627711302, T_line coefficient 0.166775444268264')

# ---- 2. omega = 0.4
rep = {1.0: (0.0518080040, 0.0050647968, 0.0437199354, 0.8994072638),
       0.1: (0.0026054524, 0.0002523557, 0.9657715792, 0.0313706127),
       0.01: (0.0000224409, 0.0000021677, 0.9997278790, 0.0002475125)}
for Np, deg in ((2048, 2), (4096, 2), (4096, 3)):
    try:
        g0, ImT, ImL = load_g0(0.4, Np=Np, deg=deg)
    except KeyError:
        continue
    print(f'--- omega = 0.4, Re from N={Np} deg={deg}')
    for s, (tq, tl, tline, R) in rep.items():
        r = solve_alpha(s * g, 0.4, g0, ImT, ImL, Kq, alpha=1.0)
        print(f"  s={s}: TqT={r['TqT']:.10f} ({tq}) TqL={r['TqL']:.10f} ({tl}) Tline={r['Tq']:.10f} ({tline}) "
              f"R={r['R']:.10f} ({R}) closure={r['closure']:.1e}")
    eta0, _ = bare_source(g, 0.4, Kq=Kq)
    print(f"  weak coefficient A_qT = {2 * np.real(np.conj(eta0) @ ImT @ eta0):.8f} (report 0.21986692)")
lim = 2 * Z / (3 * np.pi); print(f'low-frequency limit TqT/w^2 = {lim:.10f}; at w=0.4 report ratio {0.0518080040/0.16:.6f}')
