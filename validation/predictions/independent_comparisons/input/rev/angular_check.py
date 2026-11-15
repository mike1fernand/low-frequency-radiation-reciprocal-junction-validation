"""Independent check of the angular-collection numbers: power radiated by the dressed junction source at
omega = 1.2 into a cone of half-angle Theta about +x, with modes selected by GROUP-velocity direction
(and, for comparison, by wave-normal direction). My own code: shell radii from p1solver.shell_radius,
dressed source from my solver (alpha_solver.solve_alpha), coarea weight r^2/(d lambda/dr).
Usage: python3 angular_check.py [nfine ncoarse nphi]
"""
import sys, json, time
import numpy as np
sys.path.insert(0, '/home/claude/p1pred')
from p1solver import D_full, gradD_dir, phi_vectors, shell_radius, reference_profile
from alpha_solver import solve_alpha, load_g0_ref12

nfine = int(sys.argv[1]) if len(sys.argv) > 1 else 1200     # theta midpoints in [0, 60 deg]
ncoarse = int(sys.argv[2]) if len(sys.argv) > 2 else 600    # theta midpoints in [60, 180 deg]
nphi = int(sys.argv[3]) if len(sys.argv) > 3 else 48        # phi midpoints in [0, pi/4]
om = 1.2
g0, ImT, ImL = load_g0_ref12()
g = reference_profile(6)
models = {'full': (1.0, 1.0), 'X-only': (0.0, 1.0), 'R-only': (1.0, 0.0)}
etas = {}
for name, (a, x) in models.items():
    r = solve_alpha(g, om, g0, ImT, ImL, 1.69, alpha=a, xfac=x)
    etas[name] = (r['eta'].reshape(6, 2), r['TqT'], r['TqL'])

t1 = np.deg2rad(60.0)
th = np.concatenate([(np.arange(nfine) + 0.5) * t1 / nfine,
                     t1 + (np.arange(ncoarse) + 0.5) * (np.pi - t1) / ncoarse])
dth = np.concatenate([np.full(nfine, t1 / nfine), np.full(ncoarse, (np.pi - t1) / ncoarse)])
ph = (np.arange(nphi) + 0.5) * (np.pi / 4) / nphi
TH, PH = np.meshgrid(th, ph, indexing='ij')
W = (np.sin(TH) * dth[:, None] * (np.pi / 4) / nphi * 8).ravel()   # solid-angle weight, x8 for C4v sectors
ct = np.cos(TH).ravel(); st = np.sin(TH).ravel(); cp = np.cos(PH).ravel(); sp = np.sin(PH).ravel()
nh = np.stack([ct, st * cp, st * sp], -1)
cones = [5, 10, 20, 30]
res = {m: {'qT': {c: 0.0 for c in cones}, 'qL': {c: 0.0 for c in cones},
           'qT_wn': {c: 0.0 for c in cones}, 'qL_wn': {c: 0.0 for c in cones},
           'tot_qT': 0.0, 'tot_qL': 0.0, 'bwd_k_fwd_v': 0.0} for m in models}
t0 = time.time()
chunk = 60000
for i0 in range(0, len(W), chunk):
    n = nh[i0:i0 + chunk]; w = W[i0:i0 + chunk]
    for b in range(3):
        r = shell_radius(om, b, n)
        assert np.all(np.isfinite(r))
        k = r[:, None] * n
        ev, evec = np.linalg.eigh(D_full(k))
        e = evec[:, :, b]
        dirder = np.einsum('pc,pcd,pd->p', e, gradD_dir(k, n), e)
        assert np.all(dirder > 0)
        grad = np.stack([np.einsum('pc,pcd,pd->p', e, gradD_dir(k, np.tile(u, (len(k), 1))), e)
                         for u in np.eye(3)], -1)                       # grad lambda (= 2 omega v)
        vhat_x = grad[:, 0] / np.linalg.norm(grad, axis=1)
        khat_x = n[:, 0]
        jac = r**2 / dirder
        Ph = phi_vectors(k[:, 1], k[:, 2])                              # (p,3,2)
        eP = np.einsum('pc,pcb->pb', e, Ph)                             # (p,2)
        ph_s = np.exp(-1j * np.outer(k[:, 0], np.arange(6)))            # (p,6)
        for m in models:
            eta = etas[m][0]                                            # (6,2)
            S = np.einsum('ps,pb,sb->p', ph_s, eP, eta)
            dT = 2 * w * jac * np.abs(S)**2 * (np.pi / (2 * np.pi)**3)
            lab = 'qT' if b < 2 else 'qL'
            res[m]['tot_' + lab] += dT.sum()
            res[m]['bwd_k_fwd_v'] += dT[(khat_x < 0) & (vhat_x > np.cos(np.deg2rad(30)))].sum()
            for c in cones:
                cc = np.cos(np.deg2rad(c))
                res[m][lab][c] += dT[vhat_x >= cc].sum()
                res[m][lab + '_wn'][c] += dT[khat_x >= cc].sum()
    print(f'{min(i0 + chunk, len(W))}/{len(W)} directions, {time.time() - t0:.0f}s', flush=True)

print(f'grid: {nfine}+{ncoarse} theta x {nphi} phi (1/8 azimuth)')
for m in models:
    R = res[m]
    print(f"{m}: total qT {R['tot_qT']:.8f} (solver {etas[m][1]:.8f}), qL {R['tot_qL']:.8f} (solver {etas[m][2]:.8f}); "
          f"backward-k modes in 30-deg forward cone: {R['bwd_k_fwd_v']:.2e}")
    for c in cones:
        s = R['qT'][c] / (R['qT'][c] + R['qL'][c])
        swn = R['qT_wn'][c] / (R['qT_wn'][c] + R['qL_wn'][c])
        print(f"   {c:2d} deg: qT {R['qT'][c]:.4e}  qL {R['qL'][c]:.4e}  share {100*s:.3f}%   | wave-normal cone share {100*swn:.2f}%")
json.dump({m: {k: (v if not isinstance(v, dict) else {str(kk): vv for kk, vv in v.items()}) for k, v in R.items()}
           for m, R in res.items()}, open(f'angular_{nfine}_{ncoarse}_{nphi}.json', 'w'), indent=1)
