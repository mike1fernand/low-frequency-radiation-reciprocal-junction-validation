"""Build projected bulk Green blocks g0(omega) for N_c <= nmax with my independent solver:
Re part by complex-frequency extrapolation (quadrant midpoint grid, exact in k_x), Im part exactly on the
constant-frequency shells (spherical parametrization). Saves g0_<omega>.npz.
Usage: python3 build_g0_freq.py omega [nmax] [Nperp ...]
"""
import sys, time
import numpy as np
sys.path.insert(0, '..')
from p1solver import green_blocks, onshell_spherical

om = float(sys.argv[1])
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
Nlist = [int(x) for x in sys.argv[3:]] or [2048]
eps_list = (0.006, 0.008, 0.012, 0.016, 0.024)
out = {'omega': om, 'nmax': nmax, 'eps': np.array(eps_list)}
t0 = time.time()
for Np in Nlist:
    data = np.array([green_blocks(om, nmax, Np, eps=e, onshell=False, quadrant=True)['G'] for e in eps_list])
    out[f'scan_{Np}'] = data
    tt = 1j * np.array(eps_list)
    for deg in (2, 3):
        V = np.vander(tt, deg + 1, increasing=True)
        c, *_ = np.linalg.lstsq(V, data.reshape(len(eps_list), -1), rcond=None)
        out[f'c0_{Np}_deg{deg}'] = c[0].reshape(data.shape[1:])
    print(f'omega={om} Nperp={Np} scan done {time.time()-t0:.0f}s', flush=True)
for nth, nph in ((160, 320), (240, 480)):
    o = onshell_spherical(om, nmax, nth=nth, nph=nph)
    out[f'ImT_{nth}'] = o[0] + o[1]
    out[f'ImL_{nth}'] = o[2]
    print(f'omega={om} spherical {nth}x{nph} done {time.time()-t0:.0f}s', flush=True)
np.savez(f'g0_{om}.npz', **out)
# consistency: extrapolated Im vs exact on-shell Im
Imex = (out['ImT_240'] + out['ImL_240'])[nmax - 1:].real
for Np in Nlist:
    for deg in (2, 3):
        d = np.abs(out[f'c0_{Np}_deg{deg}'].imag - Imex).max()
        print(f'omega={om} N={Np} deg={deg}: max|Im extrap - on-shell| = {d:.2e}  (max|Im| {np.abs(Imex).max():.3f})')
print('sph 160 vs 240 max diff', np.abs(out['ImT_160'] + out['ImL_160'] - out['ImT_240'] - out['ImL_240']).max())
