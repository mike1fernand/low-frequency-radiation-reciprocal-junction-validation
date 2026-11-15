"""Reconstruct compact scattering inputs from the two original numerical releases."""
import argparse
from hashlib import sha256
import io
from pathlib import Path
from zipfile import ZipFile
import numpy as np

R7_HASH = "e2d3a219de57202001ccd1b20dd985079067c920959c8823368145f73d67137a"
SCREEN_HASH = "5a0f9fd7259be6a899fd904253fda62b60ea973f9c79a5e7799316d801bb9629"
R7 = "paper1-validation-integrated-r7/legacy_rc4/validation/current_numerics/analysis/resolution_execution/"
SCREEN = "Prediction_Calculations_2026-09-26/baseline/Prediction_Calculations_2026-09-26/bandwidth/"
SCREEN_N80 = (1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275)
SCREEN_N128 = (1.1125, 1.2375, 1.3)


def verify_hash(path, expected):
    digest = sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise ValueError(f"Unexpected release digest: {path}")


def arrays(z, name, keys):
    with np.load(io.BytesIO(z.read(name)), allow_pickle=False) as data:
        return [data[key] for key in keys]


def jacksons(resolution, edges, weights):
    """Return the order-256 Jackson principal-value matrix at a squared frequency."""
    lo, hi = edges[0], edges[-1]
    half = (hi-lo)*1.02 / 2
    mid = (lo+hi)/2
    x = ((edges[:-1]+edges[1:])/2 - mid)/half
    y = (resolution-mid)/half
    order = 256
    k = np.arange(order)
    theta = np.pi/(order+1)
    factors = ((order+1-k)*np.cos(k*theta)+np.sin(k*theta)/np.tan(theta))/(order+1)
    cheb_last, cheb = np.ones_like(x), x
    u_last, u = 1.0, 2*y
    hilbert = 2*factors[1]*x
    for k in range(2, order):
        cheb_last, cheb = cheb, 2*x*cheb-cheb_last
        hilbert += 2*factors[k]*u*cheb
        u_last, u = u, 2*y*u-u_last
    blocks = weights @ (hilbert/half)
    n = len(blocks)
    square = np.block([[blocks[i-j] if i >= j else blocks[j-i].conj().T
                        for j in range(n)] for i in range(n)])
    return (square.real + square.real.T)/2


def make(r7_path, screen_path, out):
    verify_hash(r7_path, R7_HASH)
    verify_hash(screen_path, SCREEN_HASH)
    result = {}
    branches = [f"Im_branch_{i}" for i in range(3)]
    with ZipFile(r7_path) as z:
        histograms = [arrays(z, R7+f"surface/histograms/sobol_seed{seed}_p23_nseg10_bins16384.npz",
                            ["edges","mass"]) for seed in (20260920, 20260921)]
        for i in range(49):
            w = (24+i)/40
            n = 10 if i == 24 else 6
            data = arrays(z, R7+f"robustness/kernels/omega{w:.6f}_N80_nseg{n}.npz", branches)
            real = np.mean([jacksons(w*w, *hist) for hist in histograms], axis=0)
            result[f"frequency_{i}"] = np.stack([real[:2*n,:2*n], *data])
        for i, ratio in enumerate([.75,.9,1.,1.1,1.25]):
            result[f"stiffness_{i}"] = np.stack(arrays(z,
                R7+f"robustness/stiffness/kdiag{ratio:g}_N80.npz", ["PV", *branches]))
    with ZipFile(screen_path) as z:
        histograms = [arrays(z, SCREEN+f"inputs/histogram_{i}.npz", ["edges","mass"])
                      for i in (0,1)]
        for w in SCREEN_N80:
            real = np.mean([jacksons(w*w, *hist) for hist in histograms],axis=0)
            data = arrays(z, SCREEN+f"inputs/shell_{w:.6f}_N80.npz", branches)
            result[f"surface_{w:.4f}"] = np.stack([real[:12,:12],*data])
        for i,w in enumerate(SCREEN_N128):
            name = f"fresh/omega{w:.6f}_N128".replace(".", "p")+".npz"
            result[f"surface_{w:.4f}"] = np.stack(arrays(z, SCREEN+name, ["PV_256", *branches]))
            bare,*shells = arrays(z,SCREEN+f"tetra_N64/omega_{i:03d}.npz",
                                  ["G0","Im0","Im1","Im2"])
            np.testing.assert_allclose(bare.imag, sum(shells), rtol=0, atol=1e-13)
            result[f"tetra_{w:.4f}"] = np.stack([bare.real,*shells])
    np.savez_compressed(out, **result)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r7",type=Path,required=True)
    ap.add_argument("--screen",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    make(args.r7, args.screen, args.output)
