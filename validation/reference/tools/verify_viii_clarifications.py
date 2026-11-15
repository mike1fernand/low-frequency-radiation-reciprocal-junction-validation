"""Recompute the branch partitions used to interpret Section VIII at omega=1.2.

Use --root PATH for the unpacked paper1 directory when running outside tools/.
The two integrations are archived inputs; the two contact formulations run here.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from verify_line_mass import (SURFACE_SHA, TETRA_SHA, TAPER, physical_response,
                              mass_weighted_schur)


def branch_partitions(root):
    methods = (
        ('surface', 'evidence/REFERENCE_KERNEL.npz', SURFACE_SHA,
         ['Im_branch_0', 'Im_branch_1', 'Im_branch_2']),
        ('tetrahedron', 'evidence/line_mass/TETRA_N192.npz', TETRA_SHA,
         ['Im0', 'Im1', 'Im2']),
    )
    results = {}
    for name, file, digest, keys in methods:
        path = root / file
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Unexpected {name} kernel')
        results[name] = {}
        with np.load(path, allow_pickle=False) as kernel:
            green = kernel['G0'][:12, :12]
            shells = [kernel[k][:12, :12] for k in keys]
            for contact, strengths in [('taper', TAPER), ('ramp', np.arange(1., 7.)/7)]:
                args = (1., strengths, green, shells)
                powers = physical_response(*args)
                check = mass_weighted_schur(*args)
                np.testing.assert_allclose(powers, check, rtol=0, atol=2e-14)
                np.testing.assert_allclose(sum(powers), 1, rtol=0, atol=2e-14)
                lower_0, lower_1, upper = powers[2:]
                if min(lower_0, lower_1) <= 0:
                    raise ValueError('Both nominal lower-branch outputs should be positive')
                bulk = lower_0 + lower_1 + upper
                qt = lower_0 + lower_1
                results[name][contact] = dict(lower_0=lower_0, lower_1=lower_1, upper=upper,
                                              bulk=bulk, qT=qt, qT_share_of_bulk=qt/bulk)
                if round(100*bulk, 1) != (13.7 if contact == 'taper' else 80.4):
                    raise ValueError('Printed incident-to-bulk fraction differs')
                if round(100*qt/bulk) != 92:
                    raise ValueError('Printed emitted qT share differs')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    print(json.dumps(dict(status='PASS', omega=1.2, incident_normalized_branch_fractions=branch_partitions(args.root),
                          scope='Nominal contacts on archived bulk matrices. Separate labels at an exact lower-branch degeneracy depend on the basis.'), indent=2))


if __name__ == '__main__':
    main()
