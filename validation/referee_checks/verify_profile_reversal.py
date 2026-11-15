"""Reproduce the equal-total, opposite-incidence and interior-permutation budgets at 1.2."""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent / "reference"
sys.path.insert(0, str(REFERENCE / "tools"))
from verify_vi_results import loaded_line_response


def main():
    records = json.loads((HERE / "profile_reversal.json").read_text(encoding="utf-8"))
    with np.load(REFERENCE / "evidence/vi_results/REFERENCE_KERNEL.npz", allow_pickle=False) as data:
        mats = [data[key][:12, :12] for key in ("G0", "Im_qT", "Im_qL")]
    ramp = np.arange(1., 7.) / 7
    cases = {
        "equal_taper": np.array([.3, .6, .6, .6, .6, .3]),
        "uniform": np.full(6, .5),
        "ramp": ramp,
        "reversed": ramp[::-1],
        "same_entrance_swap_3_and_5": ramp[[0, 1, 4, 3, 2, 5]],
    }
    computed = {}
    channels = ("reflection", "through", "qT", "qL")
    for name, profile in cases.items():
        np.testing.assert_allclose(profile.sum(), 3., rtol=0, atol=1e-15)
        powers = loaded_line_response(profile, *mats)
        expected = [records[name][key] for key in channels]
        np.testing.assert_allclose(powers, expected, rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(powers.sum(), 1., rtol=0, atol=1e-12)
        computed[name] = dict(zip(channels, map(float, powers)))
    np.testing.assert_allclose(computed["ramp"]["through"], computed["reversed"]["through"],
                               rtol=1e-12, atol=1e-14)
    print(json.dumps(computed, indent=2))
    print("PASS. Incident-power fractions at frequency 1.2 using the supplied reference kernel.")


if __name__ == "__main__":
    main()
