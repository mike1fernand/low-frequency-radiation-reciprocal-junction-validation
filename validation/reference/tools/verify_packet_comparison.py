"""Reconstruct the two-contact comparisons in Section VII from retained modal samples."""
import argparse
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def taper_samples(path, window):
    rows = read(path)["rows"]
    return [
        (row["time"], row["qT_min_over_E0"], row["qT_max_over_E0"], row["periodic_total_over_E0"])
        for row in rows if window[0] - 1e-10 <= row["time"] <= window[1] + 1e-10
    ]


def mean_energies(samples, window):
    """Integrate both allocation endpoints and total periodic bulk energy on a common grid."""
    times, low, high, bulk = np.array(samples).T
    expected_ticks = np.arange(round(window[0] / 0.2), round(window[1] / 0.2) + 1)
    np.testing.assert_allclose(times, 0.2 * expected_ticks, rtol=0, atol=1e-10)
    if not np.all((0 <= low) & (low <= high) & (high <= bulk)):
        raise ValueError("Inconsistent modal energy allocation")
    qT, total = np.trapezoid([(low + high) / 2, bulk], times, axis=1) / (times[-1] - times[0])
    return {
        "window": list(window),
        "samples": len(times),
        "mean_qT_per_initial_line": float(qT),
        "mean_bulk_per_initial_line": float(total),
        "mean_qT_divided_by_mean_bulk": float(qT / total),
    }


def stationary_predictions(ramp_dir, spectrum):
    """Use the common incident state to weight each selected contact's own qT coefficient."""
    cfg = read(ramp_dir / "configuration.json")
    q, v = np.loadtxt(ramp_dir / "initial_qv.txt").reshape(2, cfg["N_line"])
    mass, spring = cfg["masses"], cfg["scalar_line_stiffness"]
    initial_energy = (mass * (v @ v) + spring * np.sum(np.diff(q) ** 2)) / 2
    k = 2 * np.pi * np.fft.fftfreq(len(q))
    omega = 2 * np.sqrt(spring / mass) * abs(np.sin(k / 2))
    energy = mass * (abs(np.fft.fft(v)) ** 2 + abs(omega * np.fft.fft(q)) ** 2) / (2 * len(q))
    rows = next(item["rows"] for item in spectrum["families"] if item["label"] == "surface_N160")
    modes = np.array([row["positive_index"] for row in rows])
    evaluated = np.flatnonzero((k > 0) & (omega >= 0.8) & (omega <= 1.6))
    np.testing.assert_array_equal(modes, evaluated)
    np.testing.assert_array_equal([row["negative_index"] for row in rows], len(q) - modes)
    np.testing.assert_allclose([row["omega"] for row in rows], omega[modes], rtol=0, atol=1e-13)
    pair_energy = energy[modes] + energy[-modes]
    coverage = float(pair_energy.sum() / energy.sum())
    fractions = {name: float(pair_energy @ [row["states"][name]["qT"] for row in rows] / initial_energy)
                 for name in ("reference", "ramp")}
    return fractions, coverage


def verify(root):
    reference = root / "evidence" / "packet_reference"
    ramp = root / "evidence" / "ramp_packet"
    taper_window, ramp_window = (52.8, 55.4), (57.4, 64.2)
    cases = {
        "reference": (taper_samples(reference / "REFERENCE_MODAL_ANALYSIS.json", taper_window),
                      taper_samples(reference / "BOX64_MODAL_ANALYSIS.json", taper_window), taper_window),
        "ramp": (read(ramp / "CANDIDATE_PACKET_BOX80_RESULT.json")["branch_samples"],
                 read(ramp / "CONTROL_DT004_BOX96_RESULT.json")["branch_samples"], ramp_window),
    }
    predictions, covered = stationary_predictions(ramp, read(ramp / "PACKET_SPECTRUM.json"))
    np.testing.assert_allclose(100 * covered, 99.968, atol=0.0005, rtol=0)
    reported = {"reference": [13.9, 91.9, 12.76, 12.77, -0.005],
                "ramp": [78.3, 91.4, 71.65, 71.63, 0.015]}
    results = {}
    for name, (primary_rows, larger_rows, window) in cases.items():
        primary = mean_energies(primary_rows, window)
        large = mean_energies(larger_rows, window)
        qT, qT_large = primary["mean_qT_per_initial_line"], large["mean_qT_per_initial_line"]
        predicted = predictions[name]
        gap, large_gap = 100 * (qT - predicted), 100 * (qT_large - predicted)
        displayed = [round(100 * primary["mean_bulk_per_initial_line"], 1),
                     round(100 * primary["mean_qT_divided_by_mean_bulk"], 1),
                     round(100 * qT, 2), round(100 * predicted, 2), round(gap, 3)]
        np.testing.assert_allclose(displayed, reported[name], rtol=0, atol=1e-12)
        results[name] = {"primary": primary, "larger_box": large,
                         "stationary_in_band_per_initial_line": predicted,
                         "direct_minus_stationary_pp": gap,
                         "larger_box_minus_stationary_pp": large_gap,
                         "larger_box_minus_primary_pp": 100 * (qT_large - qT)}
    ramp_result = results["ramp"]
    if not ramp_result["larger_box_minus_stationary_pp"] < 0 < ramp_result["direct_minus_stationary_pp"]:
        raise ValueError("The ramp finite-box comparison no longer brackets the stationary estimate")
    np.testing.assert_allclose(100 * ramp_result["larger_box"]["mean_qT_per_initial_line"],
                               71.627, rtol=0, atol=0.0005)
    return {"status": "PASS", "evaluated_periodic_packet_energy_fraction": covered, "comparisons": results,
            "scope": "Stored modal-energy averages on each primary window and the same-window larger box, "
                     "compared with the evaluated-band prediction. The signed differences describe central "
                     "estimates; these controls supply no combined error interval."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Package root containing evidence/")
    print(json.dumps(verify(parser.parse_args().root), indent=2))
