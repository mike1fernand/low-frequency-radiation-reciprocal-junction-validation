"""Compare independent Python and C++ reference-packet measurements at shared times.

The Python trajectory contains ten integer-time snapshots at 48 through 57.
The C++ reference CSV has 46 samples over the same interval at spacing 0.2.
Pass a fresh C++ CSV and/or --python-json for newly evolved trajectories.
"""
import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FIELDS = {
    "T_current": "T_current",
    "T_coupling": "T_coupling",
    "T_bulk_open": "T_bulk_open",
    "current": "instantaneous_coupling_current_over_E0",
    "historical_shell": "boundary_shell_kinetic_fraction_width2",
    "drift": "total_energy_relative_drift",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpp_csv", nargs="?", type=Path,
                        default=ROOT / "earlier_v4/evidence/sectionVII_scalar_closure_rows.csv")
    parser.add_argument("--python-json", type=Path,
                        default=ROOT / "v32/results/target_blind_global_td_120x56_dt004_v32.json")
    args = parser.parse_args()
    stored = json.loads(args.python_json.read_text())
    python_rows = stored["snapshots"]
    expected_times = list(range(48, 58))
    if [s["time"] for s in python_rows] != expected_times:
        raise ValueError("The Python record must contain all ten declared integer-time snapshots.")
    with args.cpp_csv.open(newline="") as fp:
        cpp_rows = list(csv.DictReader(fp))
    times = [float(row["time"]) for row in cpp_rows]
    if len(times) != 46 or any(not math.isclose(t, 48 + i / 5, rel_tol=0, abs_tol=1e-10) for i, t in enumerate(times)):
        raise ValueError("The C++ rows must cover the 46 reference times once and in order.")
    errors = {}
    for cpp_key, py_key in FIELDS.items():
        differences = []
        for index, row in enumerate(python_rows):
            x = float(cpp_rows[5*index][cpp_key])
            y = float(row[py_key])
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError(f"Non-finite {cpp_key} at {row['time']}")
            differences.append(abs(x-y))
        errors[cpp_key] = max(differences)
    # Use the upward rounded maximum discrepancy reported in the v4
    # comparison of these particular records, measured in their recorded units.
    recorded_cap = 2.83e-15
    ok = max(errors.values()) <= recorded_cap
    print(json.dumps({"status": "PASS" if ok else "FAIL", "samples": len(python_rows),
                      "times": expected_times, "maximum_absolute_difference": max(errors.values()),
                      "maximum_by_field": errors, "historical_reproduction_cap": recorded_cap}, indent=2))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
