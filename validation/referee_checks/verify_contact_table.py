"""Check the five printed contact comparisons at frequency 1.2 in Table 3."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
ROOT = VALIDATION.parent
EVIDENCE = VALIDATION / "reference/evidence/vi_results"
REFERENCE = [Fraction(n, 10) for n in (1, 2, 2, 2, 2, 1)]
PROFILES = {
    "Reference taper": REFERENCE,
    "Doubled taper": [2 * g for g in REFERENCE],
    "Scaled taper": [3 * g for g in REFERENCE],
    "Uniform": [Fraction(1, 2)] * 6,
    "Increasing ramp": [Fraction(j, 7) for j in range(1, 7)],
}


def read_table(source):
    tables = re.findall(r"\\begin\{table\}.*?\\end\{table\}", source, re.DOTALL)
    table, = (table for table in tables if r"\label{tab:s6-coupling-comparison}" in table)
    header, body = table.split(r"\hline")
    headings = header.rsplit("Profile &", 1)[-1].strip().removesuffix(r"\\")
    if [part.strip() for part in headings.split("&")] != [
        "Total stiffness", "Peak link stiffness", "Reflected", "qT conversion"
    ]:
        raise ValueError("Unexpected contact-table column headings")
    rows = []
    for row in body.split(r"\end{tabular}", 1)[0].split(r"\\"):
        if row.strip():
            rows.append([cell.strip() for cell in row.split("&")])
    if [row[0] for row in rows] != list(PROFILES):
        raise ValueError("Unexpected contact-table profiles or row order")
    return rows


def check_values(rows, powers, method):
    for name, total, peak, reflected, converted in rows:
        profile = PROFILES[name]
        if Fraction(total.strip("$")) != 8 * sum(profile) or Fraction(peak.strip("$")) != max(profile):
            raise ValueError(f"Incorrect stiffness total or peak for {name}")
        reflection, conversion = powers[name]
        precision = 2 if name == "Reference taper" else 1
        expected = (f"${100 * reflection:.{precision}f}\\%$", f"${100 * conversion:.1f}\\%$")
        if (reflected, converted) != expected:
            raise ValueError(f"{name}, {method}: printed {(reflected, converted)}, expected {expected}")


def saved_powers():
    scales = json.loads((EVIDENCE / "SCALE_PROBE.json").read_text(encoding="utf-8"))
    records = {row["scale"]: row for row in scales["states"] if row["cache"] == "N192_M256"}
    powers = {}
    for name, factor in (("Reference taper", 1), ("Doubled taper", 2), ("Scaled taper", 3)):
        row = records[factor]
        powers[name] = row["Rq"], row["TqT"]
    comparison = json.loads((HERE / "profile_reversal.json").read_text(encoding="utf-8"))
    for name, key in (("Uniform", "uniform"), ("Increasing ramp", "ramp")):
        row = comparison[key]
        powers[name] = row["reflection"], row["qT"]
    return powers


def replay(rows):
    import numpy as np

    sys.path.insert(0, str(VALIDATION / "reference/tools"))
    from verify_vi_results import loaded_line_response

    # Six complete motifs use the first twelve projected corner coordinates.
    with np.load(EVIDENCE / "REFERENCE_KERNEL.npz", allow_pickle=False) as data:
        surface = (data["G0"], data["Im_qT"], data["Im_qL"])
    with np.load(VALIDATION / "reference_methods/tetra_ladder/N192.npz", allow_pickle=False) as data:
        tetrahedron = (data["G0"], data["Im0"] + data["Im1"], data["Im2"])
    for name, matrices in (("surface N192/M256", surface), ("tetrahedron N192", tetrahedron)):
        core = [matrix[:12, :12] for matrix in matrices]
        powers = {}
        for label, profile in PROFILES.items():
            fractions = loaded_line_response(np.array(profile, dtype=float), *core)
            powers[label] = fractions[0], fractions[2]  # reflected and emitted qT / incident power
        check_values(rows, powers, name)
        print(f"PASS. Five separate coupled solves on the archived {name} kernel match Table 3.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true", help="solve the five contacts with both archived kernels (requires NumPy)")
    args = parser.parse_args()
    text = (ROOT / "manuscript/source/paper1_main.tex").read_text(encoding="utf-8")
    rows = read_table(text)
    check_values(rows, saved_powers(), "saved surface results")
    print("PASS. All five Table 3 rows match the exact stiffness inputs and saved surface results.")
    if args.replay:
        replay(rows)


if __name__ == "__main__":
    main()
