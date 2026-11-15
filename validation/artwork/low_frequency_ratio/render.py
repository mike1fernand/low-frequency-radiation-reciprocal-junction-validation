"""Plot the nominal taper's sampled qT conversion divided by its acoustic limit."""
import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOW = Path("validation/predictions/baseline/Prediction_Calculations_2026-09-26/low_frequency")
SCAN = Path("validation/reference/evidence/vi_results/ACTIVE85.json")
REFERENCE = Path("validation/reference/evidence/vi_results/REFERENCE_METHODS.json")
LOW_FREQUENCIES = (0.05, 0.1, 0.2, 0.4)
CT = 2.0 * 1.3 / (3.0 * math.pi)
FIELDS = ("omega", "T_qT_modal", "normalized", "series")


@dataclass(frozen=True)
class Point:
    omega: float
    qT: float
    series: str

    @property
    def normalized(self):
        return self.qT / (CT * self.omega**2)


def extract_sources(repo):
    def load(relative):
        return json.loads((repo / relative).read_text(encoding="utf-8"))

    record = load(LOW / "inputs.json")
    assert record["units"] == dict(a=1, m_b=1, m_q=1, K_axis=1, K_diag=1, K_q=1.69)
    assert record["profile"] == [0.1, 0.2, 0.2, 0.2, 0.2, 0.1]
    fixed = {}
    for row in load(LOW / "refined/rows.json"):
        if row["s"] == 1 and row["omega"] in LOW_FREQUENCIES:
            # The scale list intentionally repeats s=1 at one frequency.
            assert fixed.setdefault(row["omega"], row["TqT"]) == row["TqT"]
    assert tuple(sorted(fixed)) == LOW_FREQUENCIES
    low = [Point(w, fixed[w], "refined") for w in LOW_FREQUENCIES]

    rows = load(SCAN)["families"]["frequency"]
    assert len(rows) == 49
    for i, row in enumerate(rows):
        assert math.isclose(row["omega"], (24 + i) / 40, abs_tol=1e-14)
        assert row["method"]["surface_grid"] == 80 and row["method"]["PV_order"] == 256
    scan = [Point(row["omega"], row["T_qT_modal"], "scan") for row in rows]

    # The paper's finer reference is a corroboration of the three-digit ratio
    # stated in the text; the marker uses the same N=80 scan as its neighbors.
    selected = Point(1.2, load(REFERENCE)["partition"]["TT_modal"], "N192")
    assert round(scan[24].normalized, 3) == round(selected.normalized, 3) == 0.316
    return low + scan


def read_points():
    with (HERE / "points.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames) != FIELDS:
            raise ValueError("Unexpected point table columns")
        points = []
        for row in reader:
            pt = Point(float(row["omega"]), float(row["T_qT_modal"]), row["series"])
            if not math.isclose(pt.normalized, float(row["normalized"]), rel_tol=2e-15):
                raise ValueError(f"Normalization differs at frequency {pt.omega}")
            points.append(pt)
    return points


def plot(points, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if destination.suffix.lower() != ".pdf":
        raise ValueError("Use a .pdf output for the paper's vector figure")
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        "font.size": 8, "axes.labelsize": 9.5, "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "xtick.major.size": 3, "ytick.major.size": 3,
    })
    fig, ax = plt.subplots(figsize=(3.36, 2.885))
    fig.subplots_adjust(left=0.17, bottom=0.22, right=0.97, top=0.94)
    ax.set(xlim=(0, 1.86), ylim=(0, 1.53),
           xticks=[0.0, 0.4, 0.8, 1.2, 1.6],
           yticks=[0.0, 0.5, 1.0, 1.5],
           xlabel=r"Frequency $\omega\sqrt{m_b/K_{\rm diag}}$",
           ylabel=r"$T_{\rm qT}^{\rm modal}/(C_T\omega^2)$")
    ax.axhline(1.0, color="0.46", linestyle=(0, (4, 2)), linewidth=0.8, zorder=1)
    ax.text(1.81, 1.05, "Low-frequency limit", color="0.30", ha="right", va="bottom",
            fontsize=7.4)
    for series, face, size in [("refined", "white", 4.7), ("scan", "#aa450e", 3.1)]:
        selected = [p for p in points if p.series == series]
        ax.plot([p.omega for p in selected], [p.normalized for p in selected],
                linestyle="none", marker="o", markersize=size, markeredgewidth=0.75,
                markeredgecolor="#aa450e", markerfacecolor=face, zorder=3)
    benchmark = next(p for p in points if p.omega == 1.2)
    ax.plot(benchmark.omega, benchmark.normalized, marker="o",
            markeredgecolor="#252525", markerfacecolor="none",
            markeredgewidth=0.8, markersize=6.9, zorder=4)
    ax.annotate(r"$\omega=1.2$", xy=(benchmark.omega, benchmark.normalized),
                xytext=(1.43, 0.51), ha="center", va="bottom", fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.7, shrinkB=3.3))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(destination, metadata={"Title": "Fixed-contact qT conversion relative to the acoustic law",
                                       "Creator": "Paper 1 figure generator",
                                       "CreationDate": None, "ModDate": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--repository", type=Path,
                        help="Validation repository root, inferred when run in the release")
    args = parser.parse_args()
    points = read_points()
    if args.check_sources:
        repo = args.repository.resolve() if args.repository else HERE.parents[2]
        if points != extract_sources(repo):
            raise ValueError("Point table differs from archived source rows")
        print(f"Verified all {len(points)} points and the finer reference at 1.2.")
    plot(points, args.output)
    print(f"{hashlib.sha256(args.output.read_bytes()).hexdigest()}  {args.output}")


if __name__ == "__main__":
    main()
