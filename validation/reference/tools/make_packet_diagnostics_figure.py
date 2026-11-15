"""Draw Figure 6(b) from the reference-taper trajectory and six-face audit."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
from matplotlib.figure import Figure
import numpy as np


def draw(evidence, output):
    def read(name):
        return json.loads((evidence / name).read_text())

    with (evidence / "REFERENCE_TRAJECTORY_ROWS.csv").open() as source:
        rows = list(csv.DictReader(source))

    def column(key):
        return np.array([float(row[key]) for row in rows])

    times = column("time")
    modal = read("REFERENCE_MODAL_ANALYSIS.json")["rows"]
    np.testing.assert_allclose(times, [r["time"] for r in modal], rtol=0, atol=1e-10)
    periodic = np.array([r["periodic_total_over_E0"] for r in modal])
    diagnostics = [
        ("coupling", column("T_coupling") / 0.0025, "s", "-"),
        ("current", abs(column("current")) / 0.0005, "^", "--"),
        ("kinetic shell", column("historical_shell") / 0.002, "v", "-."),
        ("periodic/open", abs(periodic / column("T_bulk_open") - 1) / 0.0025, "D", ":"),
        ("drift", abs(column("drift")) / 0.00001, "+", (0, (5, 2, 1, 2, 1, 2))),
    ]
    shell = read("REFERENCE_SHELL_AUDIT.json")
    accepted = read("figure6a_data.json")["accepted_interval"]
    np.testing.assert_allclose(shell["interval"], accepted, rtol=0, atol=1e-10)

    with mpl.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42,
                         "path.simplify": False}):
        # Use the same panel dimensions as 6(a). Reserve the top strip for a two-row key.
        fig = Figure(figsize=(249.276 / 72, 181.655 / 72))
        axes = fig.add_axes([0.135, 0.185, 0.844, 0.635])
        axes.axvspan(*accepted, facecolor="C0", alpha=0.12, edgecolor="none")
        for i, (label, y, marker, style) in enumerate(diagnostics):
            axes.plot(times, y, color=f"C{i}", marker=marker, markevery=(i, 5),
                      markersize=3.3, markerfacecolor="white", markeredgewidth=0.8,
                      linewidth=1.2, linestyle=style, label=label)
        # This is a separate check restricted to samples in the selected interval.
        axes.plot([r["time"] for r in shell["rows"]], [r["union"] / 0.002 for r in shell["rows"]],
                  color="C5", marker="o", markersize=2.5, linewidth=1.2, label="total shell")
        axes.axhline(1, color="0.15", linewidth=1, linestyle="--", label="threshold")
        axes.set(xlim=(48, 57), ylim=(0, 4.8), xticks=[48, 50, 52, 54, 56],
                 yticks=range(5), xlabel=r"Time $\tau$", ylabel="Normalized diagnostic")
        axes.grid(True, alpha=0.22)
        axes.set_axisbelow(True)
        axes.tick_params(pad=2)
        axes.xaxis.labelpad = 2
        axes.yaxis.labelpad = 2

        handles, labels = axes.get_legend_handles_labels()
        order = [0, 1, 2, 4, 3, 5, 6]
        fig.legend([handles[i] for i in order], [labels[i] for i in order],
                   loc="upper center", bbox_to_anchor=(0.5, 0.99), borderaxespad=0,
                   ncols=4, frameon=False, fontsize=8.1, columnspacing=0.62,
                   handlelength=1.25, handletextpad=0.32, labelspacing=0.35, borderpad=0.25)
        fig.savefig(output, format="pdf", metadata={"CreationDate": None, "ModDate": None,
                    "Title": "Reference-taper packet interval diagnostics"})


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=root / "evidence/packet_reference")
    parser.add_argument("--output", type=Path, default=root / "context/figures/paper1_sectionVII_fig6b_selector_controls_v2.pdf")
    args = parser.parse_args()
    draw(args.evidence, args.output)
