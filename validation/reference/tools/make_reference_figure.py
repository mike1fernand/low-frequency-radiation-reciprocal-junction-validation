"""Regenerate Figure 4 from the archived reference-partition and solver records."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

ARCHIVED_REFERENCE_SHA256 = "956c9f0dd9d7ae3f0c39cd59f7dc681f43e197575e0cfadecb6f1e37dcdd04b9"


def make_figure(reference, output):
    raw = reference.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ARCHIVED_REFERENCE_SHA256:
        raise ValueError("Figure input differs from the archived reference record")
    data = json.loads(raw)
    partition = data["partition"]
    if abs(sum(partition.values()) - 1) > 2e-11:
        raise ValueError("Reference channels do not balance")
    estimates = {row["label"].split("\n")[0]: row["value"] for row in data["comparisons"]}
    surface = estimates["Surface/PV"]
    if surface != partition["TT_modal"]:
        raise ValueError("Surface point in the two panels differs")
    differences = [(estimates[name] - surface) * 1e5 for name in ("Surface/PV", "Tetrahedron", "GB")]

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.0,
                         "axes.titlesize": 9.1, "axes.labelsize": 8.0,
                         "xtick.labelsize": 7.6, "ytick.labelsize": 7.6,
                         "axes.linewidth": .65, "pdf.fonttype": 42})
    fig = plt.figure(figsize=(7.0, 2.28))
    channels = fig.add_axes([.133, .265, .320, .575])
    methods = fig.add_axes([.616, .265, .360, .575])
    for ax in (channels, methods):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=.65)
        ax.set_axisbelow(True)

    channel_keys = ("Tq", "Rq", "TT_modal", "TL_modal")
    labels = ("Transmitted line", "Reflected line", "qT bulk", "qL bulk")
    percentages = [100 * partition[key] for key in channel_keys]
    colors = ("#ee9030", "#3982ae", "#bb4d45", "#409774")
    bars = channels.barh(range(4), percentages, color=colors, height=.5)
    channels.set_yticks(range(4), labels)
    channels.set_ylim(3.6, -.60)
    channels.set_xlim(0, 100)
    channels.set_xticks([0, 25, 50, 75, 100])
    channels.xaxis.grid(color="#dfdfdf", linewidth=.55)
    channels.tick_params(axis="y", length=0)
    channels.set_xlabel("Share of incident power (%)")
    labels_a = []
    for bar, value in zip(bars, percentages):
        label = f"{value:.{1 if value >= 10 else 2}f}%"
        labels_a.append(label)
        channels.annotate(label, (value, bar.get_y()+bar.get_height()/2),
                          xytext=(3.0, 0), textcoords="offset points",
                          va="center", ha="left", fontsize=7.7)
    channels.set_title("(a)  Outgoing power", loc="left", pad=8)

    x = [0, 1, 2]
    band_color = "#e0edf4"
    methods.axhspan(min(differences), max(differences), color=band_color,
                    linewidth=0, label="Range of three estimates")
    methods.axhline(0, color="#536776", linewidth=.7)
    points = methods.scatter(x, differences, s=27, color="#1b658f", zorder=3)
    methods.set_ylim(-6.35, 1.15)
    methods.set_xlim(-.38, 2.38)
    methods.set_yticks([-6, -4, -2, 0])
    methods.yaxis.grid(color="#e1e6e9", linewidth=.55)
    methods.set_xticks(x, ["Surface/PV", "Tetrahedron", "Generalized\nBloch"])
    methods.set_ylabel("Difference from surface/PV\n($10^{-5}$ of incident power)")
    methods.set_title("(b)  qT conversion by method", loc="left", pad=8)
    labels_b = []
    for i, value in enumerate(differences):
        label = "0" if value == 0 else f"{value:.2f}".replace("-", "−")
        labels_b.append(label)
        methods.annotate(label, (i, value), xytext=(0, 4.4), textcoords="offset points",
                         ha="center", va="bottom", fontsize=7.7)
    fig.legend(*methods.get_legend_handles_labels(), loc="center",
               bbox_to_anchor=(.802,.04), frameon=False, fontsize=7.6,
               handlelength=1.5, handleheight=.9)
    # Check coordinates and printed precision, including the left/right panel normalization.
    if labels_a != ["86.0%", "0.36%", "12.6%", "1.10%"] or labels_b != ["0", "−5.13", "−5.29"]:
        raise ValueError("Displayed roundings differ from the archived case")
    if [bar.get_width() for bar in bars] != percentages or list(points.get_offsets()[:, 1]) != differences:
        raise ValueError("Plotted values differ from the archived case")
    if abs((estimates["Tetrahedron"]/surface - 1) * 100 + data["surface_tetra_relative_percent"]) > 1e-13:
        raise ValueError("Tetrahedron relative difference is inconsistent")
    output.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime(2026, 11, 10, tzinfo=timezone.utc)
    fig.savefig(output, metadata={"Title": "Reference power partition and qT solver comparison",
                                  "CreationDate": stamp, "ModDate": stamp})
    plt.close(fig)
    return {"status": "PASS", "source_sha256": ARCHIVED_REFERENCE_SHA256,
            "panel_a_percentages_Tq_Rq_qT_qL": percentages,
            "panel_b_differences_in_units_of_1e-5": differences,
            "figure_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=root/"evidence/vi_results/REFERENCE_METHODS.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(make_figure(args.reference, args.output), indent=2))
