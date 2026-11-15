"""Compare the four vector curves in Fig. 5 with the archived frequency scan."""
from pathlib import Path
import json

import fitz
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CURVES = {
    "Rq": (0.185, 0.185, 0.185),
    "Tq": (0.0, 0.447, 0.698),
    "T_qT_modal": (0.835, 0.369, 0.0),
    "T_qL_modal": (0.0, 0.620, 0.451),
}


def check():
    scan = json.loads((ROOT / "evidence/vi_results/ACTIVE85.json").read_text())["families"]["frequency"]
    freq = np.array([r["omega"] for r in scan])
    figure = ROOT / "context/figures/paper1_frequency_partition.pdf"
    with fitz.open(figure) as pdf:
        page = pdf[0]
        paths = page.get_drawings()
        labels = page.get_text("words")
    # Axis tick positions come from the long, regularly spaced gray gridlines.
    # The tick labels printed on each axis independently set its physical range.
    grid = [d["rect"] for d in paths if d["color"] and
            np.max(np.abs(np.asarray(d["color"]) - 0.69)) < 0.01 and len(d["items"]) == 1]
    xgrid = sorted(g.x0 for g in grid if g.height > 100 and g.width < 0.01)
    ygrid = sorted((g.y0 for g in grid if g.width > 100 and g.height < 0.01), reverse=True)
    xticks = [float(w[4]) for w in labels if w[1] > max(ygrid) and w[4] in {f"{x/10:.1f}" for x in range(6, 19, 2)}]
    yticks = [float(w[4]) for w in labels if w[2] < min(xgrid) and w[4] in {f"{y/10:.1f}" for y in range(0, 11, 2)}]
    xvalues = sorted(xticks)
    yvalues = sorted(yticks)
    assert len(xgrid) == len(xvalues) == 7 and len(ygrid) == len(yvalues) == 6
    xscale = np.polyfit(xgrid, xvalues, 1)
    yscale = np.polyfit(ygrid, yvalues, 1)
    assert np.allclose(np.polyval(xscale, xgrid), xvalues, atol=1e-6)
    assert np.allclose(np.polyval(yscale, ygrid), yvalues, atol=1e-6)
    errors = {}
    for field, rgb in CURVES.items():
        candidates = [d["items"] for d in paths if len(d["items"]) == len(scan)-1 and
                      d["color"] and np.allclose(d["color"], rgb, atol=0.002)]
        assert len(candidates) == 1, field
        segments = candidates[0]
        assert all(op == "l" for op, _, _ in segments), field
        assert all(a[2] == b[1] for a, b in zip(segments, segments[1:])), field
        points = [segments[0][1]] + [line[2] for line in segments]
        assert np.max(np.abs(np.polyval(xscale, [p.x for p in points]) - freq)) < 1e-6
        chart = np.polyval(yscale, [p.y for p in points])
        expected = np.array([row[field] for row in scan])
        err = np.max(abs(chart - expected))
        assert err < 1e-6, (field, err)
        errors[field] = float(err)
    return {"status": "PASS", "samples_per_curve": len(scan), "curves": len(CURVES),
            "max_chart_coordinate_difference_in_power_fraction": errors,
            "scope": "Agreement of plotted vector vertices with archived samples, at PDF coordinate precision."}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
