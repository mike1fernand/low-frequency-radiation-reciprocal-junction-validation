"""Compare both Figure 6 PDFs with stored trajectory, shell and method records."""
import argparse
import csv
import json
from pathlib import Path
import re

import fitz
import numpy as np


def load(directory, name):
    return json.loads((directory / name).read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def axes(page, xlabels, ylabels):
    """Calibrate the vector coordinates from the displayed tick values and gridlines."""
    strokes = page.get_drawings()
    labels = page.get_text().splitlines()
    require(all(str(v) in labels for v in [*xlabels, *ylabels]), "Figure tick labels changed")
    grid = [d["rect"] for d in strokes if d["type"] == "s" and
            d["color"] and abs(d["color"][0] - 0.6902) < 0.001 and
            len(d["items"]) == 1 and d["items"][0][0] == "l"]
    xgrid = sorted(r.x0 for r in grid if r.width == 0)
    ygrid = sorted((r.y0 for r in grid if r.height == 0), reverse=True)
    require(len(xgrid) == len(xlabels) and len(ygrid) == len(ylabels), "Unexpected grid")
    xfit = np.polyfit(xgrid, xlabels, 1)
    yfit = np.polyfit(ygrid, ylabels, 1)
    return strokes, lambda xy: np.column_stack((
        np.polyval(xfit, np.asarray(xy)[:, 0]), np.polyval(yfit, np.asarray(xy)[:, 1])
    ))


def points(curve):
    segments = curve["items"]
    require(all(seg[0] == "l" for seg in segments), "Expected a vector polyline")
    return [segments[0][1], *(seg[2] for seg in segments)]


def window(strokes, convert, times, y_extent):
    matches = []
    for path in strokes:
        if path["type"] != "f" or not path["fill"] or not 0.4 < path["fill"][1] < 0.5:
            continue
        rect = path["rect"]
        bounds = convert([(rect.x0, rect.y1), (rect.x1, rect.y0)])
        if (np.allclose(bounds[:, 0], times, rtol=0, atol=2e-6) and
                np.allclose(bounds[:, 1], y_extent, rtol=0, atol=2e-6)):
            matches.append(rect)
    require(len(matches) == 1, "Shading must span the accepted interval and full y axis")
    return matches[0].y0


def curves(strokes, sizes):
    return [d for d in strokes if d["type"] == "s" and len(d["items"]) in sizes]


def curve_errors(convert, curve, x, y, ytol):
    data = convert(points(curve))
    np.testing.assert_allclose(data[:, 0], x, rtol=0, atol=2e-6)
    np.testing.assert_allclose(data[:, 1], y, rtol=0, atol=ytol)
    return float(max(abs(data[:, 1] - np.array(y))))


def verify(root, evidence):
    figures = root / "context/figures"
    modal = load(evidence, "REFERENCE_MODAL_ANALYSIS.json")["rows"]
    shell = load(evidence, "REFERENCE_SHELL_AUDIT.json")["rows"]
    figdata = load(evidence, "figure6a_data.json")
    with (evidence / "REFERENCE_TRAJECTORY_ROWS.csv").open() as src:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(src)]
    times = [m["time"] for m in modal]
    require(len(rows) == len(modal) == 46 and len(shell) == 14, "Trajectory row count")
    np.testing.assert_allclose(times, [r["time"] for r in rows], atol=1e-10, rtol=0)
    with fitz.open(figures / "paper1_sectionVII_fig6a_qT_window_v2.pdf") as doc:
        p = doc[0]
        ds, cv = axes(p, [48,50,52,54,56], [.122,.123,.124,.125,.126,.127,.128,.129])
        window(ds, cv, figdata["accepted_interval"], [.122, .129])
        trace, = curves(ds, [45])
        qt_error = curve_errors(cv, trace, times, [m["qT_production_over_E0"] for m in modal], 2e-9)
        horizontal = [d for d in ds if d["type"] == "s" and d["rect"].width > 195 and
                      d["rect"].height == 0 and d["color"][2] - d["color"][0] > .5]
        line, = horizontal
        central = cv(points(line))[0, 1]
        np.testing.assert_allclose(central, figdata["prediction_central"], atol=2e-9, rtol=0)
        # PyMuPDF returns the hatch tile in its own coordinate system, so read its
        # clipping rectangle directly from the PDF's page-coordinate vector commands.
        pattern = rb"/Pattern cs /H1 scn\s+([\d.]+) ([\d.]+) m\s+([\d.]+) ([\d.]+) l\s+([\d.]+) ([\d.]+) l"
        matched = re.search(pattern, p.read_contents())
        require(matched, "Missing two-method range in panel (a)")
        coords = [float(v) for v in matched.groups()]
        band = cv([(coords[0], p.rect.height-coords[1]), (coords[4], p.rect.height-coords[5])])
        np.testing.assert_allclose(sorted(band[:, 1]), figdata["prediction_range"], atol=2e-9, rtol=0)
        np.testing.assert_allclose(band[:, 0], [48,57], atol=2e-6, rtol=0)

    with fitz.open(figures / "paper1_sectionVII_fig6b_selector_controls_v2.pdf") as doc:
        p = doc[0]
        ds, cv = axes(p, [48,50,52,54,56], [0,1,2,3,4])
        plot_top = window(ds, cv, figdata["accepted_interval"], [0, 4.8])
        allcurves = curves(ds, [45,13])
        series = {
            "coupling": [r["T_coupling"] / .0025 for r in rows],
            "current": [abs(r["current"]) / .0005 for r in rows],
            "kinetic shell": [r["historical_shell"] / .002 for r in rows],
            "periodic/open": [abs(m["periodic_total_over_E0"] / r["T_bulk_open"] - 1) / .0025
                              for r, m in zip(rows, modal)],
            "drift": [abs(r["drift"]) / .00001 for r in rows],
            "total shell": [s["union"] / .002 for s in shell],
        }
        require(len(allcurves) == len(series) and all(k in p.get_text() for k in series),
                "Unexpected selector legend or missing plotted curve")
        errors = {}
        for (name, y), curve in zip(series.items(), allcurves):
            t = [s["time"] for s in shell] if name == "total shell" else times
            errors[name] = curve_errors(cv, curve, t, y, 2e-6)
        # All keys must fit on-page and clear the entire plotting area.
        text_keys = []
        for label in [*series, "threshold"]:
            matches = p.search_for(label)
            require(len(matches) == 1, "Missing or overlapping legend label " + label)
            bounds = matches[0]
            require(p.rect.contains(bounds) and bounds.y1 < plot_top - 2,
                    "Legend label obscures the data or exceeds the panel: " + label)
            text_keys.append(bounds.y1)
        guides = [d for d in ds if d["type"] == "s" and d["rect"].width > 195
                  and d["rect"].height == 0 and .14 < d["color"][0] < .16]
        guide, = guides
        np.testing.assert_allclose(cv(points(guide))[:,1], 1, rtol=0, atol=2e-6)
    return {"status": "PASS", "reference_taper": True, "direct_trace_points": len(times),
            "direct_trace_max_absolute_difference_in_energy_fraction": qt_error,
            "selector_and_shell_curves_checked": list(series),
            "selector_curve_points": 5*len(times), "post_selection_total_shell_points": len(shell),
            "largest_curve_differences_in_threshold_units": errors,
            "accepted_interval_on_both_panels": figdata["accepted_interval"],
            "panel_b_minimum_legend_clearance_points": plot_top - max(text_keys),
            "surface_line_and_two_method_band_match": True,
            "scope": "Vector artwork versus archived time-domain outputs and overlays; no additional trajectory simulation."}


if __name__ == "__main__":
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=here)
    ap.add_argument("--evidence", type=Path)
    args = ap.parse_args()
    print(json.dumps(verify(args.root, args.evidence or args.root / "evidence/packet_reference"), indent=2))
