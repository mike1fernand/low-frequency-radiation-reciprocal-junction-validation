#!/usr/bin/env python3
"""Portable reconstruction of the generalized-Bloch reference scattering result.

This script uses only relative paths, the stored generalized-Bloch bulk matrices,
and the independently coded finite line-core evaluator. It deliberately does not
import the Phase-12 reduced-core scattering implementation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from independent_core_scattering import evaluate

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"

def load_g(name: str) -> np.ndarray:
    return np.load(RAW / name, allow_pickle=True)["G"]

G10 = load_g("gb_eta010_g256_mean.npz")
G5 = load_g("gb_eta005_g384_mean.npz")
G25 = load_g("gb_eta0025_g384_mean.npz")
xs = np.array([0.01, 0.005, 0.0025])
stack = np.stack([G10, G5, G25])
X = np.column_stack([np.ones(3), xs])
coef = np.linalg.lstsq(X, stack.reshape(3, -1), rcond=None)[0]
methods = {
    "low_pair_linear_richardson": 2.0 * G25 - G5,
    "three_point_quadratic": G10 / 3.0 - 2.0 * G5 + 8.0 * G25 / 3.0,
    "three_point_linear_least_squares": coef[0].reshape(G10.shape),
}
residue = np.load(RAW / "gb_real_axis_residue_mean.npz", allow_pickle=True)
g = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1])
out = {
    "method": "generalized-Bloch bulk matrices plus independently coded finite-core scattering",
    "omega": 1.2,
    "methods": {},
}
for name, G in methods.items():
    It = (G - G.conj().T) / (2j)
    IqL = residue["Im_residue_quasi_L"]
    ILp = residue["Im_residue_longitudinal_projector"]
    parts = {
        "qL": IqL,
        "qT": It - IqL,
        "Lproj": ILp,
        "Tproj": It - ILp,
    }
    r = evaluate(1.2, g, G, parts)
    out["methods"][name] = {
        "Rq": r["Rq"],
        "Tq": r["Tq"],
        "Tbulk": r["Tbulk"],
        "closure": r["closure"],
        "TL_modal": r["partials"]["qL"],
        "TT_modal": r["partials"]["qT"],
        "TT_fraction_modal": r["partials"]["qT"] / r["Tbulk"],
        "TL_projector": r["partials"]["Lproj"],
        "TT_projector": r["partials"]["Tproj"],
    }
keys = ["Rq", "Tq", "Tbulk", "TL_modal", "TT_modal", "TT_fraction_modal"]
out["method_spread"] = {}
for key in keys:
    vals = [v[key] for v in out["methods"].values()]
    out["method_spread"][key] = {
        "min": min(vals), "max": max(vals), "range": max(vals) - min(vals)
    }
path = RAW / "generalized_bloch_reference_independent_core_v23.json"
path.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
