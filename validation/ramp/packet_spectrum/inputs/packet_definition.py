#!/usr/bin/env python3
"""Frozen physical packet preparation and bound direct comparison authority.

Preparation and comparison algebra are copied verbatim from the prior verified
producer. This module removes unrelated old frequency-curve dependencies; it
neither changes q,v nor imports old GB values into the new comparison.
"""
from pathlib import Path
import json,hashlib,math
import numpy as np
HERE=Path(__file__).resolve().parent

def require(test,message):
    if not test:raise ValueError(message)

def load_inputs(directory=HERE/'inputs'):
    directory=Path(directory);manifest=json.loads((directory/'MANIFEST.json').read_text());items={}
    for name,expected in manifest['files'].items():
        require(Path(name).name==name,'Input must be basename')
        raw=(directory/name).read_bytes();require(hashlib.sha256(raw).hexdigest()==expected,'Input hash mismatch: '+name);items[name]=json.loads(raw)
    modal=items['direct_modal.json'];selector=items['selector.json']
    require(modal['schema']=='paper1.p07.invariant-modal-analysis.v1','Wrong direct modal schema')
    require(np.allclose(modal['interval'],[52.8,55.4],rtol=0,atol=1e-12),'Direct window mismatch')
    require(np.allclose(selector['selected_interval'],[52.8,55.4],rtol=0,atol=1e-12) and selector['contract']['uses_frequency_domain_target'] is False,'Selector is not the fixed blind interval')
    lo=modal['summary']['qT_lower_time_average'];hi=modal['summary']['qT_upper_time_average']
    require(math.isfinite(lo) and math.isfinite(hi) and 0<=lo<=hi<=1,'Invalid direct allocation')
    declared=items['direct_comparison.json']['direct_midpoint']
    require(abs(float(declared)-.5*(lo+hi))<=5e-16,'Declared direct scalar differs materially from actual modal bounds')
    items['retained_comparison.json']={'direct_midpoint':declared}
    spectrum,*_=packet_spectrum(items['preparation.json'])
    require(math.isclose(modal['E0'],spectrum['open_line_energy'],rel_tol=1e-12),'Direct energy and packet preparation differ')
    return items,manifest

def packet_spectrum(config):
    """Current real packet and full-grid Parseval spectrum, including Nyquist."""
    n = config["n"]
    require(isinstance(n, int) and n >= 4 and n % 2 == 0, "An even grid with n >= 4 is required")
    for name in ("a", "m_q", "K_q", "sigma", "omega0"):
        require(math.isfinite(config[name]) and config[name] > 0, f"Invalid {name}")
    a, mass, stiffness = config["a"], config["m_q"], config["K_q"]
    require(math.isfinite(config["amplitude"]) and config["amplitude"] != 0, "Invalid amplitude")
    require(math.isfinite(config["center"]), "Invalid packet center")
    maximum_frequency = 2 * math.sqrt(stiffness / mass)
    require(config["omega0"] < maximum_frequency, "Carrier must be inside the line band")
    require(config["nyquist_velocity"] == "zero", "Unsupported Nyquist preparation")
    x = a * np.arange(n, dtype=float)
    k = 2 * np.pi * np.fft.fftfreq(n, d=a)
    omega = maximum_frequency * np.abs(np.sin(k * a / 2))
    carrier = 2 * math.asin(config["omega0"] / maximum_frequency) / a
    displacement = config["amplitude"] * np.exp(-0.5 * ((x - config["center"]) / config["sigma"]) ** 2) * np.cos(carrier * (x - config["center"]))
    direction = np.sign(k)
    direction[n // 2] = 0.0
    velocity = np.fft.ifft(-1j * direction * omega * np.fft.fft(displacement)).real
    qhat, vhat = np.fft.fft(displacement), np.fft.fft(velocity)
    weight = mass * (np.abs(vhat) ** 2 + omega ** 2 * np.abs(qhat) ** 2) / (2 * n)
    low, high = config["band"]
    require(0 <= low < high <= maximum_frequency, "Invalid comparison band")
    inside = (omega >= low) & (omega <= high)
    total = float(weight.sum())
    require(total > 0 and np.all(np.isfinite(weight)), "Invalid packet energy")
    require(float(weight[inside].sum()) > 0, "Comparison band has no packet energy")
    periodic = 0.5 * mass * float(np.dot(velocity, velocity)) + 0.5 * stiffness * float(np.sum((np.roll(displacement, -1) - displacement) ** 2))
    open_energy = 0.5 * mass * float(np.dot(velocity, velocity)) + 0.5 * stiffness * float(np.sum(np.diff(displacement) ** 2))
    coverage = float(weight[inside].sum() / total)
    summary = {
        "n": n, "carrier_kappa_a": carrier * a,
        "periodic_parseval_energy": total, "periodic_real_space_energy": periodic,
        "open_line_energy": open_energy,
        "parseval_relative_residual": abs(periodic - total) / total,
        "periodic_vs_open_relative_difference": abs(total - open_energy) / open_energy,
        "nyquist_energy_fraction": float(weight[n // 2] / total),
        "coverage": coverage, "outside_fraction": 1 - coverage,
        "comparison_band": [low, high],
        "weight_definition": "m_q*(abs(vhat)^2+omega_q^2*abs(qhat)^2)/(2*N), full DFT grid",
    }
    return summary, omega, weight, inside, displacement, velocity


def compare_prediction(prediction, coverage, direct):
    """Compare S=sum(in-band W*T)/sum(all W), not the conditional mean S/c."""
    require(math.isfinite(coverage) and 0 < coverage <= 1, "Invalid spectral coverage")
    require(math.isfinite(prediction) and 0 < prediction <= coverage, "In-band contribution must lie in (0, coverage]")
    require(math.isfinite(direct) and 0 <= direct <= 1, "Invalid direct-dynamics fraction")
    lower = prediction
    upper = lower + 1 - coverage
    return {
        "conditional_in_band_prediction": float(prediction / coverage),
        "in_band_contribution_normalized_by_full_energy": float(lower),
        "prediction_normalized_discrepancy_percent": float(100 * abs(prediction - direct) / prediction),
        "full_packet_passivity_interval": [float(lower), float(upper)],
        "endpoint_normalized_discrepancies_percent": [float(100 * abs(lower - direct) / lower), float(100 * abs(upper - direct) / upper)],
    }

