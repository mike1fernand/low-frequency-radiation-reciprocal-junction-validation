"""Check the incident amplitude in Eq. (97) and unit-power rescaling in Eqs. (100)–(101).

Solve a transparent seven-site chain with its exact outgoing lead self-energies.
Measure the mean power on each interior bond over a full vibration period.
The only numerical dependency is NumPy.
"""
import json
import numpy as np


def bond_powers(displacement, omega, spring):
    phase = np.exp(-2j * np.pi * np.arange(64) / 64)
    wave = displacement[:, None] * phase
    velocity = (-1j * omega * wave).real
    force_from_right = spring * np.diff(wave.real, axis=0)
    return np.mean(-0.5 * (velocity[:-1] + velocity[1:]) * force_from_right, axis=1)


def main():
    spring, mass, sites = 1.69, 1.0, 7
    k_values = [0.13, 0.57, 1.31, 2.63]
    stiffness = spring * (2 * np.eye(sites) - np.eye(sites, k=1) - np.eye(sites, k=-1))
    left = np.eye(sites)[0]
    endpoints = np.diag([1, 0, 0, 0, 0, 0, 1])
    max_error = 0.0
    for k in k_values:
        omega = 2 * np.sqrt(spring / mass) * np.sin(k / 2)
        gamma = 2 * spring * np.sin(k)
        sigma = -spring * np.exp(1j * k)
        operator = mass * omega**2 * np.eye(sites) - stiffness - sigma * endpoints
        chi = np.linalg.solve(operator, np.sqrt(gamma) * left)
        np.testing.assert_allclose(np.abs(chi), 1 / np.sqrt(gamma), rtol=1e-12)
        np.testing.assert_allclose(chi[1:] / chi[:-1], np.exp(1j*k), rtol=1e-12)
        np.testing.assert_allclose(bond_powers(chi, omega, spring), omega / 4, rtol=1e-12)
        physical = (2 / np.sqrt(omega)) * chi
        error = np.max(np.abs(bond_powers(physical, omega, spring) - 1))
        max_error = max(max_error, error)
        assert error < 1e-12
    print(json.dumps({
        "status": "PASS",
        "model": "transparent seven-site scalar device between identical leads",
        "dimensionless_wavenumbers": k_values,
        "bonds_per_case": sites - 1,
        "incoming_Caroli_displacement_magnitude": "1/sqrt(gamma_q)",
        "Caroli_power": "omega/4",
        "amplitude_factor_for_unit_incident_power": "2/sqrt(omega)",
        "max_unit_power_absolute_error": float(max_error)
    }, indent=2))


if __name__ == "__main__":
    main()
