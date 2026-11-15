"""Exact algebra for Section VIII similarity and the prescribed-source long-wave orders.

Requires SymPy. Run from any directory. No numerical bulk Green matrix is used.
"""
import json
from itertools import product
import sympy as sp


def check_equal(actual, expected):
    if sp.simplify(actual - expected) != 0:
        raise ValueError(f"Identity failed: {actual} versus {expected}")


def similarity():
    # These are arbitrary entries of the three operators. The statements apply
    # entry by entry to the entire graph, including its leads and contact.
    mass, bare, contact, omega = sp.symbols("mass bare contact omega", real=True)
    stiffness_scale, mass_scale = sp.symbols("alpha_K alpha_M", positive=True)
    frequency_factor = sp.sqrt(stiffness_scale / mass_scale)
    inertial = (omega * frequency_factor)**2 * mass_scale * mass
    unprimed = bare + contact - omega**2 * mass
    full = stiffness_scale * (bare + contact) - inertial
    check_equal(full, stiffness_scale * unprimed)
    fixed_contact = stiffness_scale * bare + contact - inertial
    check_equal(fixed_contact, stiffness_scale * (bare + contact / stiffness_scale - omega**2 * mass))

    # At equal displacement amplitude the force of each scaled spring changes
    # by alpha_K. Compute the velocity directly from the primed trajectory.
    tau = sp.symbols("tau", real=True)
    trajectory = sp.Function("x")
    velocity = sp.diff(trajectory(frequency_factor * tau), tau)
    expected_velocity = sp.diff(trajectory(tau), tau).subs(tau, frequency_factor*tau)
    check_equal(velocity, frequency_factor * expected_velocity)
    power_factor = stiffness_scale * frequency_factor
    check_equal(power_factor**2, stiffness_scale**3 / mass_scale)
    return {
        "complete_rescaling_operator_identity": "K' - (Omega omega)^2 M' = alpha_K (K - omega^2 M)",
        "fixed_contact_operator_identity": "K' - (Omega omega)^2 M' = alpha_K (K_bare + K_c/alpha_K - omega^2 M)",
        "equal_displacement_power_factor": "alpha_K sqrt(alpha_K/alpha_M)",
    }


def corner_projections():
    # t = rho*a. For a centered segment take local Bloch endpoint phases
    # exp(+/- i b t/2) and Fourier phase exp[-i (s_y y + s_z z)t/2].
    # g*Q is factored out. k = rho*(x,y,z), kappa = rho*b.
    x, y, z, b = sp.symbols("x y z b", real=True)
    sv = (y**2 + z**2, -x*y, -x*z)
    sh = (0, -z, y)
    # Normalizing these orthogonal-to-k vectors divides by a t-independent
    # nonzero constant for a fixed non-axial direction, preserving the orders.
    for direction in (sv, sh):
        check_equal(sum(a*c for a, c in zip((x, y, z), direction)), 0)
    check_equal(sum(a*c for a, c in zip(sv, sh)), 0)

    def projected_coefficient(sector, polarization, degree):
        coefficient = 0
        for sy, sz, endpoint in product((-1, 1), repeat=3):
            # Axial source: (2/3)*mean(q) along x at each corner.
            # Radial source: -(1/3)*(q_right-q_left)*(0,s_y,s_z).
            force = polarization[0] if sector == "X" else -endpoint*(sy*polarization[1] + sz*polarization[2])
            phase = sp.I * (endpoint*b - sy*y - sz*z) / 2
            coefficient += force * phase**degree
        return sp.expand(coefficient / (3*sp.factorial(degree)))

    cases = [
        ("X", "SV", sv, 0, 8*(y**2+z**2)/3),
        ("R", "SV", sv, 2, 2*b*x*(y**2+z**2)/3),
        ("R", "SH", sh, 4, b*y*z*(y**2-z**2)/18),
    ]
    answers = {}
    for sector, label, polarization, order, expected in cases:
        for power in range(order):
            check_equal(projected_coefficient(sector, polarization, power), 0)
        check_equal(projected_coefficient(sector, polarization, order), expected)
        answers[sector+"_"+label] = {
            "order": order,
            "coefficient_for_unnormalized_polarization": str(sp.factor(expected)),
        }
    return answers


if __name__ == "__main__":
    print(json.dumps({
        "status": "PASS",
        "similarity": similarity(),
        "corner_projection_orders": corner_projections(),
        "scope": "Exact operator and four-corner Taylor identities. A finite motif sum adds analytic phase factors. The source orders are upper orders and can vanish in special directions. They concern prescribed line motion; radiated power needs the coupled scattering solution.",
    }, indent=2))
