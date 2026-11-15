"""Stationary mass-ratio control of Section VIII at fixed springs and bulk.

The 19-coordinate solve uses physical displacements for a unit-amplitude
incident line wave. The 7-coordinate check eliminates the corners in
mass-normalized line variables with a unit Caroli input. Both use the exact
propagating-lead boundary condition and the indicated archived bulk quadrature.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from verify_scattering import hessian

OMEGA = 1.2
K = 1.69
TAPER = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1])
MASSES = (0.75, 1.0, 1.25)
SURFACE_SHA = "82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b"
TETRA_SHA = "609c5990149bfe8f5975571897d3b069c94661b7d23955e2fd513e4dea18c4ef"
LINE = 2*np.eye(7) - np.eye(7, k=1) - np.eye(7, k=-1)
ENDS = np.diag([1, 0, 0, 0, 0, 0, 1])


def close(a, b, atol=2e-12):
    np.testing.assert_allclose(a, b, rtol=0, atol=atol)


def physical_response(mass, profile, g0, shells):
    """Solve physical corner and line displacements including both identical leads."""
    phase = 2*np.arcsin(OMEGA*np.sqrt(mass/K)/2)
    if not 0 < phase < np.pi:
        raise ValueError("The incident line wave must propagate")
    gamma = 2*K*np.sin(phase)
    # hessian also checks this construction against the 48 individual bonds.
    u, v, k_contact = hessian(profile)
    line = OMEGA**2*mass*np.eye(7) - K*LINE - k_contact + K*np.exp(1j*phase)*ENDS
    joint = np.block([[line, -v.T], [g0@v, np.eye(12) + g0@u]])
    incoming = np.zeros(19, dtype=complex)
    incoming[0] = 1j*gamma
    motion = np.linalg.solve(joint, incoming)
    q, corner = motion[:7], motion[7:]
    force = -v@q-u@corner
    close(corner, g0@force)
    incident_power = OMEGA*gamma/4
    bulk = [OMEGA*np.vdot(force,shell@force).real/(2*incident_power) for shell in shells]
    return np.array([abs(q[0]-1)**2, abs(q[-1])**2, *bulk])


def mass_weighted_schur(mass, profile, g0, shells):
    """Independent 7-coordinate Schur solve using the squared paired extensions."""
    h = np.zeros((19, 19))
    for i, strength in enumerate(profile):
        q1, q2 = i, i+1
        x, radial = 7+2*i, 8+2*i
        ax = np.zeros(19)
        ax[[q1,q2,x]] = [-1/np.sqrt(mass),-1/np.sqrt(mass),1]
        rad = np.zeros(19)
        rad[[q1,q2,radial]] = [-1/np.sqrt(mass),1/np.sqrt(mass),np.sqrt(2)]
        h += (2*strength/3)*(np.outer(ax,ax)+np.outer(rad,rad))
    u, v = h[7:,7:], h[7:,:7]
    k_mass = K/mass
    phase = np.arccos(1-OMEGA**2/(2*k_mass))
    lead = -k_mass*np.exp(1j*phase)
    gamma = -2*lead.imag
    dressed = np.linalg.solve(np.eye(12)+u@g0, v)
    response = np.linalg.solve(OMEGA**2*np.eye(7) - k_mass*LINE - h[:7,:7]
                               - lead*ENDS + v.T@g0@dressed, np.eye(7)[:,0])
    field = dressed@response
    return np.array([abs(1-1j*gamma*response[0])**2, gamma**2*abs(response[-1])**2,
                     *(2*gamma*np.vdot(field,shell@field).real for shell in shells)])


def run(method, kernel, expected_sha, ref):
    if sha256(kernel.read_bytes()).hexdigest() != expected_sha:
        raise ValueError(f"Incorrect {method} archived kernel")
    with np.load(kernel,allow_pickle=False) as d:
        g0 = d["G0"][:12,:12]
        if method == "surface_N192_M256":
            shells = [d["Im_qT"][:12,:12], d["Im_qL"][:12,:12]]
        else:
            shells = [d["Im0"]+d["Im1"], d["Im2"]]
    close(sum(shells),(g0-g0.conj().T)/(2j))
    close(g0,g0.T)
    for shell in shells:
        close(shell,shell.conj().T)
        if np.linalg.eigvalsh(shell).min() < -1e-12:
            raise ValueError("Negative modal spectral weight")
    rows = []
    for mass in MASSES:
        powers = physical_response(mass,TAPER,g0,shells)
        control = mass_weighted_schur(mass,TAPER,g0,shells)
        close(powers,control)
        close(sum(powers),1)
        close(physical_response(mass,0*TAPER,g0,shells),[0,1,0,0])
        if mass == 1:
            # Archived nominal powers used a small line regulator.
            close(powers,ref,3e-10)
        if min(powers)<0 or max(powers)>1:
            raise ValueError("Passivity failed")
        rows.append(dict(mq_over_mb=mass,Rq=powers[0],Tq=powers[1],TqT=powers[2],TqL=powers[3],
                         closure=float(abs(sum(powers)-1)),
                         formulation_max_abs_difference=float(np.max(abs(powers-control)))))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    evidence = args.root/"evidence"
    a=json.loads((evidence/"REFERENCE_KERNEL.json").read_text())["scattering"]
    b=json.loads((evidence/"controls/DOWNSTREAM_TETRA_N192.json").read_text())
    if b["input_sha256"] != TETRA_SHA:
        raise ValueError("Independent tetrahedron reference digest differs")
    source=run("surface_N192_M256",evidence/"REFERENCE_KERNEL.npz",SURFACE_SHA,
               [a[k] for k in ("Rq","Tline","qT","qL")])
    tetra=run("tetra_N192",evidence/"line_mass/TETRA_N192.npz",TETRA_SHA,
              [b["reference"][k] for k in ("Rq","Tq","TT_modal","TL_modal")])
    result = {
        "status":"PASS","omega":OMEGA,"Kq":K,"mb":1,"a":1,"K_axis":1,"K_diag":1,
        "g":TAPER.tolist(),"surface_N192_M256":source,"tetra_N192":tetra,
        "max_qT_difference_percentage_points":
            100*max(abs(x["TqT"]-y["TqT"]) for x,y in zip(source,tetra)),
        "scope":"Three finite one-factor examples, selected before evaluation. All seven scalar contact sites and both scalar leads have the selected mass. New contact solves with archived numerical bare-bulk matrices. The inter-method difference is not an error bound."
    }
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
