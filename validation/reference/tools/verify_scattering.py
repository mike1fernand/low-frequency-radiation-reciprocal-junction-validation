"""Check the paper's harmonic contact and scattering equations at omega = 1.2.

Usage: python3 verify_scattering.py [reference-kernel.npz]
The selected numerical kernel is a saved input; this script does not re-integrate it.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

OMEGA = 1.2
K_LINE = 1.69
STRENGTHS = np.array([.1, .2, .2, .2, .2, .1])
PHASE = 2 * np.arcsin(OMEGA / (2 * np.sqrt(K_LINE)))
GAMMA = 2 * K_LINE * np.sin(PHASE)
BONDS = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0],
                  [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]])
DIRECTIONS = BONDS / np.linalg.norm(BONDS, axis=1)[:, None]


def assert_close(actual, expected, tolerance=2e-12):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=tolerance)


def hessian(profile):
    """Compute in reduced coordinates and independently from each of 48 springs."""
    size = len(profile)
    U = np.zeros((2 * size, 2 * size))
    V = np.zeros((2 * size, size + 1))
    K = np.zeros((size + 1, size + 1))
    microscopic = np.zeros((3 * size + 1, 3 * size + 1))
    for i, g in enumerate(profile):
        z = slice(2*i, 2*i + 2)
        q = slice(i, i + 2)
        U[z, z] = 2*g/3 * np.diag([1, 2])
        V[z, q] = 2*g/3 * np.array([[-1, -1], [-np.sqrt(2), np.sqrt(2)]])
        K[q, q] += 4*g/3 * np.eye(2)
        for sy, sz in itertools.product((-1, 1), repeat=2):
            patterns = np.array([[.5, 0, 0], [0, sy, sz]])  # axial and radial
            patterns[1] /= 2*np.sqrt(2)
            for side, dx in enumerate((1, -1)):
                bond = np.array([dx, sy, sz])/np.sqrt(3)
                row = np.zeros(3 * size + 1)
                row[z] = patterns @ bond
                row[2*size+i+side] = -bond[0]
                microscopic += g * np.outer(row, row)
    assert_close(microscopic, np.block([[U, V], [V.T, K]]))
    return U, V, K


def modal_work(force, shell):
    return OMEGA/2 * np.vdot(force, shell @ force).real


def scatter(profile, g0, shells, incident=0):
    U, V, K = hessian(profile)
    nc, nq = V.shape
    stiffness = K_LINE*(2*np.eye(nq)-np.eye(nq, k=1)-np.eye(nq, k=-1))
    end = np.zeros(nq)
    end[incident] = 1
    ends = np.diag([1] + [0]*(nq-2) + [1])
    Z = OMEGA**2*np.eye(nq) - stiffness - K + K_LINE*np.exp(1j*PHASE)*ends
    core = np.linalg.solve(np.eye(nc) + g0@U, g0)
    feedback = -V.T@core@V
    device = np.linalg.inv(Z - feedback)
    chi = device @ (np.sqrt(GAMMA)*end)
    drive = np.linalg.solve(np.eye(nc) + U@g0, V)
    eta = drive @ chi
    coupled = np.block([[Z, -V.T], [g0@V, np.eye(nc)+g0@U]])
    full = np.linalg.solve(coupled, np.r_[np.sqrt(GAMMA)*end, np.zeros(nc)])
    assert_close(full[:nq], chi)
    assert_close(full[nq:], -g0@eta)
    assert_close(-V@full[:nq]-U@full[nq:], -eta)
    widths = [2 * drive.conj().T @ shell @ drive for shell in shells]
    total_width = 1j*(feedback-feedback.conj().T)
    assert_close(sum(widths), total_width)
    assert_close(device-device.conj().T,
                 -1j * device @ (GAMMA*ends + total_width) @ device.conj().T)
    fractions = []
    phys_force = -2*eta/np.sqrt(OMEGA)
    for width, shell in zip(widths, shells):
        assert np.min(np.linalg.eigvalsh(width)) > -1e-12
        from_work = modal_work(phys_force, shell)
        assert_close(np.vdot(chi, width@chi), from_work)
        fractions.append(from_work)
    assert_close(sum(fractions), .5*np.vdot(phys_force, -1j*OMEGA*g0@phys_force).real)
    reflected = abs(1-1j*GAMMA*device[incident, incident])**2
    other = nq-1 if incident == 0 else 0
    forward = GAMMA**2*abs(device[other,incident])**2
    powers = np.array([reflected, forward, *fractions])
    assert_close(powers.sum(), 1.)
    assert np.min(powers) >= -1e-13
    return powers, chi, eta


def check_bloch_current():
    """Use the directed bond force and an independent finite-difference speed."""
    k = np.array([.73, -.48, 1.13])
    def eig(kvec):
        D = 2 * DIRECTIONS.T @ (DIRECTIONS * (1 - np.cos(BONDS @ kvec))[:,None])
        return np.linalg.eigh(D)
    values, vectors = eig(k)
    step = 2e-5
    gradient = np.array([(np.sqrt(eig(k+step*e)[0])-np.sqrt(eig(k-step*e)[0]))
                         /(2*step) for e in np.eye(3)])
    max_error = 0.
    for j, w in enumerate(np.sqrt(values)):
        e = vectors[:,j]
        phase = np.exp(1j*BONDS@k)
        relative = (phase[:,None]-1)*e
        forces = DIRECTIONS*np.sum(DIRECTIONS*relative,axis=1)[:,None]
        summed_velocity = -1j*w*(1+phase[:,None])*e
        directed = -.25*np.sum(summed_velocity.conj()*forces,axis=1).real
        transport = BONDS.T@directed
        target = .5*w**2 * gradient[:,j]
        error = np.max(np.abs(transport-target))
        max_error = max(max_error,error)
        assert_close(transport,target,3e-9)
        assert np.linalg.norm(2*transport-target) > .1
    return max_error


def check_weak_witness(chi0):
    """Locate the SM3 witness shell using bisection from its rounded coordinates."""
    k = np.array([.6318952450, -.6028746482, -.9641061600])
    def eig():
        D = 2 * DIRECTIONS.T @ (DIRECTIONS*(1-np.cos(BONDS@k))[:,None])
        return np.linalg.eigh(D)
    l, r = k[2]-.0001, k[2]+.0001
    k[2] = l
    left_sign = np.sign(eig()[0][0]-OMEGA**2)
    k[2] = r
    assert left_sign * (eig()[0][0]-OMEGA**2) < 0
    for _ in range(50):
        k[2] = (l+r)/2
        if left_sign * (eig()[0][0]-OMEGA**2) > 0:
            l = k[2]
        else:
            r = k[2]
    eigenvalues, eigenvectors = eig()
    e = eigenvectors[:,0]
    assert_close(eigenvalues[0], OMEGA**2)
    gradients = 2*BONDS.T @ (np.sin(BONDS@k) * (DIRECTIONS@e)**2)
    phi = np.zeros(12,dtype=complex)
    for i in range(6):
        for y,z in itertools.product((-1, 1), repeat=2):
            phase = np.exp(1j*k @ np.array([i+.5,y/2,z/2]))
            phi[2*i] += e[0]/2*phase
            phi[2*i+1] += (y*e[1]+z*e[2])/(2*np.sqrt(2))*phase
    eta1 = hessian(STRENGTHS)[1] @ chi0
    overlap = abs(np.vdot(phi,eta1))**2
    metrics = [np.linalg.norm(gradients),*np.diff(eigenvalues),overlap]
    assert_close(metrics,[1.8578005,.1007240,2.6995796,.7492079],8e-7)
    return metrics


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel',type=Path,nargs='?')
    args = parser.parse_args()
    local = root/'evidence'/'REFERENCE_KERNEL.npz'
    kernel = args.kernel or (local if local.exists() else root.parent/'evidence'/'REFERENCE_KERNEL.npz')
    sha = hashlib.sha256(kernel.read_bytes()).hexdigest()
    assert sha == '82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
    with np.load(kernel, allow_pickle=False) as data:
        g0 = data['G0'][:12,:12]
        shells = [data['Im_qT'][:12,:12], data['Im_qL'][:12,:12]]
        assert_close(shells[0], (data['Im_branch_0']+data['Im_branch_1'])[:12,:12])
        assert_close(shells[1], data['Im_branch_2'][:12,:12])
    assert_close(g0,g0.T)
    assert_close(sum(shells),(g0-g0.conj().T)/(2j))
    for shell in shells:
        assert np.min(np.linalg.eigvalsh(shell)) >= -1e-12
    ref, _, _ = scatter(STRENGTHS, g0, shells)
    assert_close(ref, [.00357214927124,.859745958528,.125659255392,.0110226367963],5e-10)
    blank, chi0, _ = scatter(0*STRENGTHS, g0, shells)
    assert_close(blank,[0,1,0,0])
    assert_close(np.abs(2*chi0/np.sqrt(OMEGA)), 2/np.sqrt(OMEGA*GAMMA))
    # Asymmetric contact exercises directional bulk loss in a reciprocal system.
    ramp = np.arange(1.,7.)/7
    forward = scatter(ramp,g0,shells)[0]
    reverse = scatter(ramp,g0,shells,incident=-1)[0]
    assert_close(forward[1],reverse[1])
    assert abs(forward[2:].sum()-reverse[2:].sum()) > .01
    eta1 = hessian(STRENGTHS)[1] @ chi0
    weak = 2*np.array([np.vdot(eta1, s@eta1).real for s in shells])
    assert_close(weak,[.0969632610,.0076249508],8e-8)
    errors = []
    for scale in [.004,.002,.001]:
        powers = scatter(scale*STRENGTHS,g0,shells)[0]
        errors.append(float(np.linalg.norm(powers[2:]/scale**2 - weak)))
    assert all(.47 < b/a < .53 for a,b in zip(errors,errors[1:]))
    witness = check_weak_witness(chi0)
    current_error = check_bloch_current()
    print(json.dumps({'status':'PASS','kernel_sha256':sha,
                      'reference_Rq_Tq_TqT_TqL':ref.tolist(),
                      'weak_qT_qL':weak.tolist(),
                      'SM3_gradient_lower_gap_upper_gap_squared_overlap':witness,
                      'qT_qL_weak_coefficient_errors_at_scales_0.004_0.002_0.001':errors,
                      'cell_current_finite_difference_max_absolute_error':current_error,
                      'note':'Numerical identity checks using retained quadrature. The contact Hessian, shell witness, and bond current are independently constructed.'},indent=2))


if __name__ == '__main__':
    main()
