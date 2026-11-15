#!/usr/bin/env python3
"""Check a declared 3D vector-parent realization of the P1 source motif.

This finite 6-node diagnostic is not the infinite P1 radiation problem and
does not derive the graph or its coefficients from homogeneous DST.
"""
from pathlib import Path
import hashlib, json, sys, tempfile
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
INPUTS=HERE.parent/'recheck/purity_inputs'

def main():
    provenance=json.loads((INPUTS/'INPUT_PROVENANCE.json').read_text())
    known={r['path']:r for r in provenance['files']}
    modules={n:(INPUTS/n).read_bytes() for n in ['model_api_v1_1.py','reduced_core.py']}
    for n,data in modules.items():assert hashlib.sha256(data).hexdigest()==known[n]['sha256']
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-vector-parent-') as td:
        for n,data in modules.items():(Path(td)/n).write_bytes(data)
        sys.path.insert(0,td)
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices
        # Nonunit allowed masses make an accidental unit-mass reduction fail.
        par=ModelParameters(m_spine=1.7,m_bulk=.8);g=.2
        signs=[(-1,-1),(-1,1),(1,-1),(1,1)]
        Kc=np.zeros((18,18))
        for j,(sy,sz) in enumerate(signs):
            for endpoint,n in [(0,np.array([1.,sy,sz])/np.sqrt(3)),
                               (1,np.array([-1.,sy,sz])/np.sqrt(3))]:
                b=np.zeros(18);b[3*endpoint:3*endpoint+3]=-n
                b[3*(j+2):3*(j+2)+3]=n
                Kc+=g*np.outer(b,b)
        spine=np.zeros(18);spine[0]=-1;spine[3]=1
        Kq=par.k_spine*np.outer(spine,spine)
        K=Kc+Kq
        B=np.zeros((18,4));B[0,0]=1;B[3,1]=1
        for j,(sy,sz) in enumerate(signs):
            B[3*(j+2),2]=.5
            B[3*(j+2)+1,3]=sy/(2*np.sqrt(2))
            B[3*(j+2)+2,3]=sz/(2*np.sqrt(2))
        assert np.max(abs(B.T@B-np.eye(4)))<1e-15
        group=[]
        for swap in [False,True]:
            perm=np.array([[0.,1.],[1.,0.]]) if swap else np.eye(2)
            for sy in [-1.,1.]:
                for sz in [-1.,1.]:
                    transverse=np.diag([sy,sz])@perm
                    R=np.zeros((3,3));R[0,0]=1;R[1:,1:]=transverse
                    T=np.zeros((18,18));T[:3,:3]=R;T[3:6,3:6]=R
                    for j,sign in enumerate(signs):
                        target=tuple(int(x) for x in transverse@np.array(sign))
                        dest=signs.index(target)
                        T[3*(dest+2):3*(dest+2)+3,3*(j+2):3*(j+2)+3]=R
                    group.append(T)
        P=sum(group)/len(group)
        covariance=max(float(np.max(abs(T.T@K@T-K))) for T in group)
        projector=float(np.max(abs(P-B@B.T)))
        leakage=float(np.max(abs((np.eye(18)-P)@K@P)))
        U,V,Kself=reduced_static_matrices([g],par)
        expected=np.block([[Kself,V.conj().T],[V,U]])
        reduced_error=float(np.max(abs(B.T@Kc@B-expected)))
        M=np.diag([par.m_spine]*6+[par.m_bulk]*12)
        Minv=np.linalg.inv(M)
        Mred=B.T@M@B
        symplectic=np.block([[B,np.zeros_like(B)],[np.zeros_like(B),B]])
        J18=np.block([[np.zeros((18,18)),np.eye(18)],[-np.eye(18),np.zeros((18,18))]])
        J4=np.block([[np.zeros((4,4)),np.eye(4)],[-np.eye(4),np.zeros((4,4))]])
        canonical=float(np.max(abs(symplectic.T@J18@symplectic-J4)))
        Afull=np.block([[np.zeros((18,18)),Minv],[-K,np.zeros((18,18))]])
        Ared=np.block([[np.zeros((4,4)),np.linalg.inv(Mred)],[-B.T@K@B,np.zeros((4,4))]])
        z0=np.array([.17,-.09,.03,-.02,.07,.01,-.04,.02])
        zfull0=symplectic@z0
        dynamics=max(float(np.max(abs(expm(t*Afull)@zfull0-symplectic@(expm(t*Ared)@z0))))
                     for t in [.3,1.1,2.7])
        # An arbitrary non-invariant corner displacement produces a transverse
        # line force. This rejects the stronger unrestricted-scalar claim.
        arbitrary=np.zeros(18);arbitrary[3*2+1]=1
        force=-Kc@arbitrary
        transverse_force=float(np.linalg.norm(force[[1,2,4,5]]))
        assert transverse_force>.01
        assert max(covariance,projector,leakage,reduced_error,canonical,dynamics)<1e-12
    out={
        'status':'PASS','scope':'Declared microscopic mechanical parent, one actual P1 motif; static and finite-dimensional harmonic checks only',
        'nodes':6,'full_displacement_dofs':18,'invariant_displacement_dofs':4,
        'group_elements':len(group),'contact_stiffness':g,
        'mass_fixture':{'line':par.m_spine,'bulk':par.m_bulk},
        'max_hessian_covariance_error':covariance,'fixed_projector_error':projector,
        'excluded_sector_hessian_leakage':leakage,'actual_P1_reduced_contact_block_error':reduced_error,
        'canonical_pullback_error':canonical,'finite_harmonic_flow_error':dynamics,
        'arbitrary_bulk_transverse_line_force_norm':transverse_force,
        'limitations':['No infinite-bulk radiation coefficient or scattering state computed',
                       'Not a derivation of graph, masses or stiffnesses from homogeneous H18/DST',
                       'No finite-amplitude error bound for nonlinear radial parent',
                       'Scalar vector-line sector requires C4v-invariant prepared state; arbitrary bulk is excluded'],
        'source_hashes':{n:hashlib.sha256(data).hexdigest() for n,data in modules.items()}}
    (HERE/'vector_parent_check_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
