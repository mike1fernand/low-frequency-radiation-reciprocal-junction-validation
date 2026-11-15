"""Exact bond-Hessian check and four exploratory passive endpoint-imbalance states.

Explicitly extends the equal-endpoint contact class. The same old approximate
kernels are reused; no quadrature, optimization, new bulk or trajectory is run.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import io
import json
import math
import sys
import tempfile
import numpy as np
import check_purity_relaxation as cached


class Q2:
    """Exact a+b sqrt(2), only for fixed coefficient identities."""
    def __init__(self,a=0,b=0):
        self.a,self.b=F(a),F(b)
    @staticmethod
    def coerce(x):return x if isinstance(x,Q2) else Q2(x)
    def __add__(self,x):
        x=self.coerce(x);return Q2(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return Q2(-self.a,-self.b)
    def __sub__(self,x):return self+-self.coerce(x)
    def __mul__(self,x):
        x=self.coerce(x);return Q2(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __eq__(self,x):
        x=self.coerce(x);return self.a==x.a and self.b==x.b
    def text(self):return {'rational':str(self.a),'sqrt2':str(self.b)}


SIGNS=[(-1,-1),(-1,1),(1,-1),(1,1)]


def zero(n,m):return [[Q2() for _ in range(m)] for _ in range(n)]


def exact_bond_hessian(beta):
    """g=1; integer bond extensions have squared unit normalization 1/3."""
    H=zero(14,14)
    for corner,(sy,sz) in enumerate(SIGNS):
        for endpoint,nx,stiff in [(0,1,1+F(beta)),(1,-1,1-F(beta))]:
            c=[0]*14;c[endpoint]=-nx
            c[2+3*corner:5+3*corner]=[nx,sy,sz]
            for i in range(14):
                for j in range(14):H[i][j]+=F(stiff,3)*c[i]*c[j]
    return H


def exact_energy_hessian(beta):
    H=zero(14,14)
    for corner,(sy,sz) in enumerate(SIGNS):
        A=[F(0)]*14;B=[F(0)]*14
        A[0]=A[1]=F(-1,2);A[2+3*corner]=1
        B[0]=F(-1,2);B[1]=F(1,2)
        B[3+3*corner]=sy;B[4+3*corner]=sz
        for i in range(14):
            for j in range(14):
                H[i][j]+=F(2,3)*(A[i]*A[j]+B[i]*B[j]+F(beta)*(A[i]*B[j]+B[i]*A[j]))
    return H


def exact_embedding():
    C=zero(14,4);C[0][0]=C[1][1]=Q2(1)
    for corner,(sy,sz) in enumerate(SIGNS):
        C[2+3*corner][2]=Q2(F(1,2))
        C[3+3*corner][3]=Q2(0,F(sy,4))
        C[4+3*corner][3]=Q2(0,F(sz,4))
    return C


def project(H,C):
    R=zero(4,4)
    for a in range(4):
        for b in range(4):
            R[a][b]=sum(C[i][a]*H[i][j]*C[j][b] for i in range(14) for j in range(14))
    return R


def expected_local(beta):
    b=F(beta);R=zero(4,4)
    R[0][0]=Q2(F(4,3)*(1+b));R[1][1]=Q2(F(4,3)*(1-b))
    R[2][2]=Q2(F(2,3));R[3][3]=Q2(F(4,3))
    R[2][3]=R[3][2]=Q2(0,F(2,3)*b)
    R[2][0]=R[0][2]=Q2(F(-2,3)*(1+b))
    R[2][1]=R[1][2]=Q2(F(-2,3)*(1-b))
    R[3][0]=R[0][3]=Q2(0,F(-2,3)*(1+b))
    R[3][1]=R[1][3]=Q2(0,F(2,3)*(1-b))
    return R


def exact_checks():
    C=exact_embedding()
    for a in range(4):
        for b in range(4):assert sum(C[i][a]*C[i][b] for i in range(14))==(1 if a==b else 0)
    coefficients=[]
    # H is affine in beta and linear in g. Checking its beta^0 and beta^1
    # coefficient identities proves these formulas for all real beta and g.
    for beta in [F(0),F(1)]:
        H=exact_bond_hessian(beta)
        assert H==exact_energy_hessian(beta)
        R=project(H,C);assert R==expected_local(beta)
        # The complete invariant line/X/R subspace has no Hessian leakage.
        HC=[[sum(H[i][j]*C[j][a] for j in range(14)) for a in range(4)] for i in range(14)]
        CR=[[sum(C[i][b]*R[b][a] for b in range(4)) for a in range(4)] for i in range(14)]
        assert HC==CR
        coefficients.append([[v.text() for v in row] for row in R])
    return {'arithmetic':'Exact rational coefficients in Q(sqrt(2))',
       'energy_vs_eight_physical_bonds':'PASS','orthonormal_embedding':'PASS',
       'projected_local_matrix':'PASS','source_sector_Hessian_leakage':'exact zero',
       'universal_identity_reason':'Each entry is linear in g and affine in beta; exact coefficient identity at beta=0 and beta=1.',
       'passivity_domain':'g>=0 and |beta|<=1 by sum of nonnegative spring energies; all requested points have g>0, |beta|<1.',
       'local_order':['q_left','q_right','X','R'],'beta0_matrix':coefficients[0],
       'beta1_matrix':coefficients[1]}


def direct_local(g,beta):
    H=np.zeros((14,14));C=np.zeros((14,4));C[0,0]=C[1,1]=1
    for corner,(sy,sz) in enumerate(SIGNS):
        C[2+3*corner,2]=.5;C[3+3*corner,3]=sy/(2*math.sqrt(2));C[4+3*corner,3]=sz/(2*math.sqrt(2))
        for endpoint,nx,stiff in [(0,1,g*(1+beta)),(1,-1,g*(1-beta))]:
            c=np.zeros(14);c[endpoint]=-nx;c[2+3*corner:5+3*corner]=[nx,sy,sz]
            H+=stiff*np.outer(c,c)/3
    return C.T@H@C


def matrices(g,beta):
    n=len(g);U=np.zeros((2*n,2*n),complex);V=np.zeros((2*n,n+1),complex);K=np.zeros((n+1,n+1))
    for s,gs in enumerate(g):
        h=direct_local(gs,beta)
        U[2*s:2*s+2,2*s:2*s+2]=h[2:,2:]
        V[2*s:2*s+2,s:s+2]=h[2:,:2]
        K[s:s+2,s:s+2]+=h[:2,:2]
    return U,V,K


def coupled(omega,g,beta,p,G0,lead_function):
    n=len(g);nq=n+1;U,V,Kself=matrices(g,beta)
    Gc=np.linalg.solve(np.eye(2*n)+G0@U,G0)
    Sigma=-(V.conj().T@Gc@V);Gamma=1j*(Sigma-Sigma.conj().T)
    Kq=np.diag(np.full(nq,2*p.k_spine))
    for i in range(nq-1):Kq[i,i+1]=Kq[i+1,i]=-p.k_spine
    z=complex(omega,1e-12);sl=lead_function(z,p)
    SL=np.zeros((nq,nq),complex);SR=np.zeros_like(SL);SL[0,0]=sl;SR[-1,-1]=sl
    GL=1j*(SL-SL.conj().T);GR=1j*(SR-SR.conj().T);gin=float(GL[0,0].real)
    A=p.m_spine*z*z*np.eye(nq)-Kq-Kself-Sigma-SL-SR
    rhs=np.zeros(nq,complex);rhs[0]=1;q=np.linalg.solve(A,rhs)
    f=np.linalg.solve(np.eye(2*n)+U@G0,V@q)
    return {'Rq':float(abs(1-1j*gin*q[0])**2),'Tq':float(gin*np.vdot(q,GR@q).real),
      'Tbulk':float(gin*np.vdot(q,Gamma@q).real),'q':q,'f':f,'gin':gin,
      'minimum_loaded_escape_eigenvalue':float(np.linalg.eigvalsh((Gamma+Gamma.conj().T)/2).min()),
      'transpose_symmetry_of_dressed_device':float(np.linalg.norm(A-A.T))}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',type=Path,default=Path(__file__).with_name('purity_inputs'))
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists: choose a fresh filename')
    exact=exact_checks()  # Must pass BEFORE any coupled state is evaluated.
    payloads,provenance,origin=cached.inputs_at(args.inputs,None)
    with np.load(io.BytesIO(payloads['on_shell_N80.npz']),allow_pickle=False) as d:
        BT=d['Im_branch_0'][:12,:12]+d['Im_branch_1'][:12,:12];BL=d['Im_branch_2'][:12,:12];B=d['Im_total'][:12,:12]
    assert np.linalg.norm(BT+BL-B)<1e-12
    G=np.load(io.BytesIO(payloads['pv_mean_M192.npy']),allow_pickle=False)[:12,:12]+1j*B
    assert np.linalg.norm(G-G.T)<1e-12
    reference=np.array([.1,.2,.2,.2,.2,.1]);assert abs(reference.sum()-1)<1e-14
    results={'status':'PASS','contact_class':'Explicit unequal-endpoint spring-orbit extension, not original equal-pair hypothesis',
      'omega':1.2,'g_reference':reference.tolist(),'declared_beta_values':[-.5,-.2,.2,.5],
      'input_origin':origin,'input_provenance':provenance,'exact_hessian_checks':exact,'new_states':[]}
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-endpoint-originals-') as td:
        stage=Path(td)
        for name in ['model_api_v1_1.py','reduced_core.py']:(stage/name).write_bytes(payloads[name])
        sys.path.insert(0,str(stage))
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices,scatter_from_reduced_green,analytic_q_lead_self_energy
        p=ModelParameters();U0,V0,K0=matrices(reference,0);Ua,Va,Ka=reduced_static_matrices(reference,p)
        matrix_error=max(float(np.max(np.abs(x-y))) for x,y in [(U0,Ua),(V0,Va),(K0,Ka)])
        assert matrix_error<1e-13,matrix_error
        old,old_state=scatter_from_reduced_green(1.2,reference,p,G,return_internal=True)
        control=coupled(1.2,reference,0,p,G,analytic_q_lead_self_energy)
        state_error=max(float(np.max(np.abs(control['q']-old_state['q']))),float(np.max(np.abs(control['f']-old_state['t']))),
                        *(abs(control[k]-getattr(old,k)) for k in ['Rq','Tq','Tbulk']))
        assert state_error<1e-12,state_error
        baseline_T=float(2*old_state['gamma_in']*np.vdot(old_state['t'],BT@old_state['t']).real)
        baseline_L=float(2*old_state['gamma_in']*np.vdot(old_state['t'],BL@old_state['t']).real)
        results['beta0_archived_control']={'maximum_matrix_error':matrix_error,'maximum_state_or_fraction_error':state_error,
         'incident_qT':baseline_T,'incident_qL':baseline_L,'bulk_qT_share':baseline_T/(baseline_T+baseline_L)}
        for beta in results['declared_beta_values']:
            left=reference*(1+beta);right=reference*(1-beta)
            assert np.all(left>0) and np.all(right>0)
            total=float(4*(left.sum()+right.sum()));assert abs(total-8)<1e-13
            x=coupled(1.2,reference,beta,p,G,analytic_q_lead_self_energy)
            qT=float(2*x['gin']*np.vdot(x['f'],BT@x['f']).real);qL=float(2*x['gin']*np.vdot(x['f'],BL@x['f']).real)
            component_error=abs(qT+qL-x['Tbulk']);flux_error=abs(x['Rq']+x['Tq']+qT+qL-1)
            assert component_error<1e-12 and flux_error<1e-9
            assert x['minimum_loaded_escape_eigenvalue']>-1e-9 and x['transpose_symmetry_of_dressed_device']<1e-12
            assert min(qT,qL,x['Rq'],x['Tq'])>-1e-9
            results['new_states'].append({'beta':beta,'left_spring_strengths':left.tolist(),'right_spring_strengths':right.tolist(),
             'total_physical_spring_stiffness':total,'all_48_contact_springs_positive':True,
             'Rq':x['Rq'],'Tq':x['Tq'],'incident_qT':qT,'incident_qL':qL,'bulk_qT_share':qT/(qT+qL),
             'bulk_power':x['Tbulk'],'four_channel_closure':flux_error,'branch_vs_escape_power_error':component_error,
             'minimum_loaded_escape_eigenvalue':x['minimum_loaded_escape_eigenvalue'],
             'transpose_symmetry_of_dressed_device':x['transpose_symmetry_of_dressed_device']})
    pairs=[(results['new_states'][0],results['new_states'][3]),(results['new_states'][1],results['new_states'][2])]
    mirror_transmission_error=max(abs(a['Tq']-b['Tq']) for a,b in pairs)
    assert mirror_transmission_error<1e-12
    results['mirror_and_reciprocity_control']={'line_transmission_beta_sign_error':mirror_transmission_error,
      'scope':'For mirror-symmetric reference g and identical leads, reflection plus reciprocity predicts equal line transmission at beta and -beta; one-sided bulk conversion need not agree.'}
    results['limits']=['Four predeclared new cached-matrix states after exact local-Hessian and beta0 archived checks; no optimization.',
      'Same old N80 branch/PV kernels; new-profile numerical quadrature error is NOT enclosed.',
      'Endpoint imbalance is an explicit positive-spring contact extension; original equal-endpoint theorem/certificate not extended.',
      'Reciprocity is retained. No nonreciprocity or two independently driven input amplitudes are claimed.',
      'No physical material realization, electromagnetic identification, or field-wide novelty established.']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(results,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'status':'PASS','output':str(args.output),'exact_hessian_checks':'PASS','beta0_matrix_error':matrix_error,
        'beta0_state_error':state_error,'new_states':results['new_states']},indent=2))


if __name__=='__main__':main()
