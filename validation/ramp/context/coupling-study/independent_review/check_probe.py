"""Independent physical-spring Hessian and coupled block-solve review.

Imports no archived solver. Reuses the declared bare Green kernels only.
No kernel generation, optimization, trajectory or manuscript write.
"""
from pathlib import Path
import hashlib, json, math
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
PROBE=ROOT/'numerical_probe'
OUT=Path(__file__).resolve().parent
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
INPUT=PROBE/'inputs'
expected={
 'model_api_v1_1.py':'aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980',
 'reduced_core.py':'0c6967e11c84d6c2fffeadc250c9d02f39374f30b1eade90d30660aa7144ed57',
 'REFERENCE_KERNEL.npz':'82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b',
}
for name,h in expected.items():assert H(INPUT/name)==h
meta=json.loads((INPUT/'REFERENCE_KERNEL.json').read_text())
assert meta['omega']==1.2 and meta['surface_grid']==192
assert meta['PV_selection']['Jackson_order']==256
assert meta['model_parameters']=={'a':1.,'m_bulk':1.,'m_spine':1.,'k_axis':1.,'k_diag':1.,'k_spine':1.69,'k_couple':.2}
with np.load(INPUT/'REFERENCE_KERNEL.npz') as z:
    G192=z['G0'][:12,:12]
    BT192=z['Im_branch_0'][:12,:12]+z['Im_branch_1'][:12,:12]
    BL192=z['Im_branch_2'][:12,:12]
with np.load(INPUT/'on_shell_N80.npz') as z:
    BT80=z['Im_branch_0'][:12,:12]+z['Im_branch_1'][:12,:12]
    BL80=z['Im_branch_2'][:12,:12]
G80=np.load(INPUT/'pv_mean_M192.npy')[:12,:12]+1j*(BT80+BL80)
kernels={'N192_M256':(G192,BT192,BL192),'N80_M192':(G80,BT80,BL80)}

def spring_blocks(g):
    """Assemble each of the actual eight central springs, then project."""
    n=len(g);nq=n+1;nb=12*n
    K=np.zeros((nq+nb,nq+nb));B=np.zeros((nb,2*n))
    signs=[(-1,-1),(-1,1),(1,-1),(1,1)]
    for s,gs in enumerate(g):
        for c,(sy,sz) in enumerate(signs):
            idx=nq+12*s+3*c
            B[idx-nq,2*s]=.5
            B[idx-nq+1,2*s+1]=sy/(2*math.sqrt(2))
            B[idx-nq+2,2*s+1]=sz/(2*math.sqrt(2))
            bulk=np.array([s+.5,sy*.5,sz*.5])
            for endpoint in [s,s+1]:
                direction=bulk-np.array([endpoint,0,0]);direction/=np.linalg.norm(direction)
                link=np.zeros(nq+nb);link[endpoint]=-direction[0];link[idx:idx+3]=direction
                K+=gs*np.outer(link,link)
    assert np.max(np.abs(B.T@B-np.eye(2*n)))<1e-14
    return B.T@K[nq:,nq:]@B,B.T@K[nq:,:nq],K[:nq,:nq]

records=[];states_saved={};worst=0.;case_count=0
sources=[PROBE/'SCALE_PROBE.json']
profile=PROBE/'profile_comparison/PROFILE_COMPARISON.json'
if profile.exists():sources.append(profile)
for source in sources:
    result=json.loads(source.read_text())
    rows=result['states']
    for row in rows:
        g=np.array(row['g']);assert len(g)==6 and np.min(g)>=0
        assert row.get('omega',result.get('omega'))==1.2
        n=6;nq=7;nb=12;omega=1.2;eta=1e-12;Kq=1.69
        G,BT,BL=kernels[row['cache']]
        U,V,Kself=spring_blocks(g)
        Kline=np.diag(np.full(nq,2*Kq))+np.diag(np.full(nq-1,-Kq),1)+np.diag(np.full(nq-1,-Kq),-1)
        z=complex(omega,eta)
        # Independent complex-dispersion evaluation of the outgoing lead root.
        klead=2*np.arcsin(z/(2*math.sqrt(Kq)))
        selflead=-Kq*np.exp(1j*klead);gamma=-2*selflead.imag
        assert selflead.imag<0 and gamma>0
        lead=np.zeros((nq,nq),complex);lead[0,0]=lead[-1,-1]=selflead
        Aq0=z*z*np.eye(nq)-Kline-Kself-lead
        block=np.block([[Aq0,V.conj().T@G],[-V,np.eye(nb)+U@G]])
        rhs=np.zeros(nq+nb,complex);rhs[0]=1
        solution=np.linalg.solve(block,rhs);q=solution[:nq];t=solution[nq:]
        # Convert unit-RHS Green response to unit-displacement incoming wave.
        physical_q=1j*gamma*q;physical_force=1j*gamma*t
        incident_power=omega*gamma/4
        values={
            'Rq':float(abs(physical_q[0]-1)**2),
            'Tline':float(abs(physical_q[-1])**2),
            'TqT':float((omega*.5*np.vdot(physical_force,BT@physical_force).real)/incident_power),
            'TqL':float((omega*.5*np.vdot(physical_force,BL@physical_force).real)/incident_power),
        }
        errors={k:abs(v-row[k]) for k,v in values.items()}
        largest=max(errors.values());worst=max(worst,largest)
        assert largest<2e-13,(row,values,errors)
        assert min(values.values())>=-1e-13
        closure=abs(sum(values.values())-1)
        assert closure<1e-9
        assert min(np.linalg.eigvalsh(BT))>=-1e-12 and min(np.linalg.eigvalsh(BL))>=-1e-12
        assert abs(8*sum(g)-row['physical_total_spring_stiffness'])<1e-12
        doubled_qt_closure=abs(sum(values.values())+values['TqT']-1)
        assert doubled_qt_closure>1e-5
        key=f'case{case_count}'
        for label,x in [('q',q),('t',t),('U',U),('V',V),('Kself',Kself)]:states_saved[key+'_'+label]=x
        case_count+=1
        records.append({'source':str(source.relative_to(PROBE)),'cache':row['cache'],'g':g.tolist(),'independent_partition':values,'maximum_difference':largest,'closure':closure,'double_qT_adverse_closure':doubled_qt_closure,'relative_block_residual':float(np.linalg.norm(block@solution-rhs)/(np.linalg.norm(block)*np.linalg.norm(solution)+np.linalg.norm(rhs)))})
OUT.mkdir(exist_ok=True)
np.savez_compressed(OUT/'INDEPENDENT_STATES.npz',**states_saved)
receipt={
 'status':'PASS','date':'2026-11-08','independent_case_count':case_count,
 'method':'Actual8-central-spring Hessian per active motif, independent orthonormal X/R projection and direct full coupled-block solve. No archived model/solver code imported. Bare cached Green matrices reused.',
 'incident_normalization':'Unit incident axial displacement: physical q=i*gamma*q_Green, physical bulk source=i*gamma*t_Green; P_inc=omega*gamma/4; P_qT=omega/2*fdagger*(Imbranch0+Imbranch1)*f. The branch sum already counts both lower modes.',
 'maximum_partition_disagreement':worst,
 'records':records,'result_source_sha256':{str(p.relative_to(PROBE)):H(p) for p in sources},
 'input_sha256':{k:H(INPUT/k) for k in expected},
 'scope':'Independent finite-system algebra/normalization review only. Both routes share archived approximate bare Green inputs. No new quadrature, rigorous kernel error enclosure, bandwidth, global optimum, physical fabrication or experimental confirmation.',
 'confirmed_issues':[],
}
(OUT/'REVIEW.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','independent_case_count','maximum_partition_disagreement','confirmed_issues']},indent=2))
