"""Check the P3 selected-band/P1 frequency distinction, without scattering runs."""
from pathlib import Path
import hashlib, importlib.util, json, math, sys
from fractions import Fraction
import numpy as np

sys.dont_write_bytecode = True
HERE=Path(__file__).resolve().parent
MODEL=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/model_api_v1_1.py')
P3=HERE.parent/'P3_CURRENT.txt'
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert H(MODEL)=='aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980'
before=H(MODEL)
spec=importlib.util.spec_from_file_location('bound_model_api',MODEL)
model=importlib.util.module_from_spec(spec);sys.modules[spec.name]=model;spec.loader.exec_module(model)
p=model.ModelParameters()
assert p.a==p.m_bulk==p.k_axis==p.k_diag==1
bonds=np.asarray(model.AXIAL_BONDS+model.FACE_DIAGONAL_BONDS,float)
expected=np.asarray([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
assert np.array_equal(bonds,expected)

# Independent exact rational coefficients for the six quadratic monomials,
# retaining the normalization |b|^-2. No symbolic dependency is needed.
pairs=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
C=[[[Fraction(0) for _ in range(3)] for _ in range(3)] for _ in pairs]
target=np.zeros((6,3,3),int)
for m,(i,j) in enumerate(pairs):
    if i==j:
        target[m]=np.eye(3,dtype=int);target[m,i,i]+=2
    else:target[m,i,j]=target[m,j,i]=2
for b in expected.astype(int):
    denom=int(b@b)
    for m,(i,j) in enumerate(pairs):
        coef=int(b[i]*b[j])*(1 if i==j else 2)
        for a in range(3):
            for c in range(3):C[m][a][c]+=Fraction(coef*int(b[a]*b[c]),denom)
assert all(C[m][a][c]==int(target[m,a,c]) for m in range(6) for a in range(3) for c in range(3))

def independent_D(k,normalize=True):
    D=np.zeros((3,3))
    for b in bonds:
        D+=2*(1-math.cos(float(k@b)))*np.outer(b,b)/(b@b if normalize else 1)
    return D

rng=np.random.default_rng(20261108)
points=list(rng.uniform(-math.pi,math.pi,size=(4096,3)))
for x in [-math.pi,-1.,-.01,0.,.01,1.,math.pi]:
 for y in [-math.pi,-1.,-.01,0.,.01,1.,math.pi]:
  for z in [-math.pi,-1.,-.01,0.,.01,1.,math.pi]:points.append(np.array([x,y,z]))
dirs=rng.normal(size=(128,3));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
bandpoints=[rho*v for rho in [1e-5,.001,.01,.05,.2,.5,.9,1.] for v in dirs]
points+=bandpoints
api_error=0.;loewner_min=math.inf;lower_excess=-math.inf;lower_band_max=0.
for k in points:
    rho2=float(k@k);D=independent_D(k);D2=rho2*np.eye(3)+2*np.outer(k,k)
    api_error=max(api_error,float(np.max(np.abs(D-model.bulk_dynamical_matrix(k,p)))))
    loewner_min=min(loewner_min,float(np.linalg.eigvalsh(D2-D).min()))
    lower_excess=max(lower_excess,float(np.linalg.eigvalsh(D)[1]-rho2))
for k in bandpoints:lower_band_max=max(lower_band_max,float(np.linalg.eigvalsh(independent_D(k))[1]))
assert api_error<1e-12
assert loewner_min>-1e-12
assert lower_excess<1e-12
assert lower_band_max<=1+1e-12

# Explicit adverse controls show why branch choice, unit springs and normalized
# face-diagonal vectors are load-bearing rather than optional conventions.
k=np.array([.1,0.,0.]);rho2=float(k@k)
adverse={
 'upper_branch_not_bounded_by_rho_squared':float(np.linalg.eigvalsh(independent_D(k))[2]-rho2),
 'omitting_unit_bond_normalization_breaks_coefficient':float(np.linalg.eigvalsh(independent_D(k,False))[1]-rho2),
 'quadrupling_spring_stiffness_requires_rescaling':float(np.linalg.eigvalsh(model.bulk_dynamical_matrix(k,model.ModelParameters(k_axis=4,k_diag=4)))[1]-rho2),
}
assert min(adverse.values())>0

# The continuum Maxwell readouts re-express the one existing two-polarization
# energy. Both electric-like and magnetic-like terms are part of that energy.
max_energy_error=0.;doubling_rejected=True
for k in dirs:
    rho=rng.uniform(.01,1);k=rho*k;P=np.eye(3)-np.outer(k,k)/(k@k)
    v=P@(rng.normal(size=3)+1j*rng.normal(size=3))
    vt=P@(rng.normal(size=3)+1j*rng.normal(size=3))
    E=-vt;B=1j*np.cross(k,v)
    HT=.5*(np.vdot(vt,vt).real+(k@k)*np.vdot(v,v).real)
    HEB=.5*(np.vdot(E,E).real+np.vdot(B,B).real)
    err=abs(HT-HEB)/max(HT,1e-30)
    max_energy_error=max(max_energy_error,float(err))
    doubling_rejected &= abs(2*HEB-HT)>1e-8*HT
assert max_energy_error<1e-12 and doubling_rejected
assert H(MODEL)==before

out={
 'status':'PASS','date':'2026-11-08',
 'scope':'Analytic bond-sum/Loewner deduction with bounded matrix and continuum-energy checks. No Green-kernel, scattering, optimizer or trajectory run.',
 'units':{'a':p.a,'mass':p.m_bulk,'axial_stiffness':p.k_axis,'face_diagonal_stiffness':p.k_diag},
 'model_api_sha256':H(MODEL),'P3_text_sha256':H(P3),
 'bond_families':expected.astype(int).tolist(),
 'exact_symbolic_quadratic_sum':'|k|^2 I + 2 k k^T',
 'analytic_claim':'0 <= D(k) <= |k|^2 I + 2 k k^T in Loewner order; ascending lambda0,lambda1 <= |k|^2.',
 'frequency_consequence':'At omega=1.2, either lower-branch shell requires |k|>=1.2 and is disjoint from P3 selected 0<|k|<=1. This concerns the lower pair, not upper branch.',
 'numerical_points':len(points),'selected_band_points':len(bandpoints),
 'api_max_abs_difference':api_error,'minimum_D2_minus_D_eigenvalue':loewner_min,
 'maximum_lower_eigenvalue_minus_rho_squared':lower_excess,
 'largest_checked_selected_band_lower_eigenvalue':lower_band_max,
 'adverse_controls':adverse,'continuum_transverse_energy_cases':128,
 'continuum_energy_max_relative_difference':max_energy_error,
 'extra_factor_two_rejected':bool(doubling_rejected),
 'important_limits':['Bound does not invalidate P1 lower-cluster qT power outside P3 selected band.',
 'No efficiency improvement follows from relabelling selected modes or Maxwell E/B readouts.',
 'Equal stiffness/mass/spacing units are necessary for the numerical threshold1; with commonK,m and dimensionless theta=ak, lower omega<=sqrt(K/m)|theta|.',
 'Maxwell energy identity is the leading continuum transverse subsystem, not exact finite-k Dy/current.'],
}
(HERE/'FREQUENCY_BOUND_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
