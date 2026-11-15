"""Independent static recheck of archived same-bulk packet proof; no trajectory."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import numpy as np
# Elements a+b sqrt(3), exact rational operations.
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def sub(x,y):return add(x,neg(y))
def mul(x,y):return (x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def q(a,b=1):return (F(a,b),F(0))
def mv(M,v):return [sumq([mul(a,b) for a,b in zip(row,v)]) for row in M]
def sumq(xs):
 out=q(0)
 for x in xs:out=add(out,x)
 return out
M=[[q(15,2),q(3,2),(F(0),F(1))],[q(3,2),q(11,2),(F(0),F(1))],[(F(0),F(1)),(F(0),F(1)),q(6)]]
v=[q(4),q(2),q(3)];w=mv(M,v);z=mv(M,w)
det=sumq([mul(v[i],sub(mul(w[(i+1)%3],z[(i+2)%3]),mul(w[(i+2)%3],z[(i+1)%3]))) for i in range(3)])
assert det==(F(432),F(-236)),det
# sqrt(3) irrational and coefficient nonzero => exact cyclic witness !=0.
a=[F(4),F(27,8),F(13,4)];vand=(a[1]-a[0])*(a[2]-a[0])*(a[2]-a[1]);assert vand==F(-15,256)
bonds=np.asarray([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
def symbol(k):
 return np.einsum('b,bi,bj->ij',4*np.sin(bonds@k/2)**2,bonds,bonds/np.sum(bonds*bonds,axis=1)[:,None])
def component_symbol(k):
 c=np.cos(k);s=np.sin(k);D=2*np.outer(s,s)
 np.fill_diagonal(D,6-2*c*(1+c.sum()-c))
 return D
max_symbol_residual=0.
port_radii=[.01,.05,.1,.2,.4]
port_directions=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,1,1),(1,2,3),(3,-2,7)]
port_cases=0
for radius in port_radii:
 for direction in port_directions:
  n=np.asarray(direction,dtype=float);n/=np.linalg.norm(n);k=radius*n
  D=symbol(k);residual=float(np.max(abs(D-component_symbol(k))));max_symbol_residual=max(max_symbol_residual,residual)
  assert residual<1e-12,residual
  lam,V=np.linalg.eigh(D);assert lam[2]-lam[1]>0
  assert np.linalg.norm(V[:,:2]@V[:,:2].T@n)<=radius*radius/4+2e-12
  port_cases+=1
rng=np.random.default_rng(61106);max_bound_ratio=0.;max_operator_defect=0.;min_gap=100.
for rho in np.geomspace(.001,1.8,20):
 for _ in range(16):
  n=rng.normal(size=3);n/=np.linalg.norm(n);k=rho*n;D=symbol(k);D2=rho*rho*(np.eye(3)+2*np.outer(n,n));E=D2-D
  lo,hi=np.linalg.eigvalsh(E)[[0,-1]]
  assert lo>=-1e-12 and hi<=rho**4/2+1e-12,(rho,lo,hi)
  lam,V=np.linalg.eigh(D);P=V[:,:2]@V[:,:2].T;leak=np.linalg.norm(P@n)
  assert leak<=rho*rho/4+2e-12,(rho,leak)
  max_bound_ratio=max(max_bound_ratio,float(leak/(rho*rho/4)))
  max_operator_defect=max(max_operator_defect,float(max(0,-lo,hi-rho**4/2)))
  min_gap=min(min_gap,float((lam[2]-lam[1])/(rho*rho)))
k=np.array([np.pi/2,np.pi/3,np.pi/3]);e=np.array([0.,1.,-1.])/np.sqrt(2);assert np.linalg.norm(symbol(k)@e-3*e)<1e-12 and abs(e@k)<1e-12
report={'status':'PASS','static_numeric_cases':320,'explicit_portability_cases':port_cases,'portability_radii':port_radii,'maximum_component_symbol_residual':max_symbol_residual,'exact_wave_normal_witness':{'rational':str(det[0]),'sqrt3_coefficient':str(det[1])},'exact_centered_difference_vandermonde':str(vand),'spectral_radius_range':[.001,1.8],'maximum_observed_to_archived_amplitude_bound_ratio':max_bound_ratio,'maximum_numerical_operator_defect':max_operator_defect,'minimum_scaled_lower_upper_gap':min_gap,'scope':'Independent static confirmation of pre-existing archived identities and sampled bound, not a proof of global novelty, trajectory propagation, or tight reference-frequency polarization certificate. Original packet_theorem checker not found in inspected restored readers/archives; this is new check code.', 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
