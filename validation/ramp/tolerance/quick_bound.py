from pathlib import Path
import itertools,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
INPUT=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs')
sys.path.insert(0,str(INPUT))
from model_api_v1_1 import ModelParameters
from reduced_core import reduced_static_matrices,analytic_q_lead_self_energy
p=ModelParameters();n=6;kstar=p.k_spine;omega=1.2;eta=1e-12
z=complex(omega,eta);sig=analytic_q_lead_self_energy(z,p);gamma=-2*sig.imag
with np.load(INPUT/'REFERENCE_KERNEL.npz') as f:G=f['G0'][:12,:12].copy();Q=f['Im_branch_0'][:12,:12]+f['Im_branch_1'][:12,:12]
Kq=np.diag(np.full(7,2*p.k_spine));Kq+=np.diag(np.full(6,-p.k_spine),1)+np.diag(np.full(6,-p.k_spine),-1)
SL=np.zeros((7,7),complex);SL[0,0]=SL[-1,-1]=sig
GREF=np.arange(1,7)/7

def matrices(g):
 U,V,K=reduced_static_matrices(g,p)
 h=z*z*np.eye(7)-Kq-K-SL
 A=np.block([[h/kstar,V.conj().T@G],[-V/kstar,np.eye(12)+U@G]])
 D=np.vstack([V.conj().T,U])/kstar
 return A,U,V,D
A,U0,V,D0=matrices(GREF);R=np.linalg.inv(A);W0=np.eye(12)+U0@G;W0inv=np.linalg.inv(W0)
b=np.zeros(19,complex);b[0]=math.sqrt(gamma/kstar);y0=R@b
qlam,qv=np.linalg.eigh(Q);F=np.zeros((12,19),complex);F[:,7:]=math.sqrt(2*kstar)*(qv*np.sqrt(np.maximum(qlam,0)))@qv.conj().T
norm=np.linalg.norm
T0=norm(F@y0)**2

def compute(e):
 vals={x:0. for x in ['r_c','r','v','w','d_out','h_g','h_c','d_eta','w_eta']}
 for signs in itertools.product([-1,1],repeat=6):
  Av,Uv,Vv,Dv=matrices(GREF*(1+e*np.array(signs)))
  B=R@(Av-A)
  vals['r_c']=max(vals['r_c'],norm(W0inv@(Uv-U0)@G,2))
  vals['r']=max(vals['r'],norm(B,2));vals['v']=max(vals['v'],norm(B@y0))
  vals['w']=max(vals['w'],norm(F@B,2));vals['d_out']=max(vals['d_out'],norm(F@B@y0))
  vals['h_g']=max(vals['h_g'],norm(R@Dv,2));vals['h_c']=max(vals['h_c'],norm(W0inv@Uv/kstar,2))
  vals['d_eta']=max(vals['d_eta'],norm((B@y0)[7:]));vals['w_eta']=max(vals['w_eta'],norm(B[7:,:],2))
 r=vals['r'];L=math.sqrt(T0)-vals['d_out']-vals['w']*vals['v']/(1-r)
 Z=norm(y0[7:])+vals['d_eta']+vals['w_eta']*vals['v']/(1-r)
 eG=.02;eQ=.005;rho=r+vals['h_g']*eG;sigma=vals['h_g']*eG*Z/(1-rho)
 floor=(max(0,L-norm(F,2)*sigma))**2-2*eQ*(Z+sigma)**2
 upper=(math.sqrt(T0)+vals['d_out']+vals['w']*vals['v']/(1-r)+norm(F,2)*sigma)**2+2*eQ*(Z+sigma)**2
 return {'e':e,'T0':T0,**vals,'L':L,'Z':Z,'rho':rho,'sigma':sigma,'F_norm':norm(F,2),'finite_floor':max(0,L)**2,'infinite_floor_conditional':floor,'infinite_upper_conditional':upper,'core_total':vals['r_c']+vals['h_c']*eG}
result=[compute(0),compute(.05)]
Path(__file__).with_name('QUICK_BOUND.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
