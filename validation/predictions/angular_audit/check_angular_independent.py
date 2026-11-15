"""Independent bond-symbol and literal Fourier finite-cone check.
Inputs: frozen dressed sources rebuilt by the archived source core. No archived
symbol, root, gradient, density, symmetry sector, or cone integrator is used.
The new integration uses a full 2pi midpoint azimuth grid and scalar brent roots.
It verifies observed local connected cone patches, not global preimage completeness.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.optimize import brentq
BASE=Path(__file__).resolve().parent.parent/'baseline/Prediction_Calculations_2026-09-26/angular'
sys.path.insert(0,str(BASE))
from core import states
OUT=Path(__file__).resolve().parent
src,ref=states();eta=src['full']
bonds=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],float)
e=bonds/np.linalg.norm(bonds,axis=1)[:,None];bb=np.einsum('bi,bj->bij',e,e)
pos=[];force=[]
for seg in range(6):
 for sy in [-1,1]:
  for sz in [-1,1]:
   pos.append([seg,sy/2,sz/2]);force.append([eta[2*seg]/2,sy*eta[2*seg+1]/(2*np.sqrt(2)),sz*eta[2*seg+1]/(2*np.sqrt(2))])
pos=np.array(pos);force=np.array(force)
def root(theta,phi,branch,sign):
 n=np.array([sign*np.cos(theta),np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi)])
 def calc(r):
  ph=bonds@(r*n);D=np.einsum('b,bij->ij',2*(1-np.cos(ph)),bb);return np.linalg.eigh(D)
 r=brentq(lambda r:calc(r)[0][branch]-1.44,0,1.32/np.max(abs(n)),xtol=3e-15)
 vals,ev=calc(r);v=ev[:,branch];k=r*n
 grad=np.einsum('b,ba,b->a',2*np.sin(bonds@k),bonds,(e@v)**2)
 f=np.exp(-1j*(pos@k))@force
 rho=r*r*abs(v@f)**2/(4*np.pi*np.pi*(grad@n))
 return rho,grad

def cones(angle,nphi=64,nmu=28):
 th=np.deg2rad(angle);threshold=np.cos(th);u,w=np.polynomial.legendre.leggauss(nmu);ans={}
 # Full circle; offset prevents exact degenerate planes and dihedral diagonals.
 phis=(np.arange(nphi)+.371)*2*np.pi/nphi
 for sign,side in [(1,'forward'),(-1,'backward')]:
  terms=[]
  for j in range(3):
   total=0.
   for ph in phis:
    def fun(t):
     _,g=root(t,ph,j,sign);return sign*g[0]/np.linalg.norm(g)-threshold
    boundary=brentq(fun,1e-6,max(.15,th*2),xtol=2e-13)
    lo=np.cos(boundary);mus=lo+(u+1)*(1-lo)/2
    for mu,wt in zip(mus,w*(1-lo)/2):total+=wt*root(np.arccos(mu),ph,j,sign)[0]*2*np.pi/nphi
   terms.append(total)
  ans[side]={'qT':sum(terms[:2]),'qL':terms[2]}
 return ans

aT=.72;cT=.28;kT=np.arccos(cT);qT=2*np.sin(kT);detT=4-aT*aT
cL=.76;kL=np.arccos(cL);qL=6*np.sin(kL);bL=1+2*cL
coef={}
for sign,side in [(1,'forward'),(-1,'backward')]:
 def amps(k):return np.exp(-1j*k*np.arange(6))@eta.reshape(6,2)
 ax,ar=amps(sign*kT);C=-1/np.tan(sign*kT/2)*ax-1j/np.sqrt(2)*ar
 F=2*amps(sign*kL)[0]
 AT=abs(C)**2*qT**3/(16*np.pi*detT)*(2/(2+aT)**2)
 AL=abs(F)**2*qL/(16*np.pi*bL*bL)
 coef[side]={'qT_theta4':float(AT),'qL_theta2':float(AL),'share_theta2':float(AT/AL),'qT_amplitude_C_abs2':float(abs(C)**2),'qL_axial_f_abs2':float(abs(F)**2)}
rows=[]
for a in [.25,.5,1,10,20,30]:
 val=cones(a)
 for side in val:
  val[side]['qT_relative_to_local_law']=val[side]['qT']/(coef[side]['qT_theta4']*np.deg2rad(a)**4)-1
  val[side]['qL_relative_to_local_law']=val[side]['qL']/(coef[side]['qL_theta2']*np.deg2rad(a)**2)-1
 rows.append({'half_angle_deg':a,**val})
 print(json.dumps(rows[-1]),flush=True)
fine=json.loads((BASE/'integral_fine.json').read_text());comparisons=[]
for row in rows:
 a=row['half_angle_deg']
 if a==.25:continue
 for side in ['forward','backward']:
  old=fine['cone_branch_powers'][f'group_{side}_{a:g}']['full']
  comparisons.append({'angle':a,'side':side,'qT_relative_difference':row[side]['qT']/sum(old[:2])-1,'qL_relative_difference':row[side]['qL']/old[2]-1})
out={'method':__doc__,'baseline_core_sha256':hashlib.sha256((BASE/'core.py').read_bytes()).hexdigest(),'local_coefficients':coef,'cones':rows,'archive_comparison':comparisons}
(OUT/'INDEPENDENT_ANGULAR_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print('DONE')
