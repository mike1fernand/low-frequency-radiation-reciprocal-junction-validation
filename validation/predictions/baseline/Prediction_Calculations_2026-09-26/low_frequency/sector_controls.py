#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import low_frequency as m

def sector(g,kind):
 U,V,Ks,L=m.matrices(g)
 # Replace the complete local Ks; U,V have disjoint X/R row sectors.
 Ks[:]=0
 if kind=='X':U[1::2,:]=0;V[1::2]=0;loc=np.array([[1.,1.],[1.,1.]])
 else:U[0::2,:]=0;V[0::2]=0;loc=np.array([[1.,-1.],[-1.,1.]])
 for i,x in enumerate(g):Ks[i:i+2,i:i+2]+=2*x/3*loc
 return U,V,Ks,L

out=[]
for kind in ['X','R']:
 U,V,Ks,L=sector(m.GREF,kind)
 z=np.load(m.ROOT/'validation/kernel_0.00000000.npz');P=z['PV'];L0=L.copy();L0[0,0]-=m.KQ;L0[-1,-1]-=m.KQ
 static=L0+Ks-V.T@np.linalg.solve(np.eye(12)+P@U,P)@V
 for w in [.2,.05,.0125,.003125,.00078125]:
  z=np.load(m.ROOT/f'validation/kernel_{w:.8f}.npz');P=z['PV'];bs=[z[f'Q{j}'] for j in range(3)];G=P+1j*sum(bs)
  kap=2*np.arcsin(w/(2*m.CQ));gam=2*m.KQ*np.sin(kap);gc=np.linalg.solve(np.eye(12)+G@U,G)
  A=w*w*np.eye(7)-L-Ks+V.T@gc@V;A[0,0]+=m.KQ*np.exp(1j*kap);A[-1,-1]+=m.KQ*np.exp(1j*kap)
  chi=np.linalg.solve(A,np.eye(7)[:,0]*np.sqrt(gam));eta=np.linalg.solve(np.eye(12)+U@G,V@chi)
  T=[float(2*np.real(eta.conj()@b@eta)) for b in bs];R=abs(1-1j*np.sqrt(gam)*chi[0])**2;line=gam*abs(chi[-1])**2
  out.append({'kind':kind,'omega':w,'static_min_eigenvalue':float(np.linalg.eigvalsh(static).min()),'TqT':T[0]+T[1],'TqL':T[2],'TqT_over_omega2':(T[0]+T[1])/w**2,'TqT_over_omega4':(T[0]+T[1])/w**4,'qT_share':(T[0]+T[1])/sum(T),'R':R,'Tline':line,'closure':R+line+sum(T)-1})
m.dump(m.ROOT/'sector_controls.json',out)
G=np.load(m.ROOT/'validation/kernel_0.00000000.npz')['PV']
C=G[1::2,1::2];S=np.diag(4*m.GREF/3)
W=.5*np.linalg.inv(np.linalg.inv(S)+C)
h=np.ones(6)@W@np.linalg.solve(m.KQ*np.eye(6)+W,np.ones(6))
M2=2*m.CQ*h*h
pred={'h':float(h),'M_squared':float(M2),
      'TqT_over_omega4':float(M2/(60*np.pi)),
      'TqL_over_omega4':float(M2/(15*np.pi*3**2.5)),
      'qT_share':float(3**2.5/(4+3**2.5))}
m.dump(m.ROOT/'sector_predictions.json',pred)
checks=[]
U,V,Ks,L=sector(m.GREF,'R')
for w in [.05,.0125,.003125,.00078125]:
 vals=[]
 for tag in ['refined','validation']:
  z=np.load(m.ROOT/tag/f'kernel_{w:.8f}.npz')
  bs=[z[f'Q{j}'] for j in range(3)];G=z['PV']+1j*sum(bs)
  kap=2*np.arcsin(w/(2*m.CQ));gam=2*m.KQ*np.sin(kap)
  A=w*w*np.eye(7)-L-Ks+V.T@np.linalg.solve(np.eye(12)+G@U,G)@V
  A[0,0]+=m.KQ*np.exp(1j*kap);A[-1,-1]+=m.KQ*np.exp(1j*kap)
  chi=np.linalg.solve(A,np.eye(7)[:,0]*np.sqrt(gam))
  eta=np.linalg.solve(np.eye(12)+U@G,V@chi)
  vals.append([float(2*np.real(eta.conj()@(bs[0]+bs[1])@eta)),
               float(2*np.real(eta.conj()@bs[2]@eta))])
 checks.append({'omega':w,
  'qT_grid_relative_difference':abs(vals[1][0]/vals[0][0]-1),
  'qL_grid_relative_difference':abs(vals[1][1]/vals[0][1]-1),
  'qT_coefficient_relative_error':vals[1][0]/w**4/pred['TqT_over_omega4']-1,
  'qL_coefficient_relative_error':vals[1][1]/w**4/pred['TqL_over_omega4']-1})
m.dump(m.ROOT/'sector_refinement.json',checks)
print(pred)
