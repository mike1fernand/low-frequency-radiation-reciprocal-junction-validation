"""Independent rational shell exclusion by distance to 2πZ on a different mesh; supplemental numerical spot checks."""
import json, math, itertools, sys
from functools import cache
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
B = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
# Actual period is contained in this rational interval, and our cube is wider than the true zone.
mults=[(-2*F(355,113),-2*F(333,106)),(F(0),F(0)),(2*F(333,106),2*F(355,113))]
M=10**7

def distance(a,b,c,d):return max(0,c-b,a-d)

def cert(sub,threshold):
 step=F(63,20)/sub
 @cache
 def weight(i,j):
  # Fold trigonometric intervals to distance from the nearest periodic zero, 0 <= d <= pi.
  d=min(distance(step*i,step*j,x,y) for x,y in mults)
  # The rationally computed d is <= the min genuine distance to the zero lattice; use degree 12 underbound for 2(1-cos).
  assert 0<=d<=F(355,113)
  f=2*sum((-1)**(n+1)*d**(2*n)/math.factorial(2*n) for n in range(1,7));return max(0,(f*M).__floor__())
 Q=2*M*threshold
 assert Q.denominator==1
 interior=count=0
 # Independent assembly from nine literal bond vectors, using reversed principal-minor order.
 for cell in itertools.combinations_with_replacement(range(sub),3):
  if sum((i+1)**2 for i in cell)*step**2 < F(9,4):interior+=1;continue
  L=[[0,0,0],[0,0,0],[0,0,0]]
  for b in B:
   t=sum(u*v for u,v in zip(cell,b)); lo=t+sum(min(z,0) for z in b); hi=t+sum(max(z,0) for z in b)
   p=weight(lo,hi)*(2//sum(v*v for v in b))
   for i in range(3):
    for j in range(i,3): L[i][j]+=p*b[i]*b[j]
  for i in range(3):L[i][i]-=int(Q)
  a,e,f=L[0];_,b,g=L[1];c=L[2][2]
  assert c>0 and c*b-g*g>0 and (a*b*c+2*e*f*g-a*g*g-b*f*f-c*e*e)>0,(cell,threshold)
  count+=1
 return {'sub':sub,'radius':1.5,'threshold':str(threshold),'inside':interior,'exterior':count}

def numerical():
 rng=np.random.default_rng(373);b=np.array(B); dyads=np.einsum('ni,nj->nij',b,b)/np.sum(b*b,axis=1)[:,None,None]; dirs=np.vstack([np.eye(3),[1.,1.,0.],[1.,1.,1.],rng.normal(size=(120,3))]); dirs=dirs/np.linalg.norm(dirs,axis=1)[:,None]; loc=[];gaps=[];slopes=[]
 for n in dirs:
  q=b@n
  for br in range(3):
   def eig(r):return np.linalg.eigvalsh(np.einsum('i,ijk->jk',2*(1-np.cos(r*q)),dyads))
   root=brentq(lambda r:eig(r)[br]-1.44,0,1.75); ev=eig(root);dprime=np.einsum('i,ijk->jk',2*q*np.sin(root*q),dyads)
   loc.append(root);gaps.append(ev[2]-ev[1]);slopes.append(np.linalg.eigvalsh(dprime)[0]);assert np.linalg.eigvalsh(dprime - (root/2)*np.eye(3))[0]>0
 return {'root_count':len(loc),'maximum_sample_radius':max(loc),'minimum_sample_upper_gap':min(gaps),'minimum_sample_radial_matrix_eig':min(slopes)}
if __name__=='__main__':
 s=F(12,25);u=F(9,4);f=lambda x:2*x-x*x/2;assert f(s)==F(528,625) and f(u)>f(s)
 print(json.dumps({'cert':cert(54,F(29,20)),'numerical':numerical(),'analytic_gap':str(f(s))},indent=2))
