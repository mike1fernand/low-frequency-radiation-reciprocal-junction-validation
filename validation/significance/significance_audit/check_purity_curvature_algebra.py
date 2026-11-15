#!/usr/bin/env python3
"""Exact rational sphere/local-corner checks; no external symbolic dependency.

This validates the displayed algebra on exact polynomial integrals and signed
force fixtures. It computes no Green kernel, static loaded state or scattering.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json

class Poly:
    def __init__(self,value=0):
        self.t=value if isinstance(value,dict) else ({(0,0,0):F(value)} if value else {})
        self.t={k:F(v) for k,v in self.t.items() if v}
    def __add__(self,other):
        other=aspoly(other);d=dict(self.t)
        for k,v in other.t.items():d[k]=d.get(k,F(0))+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.t.items()})
    def __sub__(self,other):return self+-aspoly(other)
    def __rsub__(self,other):return aspoly(other)+-self
    def __mul__(self,other):
        other=aspoly(other);d={}
        for k,v in self.t.items():
            for l,w in other.t.items():
                a=tuple(x+y for x,y in zip(k,l));d[a]=d.get(a,F(0))+v*w
        return Poly(d)
    __rmul__=__mul__
    def __truediv__(self,other):return self* (F(1)/other)
    def __pow__(self,n):
        out=Poly(1)
        for _ in range(n):out=out*self
        return out
    def __eq__(self,other):return self.t==aspoly(other).t
def aspoly(x):return x if isinstance(x,Poly) else Poly(x)
def df(n):
    out=1
    for i in range(n,0,-2):out*=i
    return out
def avg(p):
    total=F(0)
    for powers,c in aspoly(p).t.items():
        if any(i%2 for i in powers):continue
        total+=c*F(df(powers[0]-1)*df(powers[1]-1)*df(powers[2]-1),df(sum(powers)+1))
    return total
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def outer(a,b):return [[x*y for y in b] for x in a]
def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def mv(a,b):return [dot(row,b) for row in a]
def madd(a,b):return [[a[i][j]+b[i][j] for j in range(3)] for i in range(3)]
def scale(a,c):return [[c*v for v in row] for row in a]
def zeros():return [[Poly() for _ in range(3)] for _ in range(3)]

def main():
    n=[Poly({tuple(int(i==j) for i in range(3)):F(1)}) for j in range(3)]
    x,y,z=n;L=outer(n,n);I=[[Poly(int(i==j)) for j in range(3)] for i in range(3)]
    T=madd(I,scale(L,-1));D2=zeros();W=zeros()
    bonds=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
    for b in bonds:
        proj=[[F(b[i]*b[j],sum(v*v for v in b)) for j in range(3)] for i in range(3)]
        phase=dot(n,b);D2=madd(D2,scale(proj,phase**2));W=madd(W,scale(proj,-phase**4/12))
    assert D2==madd(scale(I,dot(n,n)),scale(L,2))
    explicit=zeros()
    for i in range(3):
        explicit[i][i]=-n[i]**4/4-sum(n[j]**4/12+n[i]**2*n[j]**2/2 for j in range(3) if i!=j)
        for j in range(3):
            if i!=j:explicit[i][j]=-n[i]*n[j]*(n[i]**2+n[j]**2)/3
    assert W==explicit
    wL=dot(n,mv(W,n));L2=scale(madd(mm(mm(T,W),L),mm(mm(L,W),T)),F(1,2))
    vals={'trace_W':avg(sum(W[i][i] for i in range(3))), 'wL':avg(wL),
          'L2_xx':avg(L2[0][0]),'wL_nx2':avg(wL*x*x),'TWT_xx':avg(mm(mm(T,W),T)[0][0])}
    assert vals=={'trace_W':-F(9,20),'wL':-F(9,28),'L2_xx':F(0),'wL_nx2':-F(3,28),'TWT_xx':-F(3,70)}
    rL=3*(vals['L2_xx']/3-F(5,18)*vals['wL_nx2']);rT=F(3,2)*(-vals['L2_xx']-F(5,2)*vals['TWT_xx'])
    assert (rL,rT)==(F(5,56),F(9,56))
    fixtures=[([F(1,2)],[F(1,2)],[F(j,3)]) for j in [-3,-1,0,2]]
    fixtures += [([F(1,2),F(3,2),F(9,2)],[F(1,3),-F(1,6),F(1,3)],[F(2,7),-F(1,5),F(3,8)]),
                 ([-F(1,2),F(3,2)],[-F(1,4),F(3,4)],[-F(4,9),F(5,11)])]
    rows=[]
    for xi,XX,rr in fixtures:
        assert sum(XX)==F(1,2)
        f0=[Poly() for _ in range(3)];f1=list(f0);f2=list(f0)
        for c,X,r in zip(xi,XX,rr):
            for sy in [-1,1]:
                for sz in [-1,1]:
                    pos=[c,F(sy,2),F(sz,2)];force=[X/2,r*sy/2,r*sz/2]
                    phase=dot(n,pos)
                    f0=[a+b for a,b in zip(f0,force)]
                    f1=[a+b*phase for a,b in zip(f1,force)]
                    f2=[a+b*phase**2 for a,b in zip(f2,force)]
        M=2*sum(c*X for c,X in zip(xi,XX));A=2*sum(c*c*X for c,X in zip(xi,XX));N=sum(rr);P=sum(c*r for c,r in zip(xi,rr))
        assert f0==[Poly(1),Poly(0),Poly(0)]
        assert f1==[M*x,N*y,N*z]
        assert f2==[A*x*x+(y*y+z*z)/4,2*P*x*y,2*P*x*z]
        ell=avg(dot(n,f1)**2-x*dot(n,f2));tt=avg(dot(f1,f1)-f2[0])-ell
        assert ell==(3*(M*M-A)+4*M*N-F(1,2)-4*P+8*N*N)/15
        assert tt==(2*(M*M-A)+2*N*N-4*M*N+4*P-2)/15
        beta=F(3,2)*tt-ell+rT-rL
        expected=-F(2,21)+2*(P-M*N)/3-N*N/3
        assert beta==expected
        d=F(7,3);assert expected==-F(2,21)+2*(P+d*N-(M+d)*N)/3-N*N/3
        if len(xi)==1:assert beta==-F(2,21)-N*N/3 and beta<0
        rows.append({'motifs':len(xi),'xi':list(map(str,xi)),'X':list(map(str,XX)),
                     'R_over_sqrt2':list(map(str,rr)),'beta':str(beta),'origin_invariance':True})
    assert F(1,14)!=-F(2,21) # omitting four-corner transverse extent changes result
    out={'status':'PASS','method':'standard-library exact rational polynomial sphere averages',
         'scope':'Actual nine-bond quartic/local-corner/angular algebra; no threshold theorem, Green kernel or loaded static source computed.',
         'algebraic_force_fixtures':rows,'fixtures_are_physically_realized_profiles':False,
         'sphere_averages':{k:str(v) for k,v in vals.items()},
         'relative_material_corrections':{'qL':str(rL),'qT':str(rT)},
         'beta':'-2/21+2*(P-M*N)/3-N*N/3',
         'source_hash':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    here=Path(__file__).resolve().parent
    (here/'PURITY_CURVATURE_ALGEBRA.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='algebraic_force_fixtures'},indent=2))
if __name__=='__main__':main()
