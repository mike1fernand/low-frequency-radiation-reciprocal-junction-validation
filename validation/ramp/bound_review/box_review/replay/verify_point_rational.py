"""Exact-rational residual certificate for the specified P1 ramp point.

NumPy only proposes an inverse. Every accepted inequality is verified using
Fraction arithmetic and rational intervals for sqrt(2) and sqrt(133).
No interval package or floating-point error assumption is used for acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib, json, time
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/REFERENCE_KERNEL.npz')
HASH='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==HASH
START=time.perf_counter()

class C:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):self.r=F(r);self.i=F(i)
    def __add__(self,b):
        if not isinstance(b,C):b=C(b)
        return C(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,b):return self+-asC(b)
    def __mul__(self,b):
        b=asC(b);return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):
        assert not isinstance(b,C);return C(self.r/F(b),self.i/F(b))
    def conj(self):return C(self.r,-self.i)
    def sq(self):return self.r*self.r+self.i*self.i
    def l1(self):return abs(self.r)+abs(self.i)
    def fl(self):return complex(float(self.r),float(self.i))

def asC(z):return z if isinstance(z,C) else C(z)
def zeros(m,n):return [[C() for _ in range(n)] for _ in range(m)]
def eye(n):return [[C(int(i==j)) for j in range(n)] for i in range(n)]
def minus(A,B):return [[a-b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def mm(A,B):
    bt=list(zip(*B));return [[sum((a*b for a,b in zip(row,col)),C()) for col in bt] for row in A]
def fro2(A):return sum(x.sq() for row in A for x in row)
def vec2(x):return sum(z.sq() for z in x)
def sqrt_bounds(n,digits=80):
    den=10**digits;lo=F(isqrt(n*den*den),den);return lo,lo+F(1,den)
S2L,S2U=sqrt_bounds(2);S13L,S13U=sqrt_bounds(133)
S2=(S2L+S2U)/2;S13=(S13L+S13U)/2
e2=(S2U-S2L)/2;e13=(S13U-S13L)/2
k=F(169,100);om=F(6,5);g=[F(j,7) for j in range(1,7)]
with np.load(SOURCE,allow_pickle=False) as f:
    gf=f['G0'][:12,:12].copy();qf=f['Im_qT'][:12,:12].copy()
    assert np.array_equal(qf,f['Im_branch_0'][:12,:12]+f['Im_branch_1'][:12,:12])
G=[[C(F(float(z.real)),F(float(z.imag))) for z in row] for row in gf]
Q=[[F(float(z)) for z in row] for row in qf]
assert all(Q[i][j]==Q[j][i] for i in range(12) for j in range(12))

# Exact LDL proves that the dyadic output metric is positive definite.
L=[[F(int(i==j)) for j in range(12)] for i in range(12)];piv=[]
for j in range(12):
    d=Q[j][j]-sum(L[j][h]**2*piv[h] for h in range(j))
    assert d>0;piv.append(d)
    for i in range(j+1,12):
        L[i][j]=(Q[i][j]-sum(L[i][h]*L[j][h]*piv[h] for h in range(j)))/d
assert all(sum(L[i][h]*piv[h]*L[j][h] for h in range(12))==Q[i][j] for i in range(12) for j in range(12))

def model(s2,s13):
    U=zeros(12,12);V=zeros(12,7);Kself=[F(0)]*7
    for j,gj in enumerate(g):
        U[2*j][2*j]=C(2*gj/3);U[2*j+1][2*j+1]=C(4*gj/3)
        V[2*j][j]=V[2*j][j+1]=C(-2*gj/3)
        V[2*j+1][j]=C(-2*s2*gj/3);V[2*j+1][j+1]=C(2*s2*gj/3)
        Kself[j]+=4*gj/3;Kself[j+1]+=4*gj/3
    A=zeros(19,19);D=zeros(19,12)
    for i in range(7):
        for j in range(7):
            h=C(om*om-2*k-Kself[i]) if i==j else C(k if abs(i-j)==1 else 0)
            if i==j and i in [0,6]:h+=C(F(97,100),3*s13/25)
            A[i][j]=h/k
    for i in range(7):
        for j in range(12):
            A[i][7+j]=sum(V[h][i].conj()*G[h][j] for h in range(12))
            D[i][j]=V[j][i].conj()/k
    for i in range(12):
        for j in range(7):A[7+i][j]=-V[i][j]/k
        for j in range(12):
            A[7+i][7+j]=C(int(i==j))+sum(U[i][h]*G[h][j] for h in range(12))
            D[7+i][j]=U[i][j]/k
    return A,D

Ac,Dc=model(S2,S13);A00,D00=model(F(0),F(0))
A2,D2=model(F(1),F(0));A2=minus(A2,A00);D2=minus(D2,D00)
A13,_=model(F(0),F(1));A13=minus(A13,A00)
deltaA=e2*sum(z.l1() for row in A2 for z in row)+e13*sum(z.l1() for row in A13 for z in row)
deltaD=e2*sum(z.l1() for row in D2 for z in row)
assert deltaA<F(1,10**60) and deltaD<F(1,10**60)

Rf=np.linalg.inv(np.array([[z.fl() for z in row] for row in Ac],complex))
R=[[C(F(float(z.real)),F(float(z.imag))) for z in row] for row in Rf]
E=minus(eye(19),mm(R,Ac));eps_c=F(1,10**13);eps=F(1,10**12)
assert fro2(E)<eps_c**2
assert fro2(R)<36
assert eps_c+6*deltaA<eps
mu=eps/(1-eps)
RD=mm(R,Dc);Hhat=RD[7:]
assert fro2(RD)<4
assert fro2(Hhat)<F(1893,2000)**2 # 0.9465
col=[R[i][0] for i in range(19)];v=col[7:]
assert vec2(col)<4 and vec2(v)<F(889,2000)**2 # 0.4445

def quadratic(v):
    out=sum((v[i].conj()*Q[i][j]*v[j] for i in range(12) for j in range(12)),C())
    assert out.i==0;return out.r

Jhat2=2*k*sum(quadratic(list(col)) for col in zip(*Hhat))
assert Jhat2<F(599,1000)**2
Fcap=F(2253,1000)
assert 2*k*sum(Q[i][i] for i in range(12))<Fcap**2
dc=mu*2+6*deltaD/(1-eps)
assert F(1893,2000)+dc<F(947,1000)
assert F(599,1000)+Fcap*dc<F(3,5)

gamlo=6*S13L/25;gamhi=6*S13U/25;c_cap=F(32,25)
assert gamhi/k<c_cap*c_cap
assert c_cap*(F(889,2000)+mu*2)<F(569,1000)
ahat_lo2=2*gamlo*quadratic(v);ahat_hi2=2*gamhi*quadratic(v)
assert ahat_lo2>F(42891,50000)**2 # .85782
assert ahat_hi2<F(17157,20000)**2 # .85785
ae=Fcap*c_cap*mu*2
assert F(42891,50000)-ae>F(4289,5000)
assert F(17157,20000)+ae<F(8579,10000)

def point_interval(eg,eq):
    z=F(569,1000)/(1-F(947,1000)*eg);d=F(3,5)*eg*z
    low=(F(4289,5000)-d)**2-2*eq*z*z
    high=(F(8579,10000)+d)**2+2*eq*z*z
    return {'eG':str(eg),'eQ':str(eq),'lower_rational':str(low),'upper_rational':str(high),
            'lower_display':float(low),'upper_display':float(high)},low,high
coarse,lo,hi=point_interval(F(1,50),F(1,200));assert lo>F(18,25)
tight,lo2,hi2=point_interval(F(1,250),F(1,500));assert lo2>F(183,250) and hi2<F(37,50)
receipt={
 'status':'PASS','date':'2026-11-08',
 'acceptance_arithmetic':'Exact Fraction comparisons; NumPy proposes inverse only. sqrt2/sqrt133 enclosed by integer-isqrt rational intervals of width1e-80.',
 'physical_parameters':{'g':['1/7','2/7','3/7','4/7','5/7','6/7'],'omega':'6/5','Kq':'169/100','a':1,'m_bulk':1,'m_spine':1,'K_axis':1,'K_diag':1,'lead_eta':0},
 'source_kernel_sha256':HASH,'source_kernel_path':str(SOURCE),
 'cached_metric':'Im_qT leading12; entrywise identical to saved floating sum Im_branch_0+Im_branch_1',
 'exact_Qc_LDL':'All12rationalpivotsstrictlypositive; exactreconstruction verified',
 'preconditioner_frobenius_cap':6,'residual_frobenius_cap':'1/1000000000000',
 'operator_caps':{'H_Frobenius':'947/1000','J_Frobenius':'3/5','eta_c_norm':'569/1000','a_c_lower':'4289/5000','a_c_upper':'8579/10000'},
 'coarse_conditional_interval':coarse,'tighter_conditional_interval':tight,
 'scope':'Rigorous candidate-specific propagation of assumed/authenticated kernel-error norm caps, fixed contact/real frequency. This script does not itself authenticate or rerun the contour integration. Not a bandwidth, packet or coupling-box certificate; not proof that exact value rounds to73.6%.',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.perf_counter()-START,
}
(HERE/'POINT_RATIONAL_CERTIFICATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_kernel_path','scope','script_sha256']},indent=2))
