"""P1 contact covariance + diagnostic response selection rule. Not a P1 radiation solve."""
from pathlib import Path
import itertools, json, hashlib
import numpy as np

corners=list(itertools.product((-1,1),repeat=2))
dim=14
vectors=[]
for endpoint,sigma in [(0,1),(1,-1)]:
    for c,(sy,sz) in enumerate(corners):
        v=np.zeros(dim);v[endpoint]=-sigma/np.sqrt(3)
        v[2+3*c:5+3*c]=np.array([sigma,sy,sz])/np.sqrt(3)
        vectors.append(v)
links=np.array([np.outer(v,v) for v in vectors])
group=[]; permutations=[]
for swap in (False,True):
    for signs in itertools.product((-1,1),repeat=2):
        S=np.diag(signs)@(np.array([[0,1],[1,0]]) if swap else np.eye(2))
        R=np.zeros((dim,dim));R[:2,:2]=np.eye(2)
        M=np.eye(3);M[1:,1:]=S
        perm=[]
        for c,yz in enumerate(corners):
            target=corners.index(tuple(int(x) for x in S@yz))
            R[2+3*target:5+3*target,2+3*c:5+3*c]=M
            perm.append(target)
        group.append(R);permutations.append(perm)
Pi=sum(group)/len(group);Q=np.eye(dim)-Pi
covariance=max(np.linalg.norm(R@links[e*4+c]@R.T-links[e*4+perm[c]])
               for R,perm in zip(group,permutations) for e in (0,1) for c in range(4))
assert covariance<1e-14
assert np.linalg.norm(Pi@Pi-Pi)<1e-14
# An intentionally finite diagnostic system: same physical contact Hessian,
# isotropic supports and passive damping. This is NOT the infinite P1 model.
Kbase=np.zeros((dim,dim));line=np.array([[2.5,-2.],[-2.,2.5]])
Kbase[:2,:2]=line;Kbase[2:,2:]=3*np.eye(12)
g0=np.full(8,.2);K0=Kbase+np.einsum('i,ijk->jk',g0,links)
Gamma=.2*np.eye(dim);omega=.7
A=K0-omega**2*np.eye(dim)-1j*Gamma
G=np.linalg.inv(A);b=np.zeros(dim);b[0]=1;x0=G@b
W=np.diag([0.,0.]+[1.]*12)
Wpatch=np.zeros((dim,dim));Wpatch[2:5,2:5]=np.eye(3)
assert max(np.linalg.norm(R@A@R.T-A) for R in group)<1e-14
assert max(np.linalg.norm(R@W@R.T-W) for R in group)<1e-14
# Both endpoint orbits separately sum to zero.
delta=.2*np.array([1.,-1.,0.,0.,0.,1.,0.,-1.])
dK=np.einsum('i,ijk->jk',delta,links)
mean_dK=sum(R@dK@R.T for R in group)/len(group)
dx=-G@dK@x0
power=lambda x,W:float(np.vdot(x,W@x).real)
derivative=lambda dx,W:float(2*np.vdot(x0,W@dx).real)
whole_slope=derivative(dx,W);patch_slope=derivative(dx,Wpatch)
assert np.linalg.norm(mean_dK)<1e-14
assert abs(whole_slope)<1e-13
assert abs(patch_slope)>1e-6
assert np.linalg.norm(Q@x0)<1e-14
assert np.linalg.norm(Pi@dx)<1e-14
rows=[]
for eps in [.1,.05,.025,.0125]:
    xp=np.linalg.solve(A+eps*dK,b);xm=np.linalg.solve(A-eps*dK,b)
    rows.append({'epsilon':eps,'whole_change_plus':power(xp,W)-power(x0,W),
                 'whole_central_slope':(power(xp,W)-power(xm,W))/(2*eps),
                 'whole_even_change_over_epsilon_squared':(.5*(power(xp,W)+power(xm,W))-power(x0,W))/eps**2,
                 'new_irrep_power_over_epsilon_squared':power(Q@xp,W)/eps**2,
                 'patch_central_slope':(power(xp,Wpatch)-power(xm,Wpatch))/(2*eps)})
assert abs(rows[-1]['whole_central_slope']) < 2e-9
assert all(3.9 < abs(a['whole_central_slope']/b['whole_central_slope']) < 4.1 for a,b in zip(rows,rows[1:]))
assert abs(rows[-1]['patch_central_slope']-patch_slope)<1e-7
# Total sum zero is NOT orbit-wise zero.
misleading=.2*np.array([1.]*4+[-1.]*4)
dK_m=np.einsum('i,ijk->jk',misleading,links)
mean_m=sum(R@dK_m@R.T for R in group)/len(group)
control_slope=derivative(-G@dK_m@x0,W)
assert misleading.sum()==0 and np.linalg.norm(mean_m)>1e-3
assert abs(control_slope)>1e-5
report={'status':'PASS','group':'C4v signed transverse permutations, order8',
        'one_motif_endpoint_orbits':[4,4],
        'six_motif_link_parameters':48,'six_motif_invariant_orbit_means':12,
        'maximum_contact_hessian_covariance_residual':covariance,
        'pure_breaking_Reynolds_hessian_norm':float(np.linalg.norm(mean_dK)),
        'diagnostic_whole_power_first_derivative':whole_slope,
        'diagnostic_single_corner_power_first_derivative':patch_slope,
        'total8sum_zero_but_orbit_mean_nonzero_control_derivative':control_slope,
        'refinement_rows':rows,
        'scope':'Actual P1 one-motif central-link Hessians checked exactly up to floating arithmetic. Response numbers use a finite supported/damped14DOF diagnostic, NOT P1 infinite-lattice qT radiation or the claimed robustness floor. Group proof supplies selection rule under its explicit regularity assumptions.',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
