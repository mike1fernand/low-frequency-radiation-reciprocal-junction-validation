"""Check explicit consequences of existing P1 laws; no scattering/trajectory run."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math
import numpy as np

Zq=Q(13,10);ksum=Q(8,3);omega_star=ksum/(2*Zq)
assert omega_star==Q(40,39)
CT=2*float(Zq)/(3*math.pi)
CL=float(Zq)/(3*math.pi*(math.sqrt(3)**3))
ratio=CT/CL
eta=CT/(CT+CL)
eta_exact=(108-6*math.sqrt(3))/107
assert abs(ratio-6*math.sqrt(3))<3e-15
assert abs(eta-eta_exact)<2e-16
curves=[]
for x in [Q(1,4),Q(1,2),Q(1),Q(2),Q(4)]:
 lam=1/(x*omega_star)
 F=(lam*ksum)**2/((lam*ksum)**2+(2*Zq)**2)
 master=1/(1+x*x)
 assert F==master
 curves.append({'x_omega_over_s_omega_star':str(x),'lambda_s_over_omega':str(lam),'exact_R_and_normalized_Tbulk':str(F),'numeric':float(F)})
cones=[]
for degrees in [10,20,30]:
 c=math.cos(math.radians(degrees))
 fT=(2-3*c+c**3)/4;fL=(1-c**3)/2
 share=CT*fT/(CT*fT+CL*fL)
 cones.append({'degrees':degrees,'limiting_qT_share':share})
# Constant-mass Fourier lower-band velocity mixtures obey the kinetic-energy
# corollary. This checks sampled instances; the projector norm proof is analytic.
bonds=np.array([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
unit=bonds/np.linalg.norm(bonds,axis=1)[:,None]
rng=np.random.default_rng(20261106);worst_ratio=0.;samples=0
for kappa in [.05,.1,.2,.4]:
 for _ in range(80):
  n=rng.normal(size=3);n/=np.linalg.norm(n);rho=rng.uniform(.5*kappa,kappa);k=rho*n
  D=np.einsum('b,bi,bj->ij',4*np.sin(bonds@k/2)**2,unit,unit)
  eig,V=np.linalg.eigh(D);assert eig[2]>eig[1]
  # A complex two-mode velocity, not a single eigenvector.
  velocity=V[:,:2]@(rng.normal(size=2)+1j*rng.normal(size=2))
  kinetic_fraction=abs(n@velocity)**2/np.vdot(velocity,velocity).real
  bound=kappa**4/16
  assert kinetic_fraction<=bound+3e-15
  worst_ratio=max(worst_ratio,float(kinetic_fraction/bound));samples+=1
scalings=[]
for alphaM,alphaK in [(4,1),(1,4),(4,4)]:
 omega_scale=math.sqrt(alphaK/alphaM)
 scalings.append({'alpha_M':alphaM,'alpha_K':alphaK,'frequency_scale':omega_scale,'fixed_displacement_power_scale':alphaK*omega_scale,'bare_only_effective_source_scale':1/alphaK})
report={'status':'PASS','scope':'Arithmetic/exact identities and sampled static projector corollary of already recorded P1/P3 laws. No new conversion efficiency, full scattering, proof of novelty or empirical validation.',
 'CT':CT,'CL':CL,'qT_to_qL_low_frequency_ratio':ratio,'qT_low_frequency_bulk_share':eta,'reference_omega_star_exact':str(omega_star),'normalized_crossover_exact_cases':curves,'existing_low_frequency_cone_shares':cones,
 'kinetic_corollary':{'samples':samples,'maximum_sampled_fraction_to_bound_ratio':worst_ratio,'at_kappa_0_2_maximum_fraction':.2**4/16,'at_kappa_0_2_maximum_percent':100*.2**4/16,'limit':'Global Fourier wave-normal projected kinetic energy of lower-band velocities with support |k|<=kappa, same constant-mass bulk. Not total elastic energy or incident conversion.'},'complete_mass_stiffness_scaling':scalings,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(report,indent=2,allow_nan=False))
