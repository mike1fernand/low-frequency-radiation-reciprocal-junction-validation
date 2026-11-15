"""Compare the physical spring forces on a copied P1 state at the instant of release."""
import itertools
import json
import numpy as np
from math import asin,sqrt

def test():
    n=120
    k=2*asin(6/13)
    m=np.arange(n)-18
    q=.002*np.exp(-m*m/98)*np.cos(k*m)
    phase=2*np.pi*np.fft.fftfreq(n)
    omega=2.6*np.abs(np.sin(phase/2))
    temporal=-1j*np.sign(phase)*omega
    temporal[n//2]=0.
    v=np.fft.ifft(np.fft.fft(q)*temporal).real
    # Embed identical axial samples on the central row. Transverse boundaries are far from
    # the only nonzero row, so transverse wrap plays no role in the forces here.
    u=np.zeros((n,7,7,3))
    u[:,3,3,0]=q
    force=np.zeros_like(u)
    potential=0.
    for r in itertools.product((-1,0,1),repeat=3):
        norm=sum(x*x for x in r)
        if norm not in (1,2) or next(x for x in r if x)<0: continue
        axis=np.array(r)/sqrt(norm)
        slip=np.roll(u,tuple(-x for x in r),axis=(0,1,2))-u
        extension=slip@axis
        f=extension[...,None]*axis
        force+=f-np.roll(f,r,axis=(0,1,2))
        potential+=float(np.sum(extension**2))/2
    lap=np.roll(q,1)+np.roll(q,-1)-2*q
    assert np.allclose(force[:,3,3,0],lap-4*q,atol=2e-17,rtol=0)
    assert np.allclose(force[:,4,3,1],(np.roll(q,1)-np.roll(q,-1))/2,atol=2e-17,rtol=0)
    scalar_kinetic=(v@v)/2
    scalar_energy=scalar_kinetic+1.69*np.sum((np.roll(q,-1)-q)**2)/2
    assert abs(potential-(2*q@q+np.sum((np.roll(q,-1)-q)**2)/2))<1.e-16
    return {'scalar_initial_energy_periodic':scalar_energy,'copy_same_coordinates_energy':potential+scalar_kinetic,
            'copy_to_scalar_energy_ratio':(potential+scalar_kinetic)/scalar_energy,
            'central_site':18,'scalar_initial_axial_acceleration':1.69*lap[18],
            'ordinary_row_initial_axial_acceleration_with_other_nodes_at_rest':force[18,3,3,0],
            'off_strand_nodes_laterally_respond_maximum_initial_acceleration':float(np.max(np.abs(force[:,4,3,1])))}
if __name__=='__main__':print(json.dumps(test(),indent=2))
