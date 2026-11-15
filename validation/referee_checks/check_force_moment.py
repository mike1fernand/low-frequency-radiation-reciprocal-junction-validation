"""Independent check of the first share correction using prescribed static force fixtures."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from itertools import product

bonds = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]])
bond_units = bonds/np.linalg.norm(bonds,axis=1)[:,None]
mu, weights=leggauss(12)
eta0=2*3**1.5/(1+2*3**1.5)
fixtures = [([0.5],[0.17]),([0.1,0.4],[0.3,-0.2]),([-0.1,0.2,0.4],[-0.05,0.15,-0.25])]

def force_fixture(xforce,rforce):
    positions, forces=[], []
    for s,(x,r) in enumerate(zip(xforce,rforce)):
        for sy,sz in product((-1,1),repeat=2):
            positions.append([s+0.5,sy/2,sz/2]);forces.append([x/2,sy*r/np.sqrt(8),sz*r/np.sqrt(8)])
    return np.array(positions),np.array(forces)

def integrate(omega):
    out = np.zeros((len(fixtures),3))
    sources = [force_fixture(*q) for q in fixtures]
    for cosine,weight in zip(mu,weights):
        for angle in (np.arange(24)+0.37)*(2*np.pi/24):
            sine = np.sqrt(1-cosine**2)
            normal = np.array([cosine,sine*np.cos(angle),sine*np.sin(angle)])
            phases = bonds@normal
            def dyad_sum(coefficients):
                return (bond_units.T*coefficients)@bond_units
            def eigenvalue(radius,branch):
                return np.linalg.eigvalsh(dyad_sum(4*np.sin(radius*phases/2)**2))[branch]
            for branch in range(3):
                radius = brentq(lambda radius:eigenvalue(radius,branch)-omega**2,0.4*omega,1.1*omega,xtol=1.e-14)
                _,basis = np.linalg.eigh(dyad_sum(4*np.sin(radius*phases/2)**2))
                polarization=basis[:,branch]
                radial_derivative=polarization@dyad_sum(2*phases*np.sin(radius*phases))@polarization
                for i,(sites,forces) in enumerate(sources):
                    fourier_force = np.exp(-1j*radius*(sites@normal))@forces
                    out[i,branch] += weight * radius**2 / radial_derivative * abs(polarization@fourier_force)**2
    shares=out[:,:2].sum(axis=1)/out.sum(axis=1)
    return (shares-eta0)/(omega**2*eta0*(1-eta0))

def predicted(xforce,rforce):
    s=np.arange(len(xforce))+0.5
    mu_x=2*s@xforce
    nu_r=np.sum(rforce)/np.sqrt(2)
    nu_xr=s@rforce/np.sqrt(2)
    return -2/21 + 2*(nu_xr - mu_x*nu_r)/3 - nu_r**2/3

if __name__=='__main__':
    targets=np.array([predicted(*f) for f in fixtures]);vals = [integrate(om) for om in (0.04,0.02)]
    print('static-moment predictions:',targets)
    print('relative-share coefficients at .04:',vals[0]);print('at .02:',vals[1]);print('Richardson errors:',(4*vals[1]-vals[0])/3-targets)
    np.testing.assert_allclose((4*vals[1]-vals[0])/3,targets,atol=3e-7,rtol=0)
