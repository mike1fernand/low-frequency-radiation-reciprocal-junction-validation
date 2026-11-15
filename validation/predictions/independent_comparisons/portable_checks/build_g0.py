import ast
import numpy as np
from p1solver import *
def g0_blocks_from_scan(scanfile, sphfile, eps_list=(0.006,0.008,0.012,0.016,0.024), deg=2, Np=2048):
    st={ast.literal_eval(k):v for k,v in np.load(scanfile,allow_pickle=True).item().items()}
    sph=np.load(sphfile,allow_pickle=True).item(); nmax=sph['nmax']
    data=np.array([st[(e,Np)] for e in eps_list]); t=1j*np.array(eps_list)
    V=np.vander(t,deg+1,increasing=True)
    coef,*_=np.linalg.lstsq(V,data.reshape(len(eps_list),-1),rcond=None)
    c0=coef[0].reshape(data.shape[1:])
    ImT=sph['T'].real; ImL=sph['L'].real
    Gb=c0.real+1j*(ImT+ImL)[nmax-1:]          # n>=0 blocks: Re extrapolated, Im exact
    return Gb, ImT, ImL, nmax
