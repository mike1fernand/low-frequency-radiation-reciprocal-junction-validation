"""Diagnostic replacement for scalar linear-tetrahedron integral only.

No quadrature, Hamiltonian, source, spectral-weight or scattering change.
Uses a geometrically convergent simplex-moment series away from the pole;
else long-double divided differences with exact coincidences only, with
90-decimal fallback for clustered knots. Fixed near-knots are never merged.
"""
import math
import numpy as np
import mpmath as mp
COUNTS={'series':0,'long_double':0,'mp90':0}
# Some platforms alias long double to double; those use the high-precision path.
HAS_EXTENDED_PRECISION=np.finfo(np.longdouble).eps < np.finfo(float).eps

def _mp_reference(e,z):
    with mp.workdps(90):
        e=[mp.mpf(float(x)) for x in e];z=mp.mpc(z)
        def f(x,n):
            y=x-z
            return [y*y*mp.log(y)/2,y*mp.log(y)+y/2,mp.log(y)+mp.mpf('1.5'),1/y][n]
        dd=[f(x,0) for x in e]
        for n in range(1,4):
            dd=[(dd[i+1]-dd[i])/(e[i+n]-e[i]) if e[i+n]!=e[i] else f(e[i],n)/math.factorial(n) for i in range(4-n)]
        return complex(dd[0])

def stable_tetrahedron_resolvent(energies,z_energy,*,confluent_rtol=None):
    array=np.asarray(energies,dtype=float)
    if array.ndim<1 or array.shape[-1]!=4:raise ValueError('energies must have last dimension four')
    # confluent_rtol is accepted solely for compatibility; no approximate knots merge.
    shape=array.shape[:-1];e=np.sort(array,axis=-1).reshape(-1,4)
    if np.any(~np.isfinite(e)) or not np.isfinite(z_energy) or z_energy.imag<=0:raise ValueError('finite energies and retarded z required')
    midpoint=(e[:,0]+e[:,-1])/2;d=midpoint-z_energy
    q=(e[:,-1]-e[:,0])/(2*np.abs(d));far=q<=.25
    out=np.empty(len(e),dtype=np.complex128)
    if np.any(far):
        # I = d^-1 sum_n h_n(x)/[(n+1)(n+2)(n+3)], x_i=-(e_i-midpoint)/d.
        # Truncation after n=27: |tail| <= q^28/[6|d|(1-q)].
        x=-(e[far]-midpoint[far,None])/d[far,None]
        s1=x.sum(axis=1);s2=sum(x[:,i]*x[:,j] for i in range(4) for j in range(i+1,4))
        s3=sum(x[:,i]*x[:,j]*x[:,k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4));s4=x.prod(axis=1)
        hm3=np.zeros(len(x),complex);hm2=hm3.copy();hm1=hm3.copy();h0=np.ones(len(x),complex);total=h0/6
        for n in range(1,28):
            hn=s1*h0-s2*hm1+s3*hm2-s4*hm3
            total+=hn/((n+1)*(n+2)*(n+3));hm3,hm2,hm1,h0=hm2,hm1,h0,hn
        out[far]=total/d[far];COUNTS['series']+=int(far.sum())
    near=np.flatnonzero(~far)
    if len(near):
        en=e[near];span=en[:,-1]-en[:,0];gaps=np.diff(en,axis=1)
        # Cancellation risk is tied to relative knot separation, not |energy|.
        hard=(np.min(gaps,axis=1)<1e-4*span)|(np.min(abs(en-z_energy),axis=1)<1e-8*span)
        if not HAS_EXTENDED_PRECISION:hard[:]=True
        regular=near[~hard]
        if len(regular):
            ee=e[regular].astype(np.longdouble);scale=(ee[:,-1]-ee[:,0]);cen=np.longdouble(z_energy.real)
            xx=(ee-cen)/scale[:,None];zz=np.clongdouble(1j*z_energy.imag)/scale
            y=xx-zz[:,None]
            # Below the pole, log(y)=log(-y)-i*pi. The removed term is
            # quadratic in e, hence its third divided difference vanishes.
            # Removing it analytically preserves tiny positive off-shell Im.
            below=ee[:,-1] < np.longdouble(z_energy.real)
            logs=np.log(y)
            logs[below]=np.log(-y[below])
            dd=y*y*logs/2
            for order in range(1,4):
                dd=(dd[:,1:]-dd[:,:-1])/(xx[:,order:]-xx[:,:-order])
            out[regular]=np.asarray(dd[:,0]/scale,dtype=np.complex128);COUNTS['long_double']+=len(regular)
        for index in near[hard]:out[index]=_mp_reference(e[index],z_energy)
        COUNTS['mp90']+=int(hard.sum())
    if np.any(~np.isfinite(out)):raise FloatingPointError('nonfinite tetrahedron scalar')
    return out.reshape(shape)
