"""Basis-invariant diagnostic residue rule for backward-error centroid ties.

This defines a quadrature representative at a tied sample, not a new branch
identity away from a tie. No eigenvalues or scalar cell integrals are changed.
"""
import numpy as np

def invariant_centroid_weights(vals,vec,F,khat):
    c=np.einsum('...am,...aA->...mA',vec.conj(),F,optimize=True)
    d=np.einsum('...am,...a->...m',vec.conj(),khat,optimize=True)
    W=np.einsum('...mA,...mB->...mAB',c.conj(),c,optimize=True)
    WL=W*abs(d)[...,None,None]**2
    tol=64*np.finfo(float).eps*np.maximum(1,np.max(abs(vals),axis=-1))
    tie01=abs(vals[...,1]-vals[...,0])<=tol
    tie12=abs(vals[...,2]-vals[...,1])<=tol
    counts={}
    for group,mask in [((0,1),tie01&~tie12),((1,2),tie12&~tie01),((0,1,2),tie01&tie12)]:
        rank=len(group);counts[str(group)]=int(mask.sum())
        if not np.any(mask):continue
        # F*P F / rank is invariant under any rotation of the tied eigenspace.
        wg=W[mask][:,group].sum(axis=1)/rank
        cg=c[mask][:,group];dg=d[mask][:,group]
        # (khat*P F)*(khat*P F) / rank, including cross terms within P.
        longitudinal=np.einsum('nm,nmA->nA',dg.conj(),cg,optimize=True)
        lg=np.einsum('nA,nB->nAB',longitudinal.conj(),longitudinal)/rank
        for branch in group:
            W[...,branch,:,:][mask]=wg
            WL[...,branch,:,:][mask]=lg
    return W,WL,counts
