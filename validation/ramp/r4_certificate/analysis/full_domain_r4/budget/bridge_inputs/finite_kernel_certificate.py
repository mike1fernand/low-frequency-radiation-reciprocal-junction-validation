"""Outward-validated whole-box result for the specified cached bulk kernels.

This certifies finite arithmetic, not the exact infinite-lattice kernel error.
The optional error-budget conclusion remains conditional on external operator
enclosures for G=kstar*g0 and Qstar=kstar*Q.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time

import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--dependency-dir',type=Path)
    parser.add_argument('--precision',type=int,default=160)
    args=parser.parse_args()
    repo,out=args.repo_root.resolve(),args.output_dir.resolve()
    if out==repo or repo in out.parents:
        raise ValueError('Output must remain outside source repository')
    out.mkdir(parents=True,exist_ok=True)
    if args.dependency_dir:
        sys.path.insert(0,str(args.dependency_dir.resolve()))
    from flint import acb,acb_mat,arb,ctx,fmpq
    import flint
    ctx.prec=args.precision
    start=time.monotonic()
    failures=[]
    checked_pivots=0

    def require(condition,message):
        if not condition:
            raise ArithmeticError(message)

    def exact_float(x):
        return arb(fmpq(*float(x).as_integer_ratio()))

    def from_numpy(a):
        a=np.asarray(a)
        return acb_mat(a.shape[0],a.shape[1],[acb(exact_float(z.real),exact_float(z.imag)) for z in a.flat])

    def adj(a):
        return a.conjugate().transpose()

    def ident(n):
        a=acb_mat(n,n)
        for i in range(n): a[i,i]=1
        return a

    def approximation(a):
        return np.array([[complex(float(a[i,j].real.mid()),float(a[i,j].imag.mid())) for j in range(a.ncols())] for i in range(a.nrows())])

    def positive_definite(h):
        # Underlying H is Hermitian by construction. Real diagonal restriction
        # is therefore exact; intervals enclose its LDL* factorization.
        nonlocal checked_pivots
        n=h.nrows()
        l=[[acb(0) for j in range(n)] for i in range(n)]
        ds=[]
        for j in range(n):
            d=h[j,j].real
            for k in range(j):
                d-=(l[j][k]*l[j][k].conjugate()).real*ds[k]
            if not d>0:
                return False
            ds.append(d)
            checked_pivots+=1
            l[j][j]=1
            for i in range(j+1,n):
                z=h[i,j]
                for k in range(j):
                    z-=l[i][k]*l[j][k].conjugate()*ds[k]
                l[i][j]=z/d
        return True

    def gram_norm_upper(h,label):
        # Approximation only proposes a rational threshold. The acceptance test
        # below is an outward interval proof of threshold^2 I - H > 0.
        approx=approximation(h)
        scale=max(0.,float(np.linalg.eigvalsh((approx+approx.conj().T)/2)[-1]))**.5
        denominator=10**7
        numerator=max(1,math.ceil(scale*denominator)+1)
        for trial in range(5):
            upper=arb(fmpq(numerator,denominator))
            if positive_definite(ident(h.nrows())*(upper*upper)-h):
                return fmpq(numerator,denominator)
            failures.append(dict(stage='norm_threshold_proposal',label=label,numerator=numerator,denominator=denominator))
            numerator+=max(1,math.ceil(scale*denominator*10**(-5+trial)))
        raise ArithmeticError('Unable to verify norm bound for '+label)

    def op_upper(a,label):
        return gram_norm_upper(adj(a)*a,label)

    def vector_upper(x,label,metric=None):
        square=(adj(x)*(metric*x if metric is not None else x))[0,0].real
        require(square>=0,label+': norm square not enclosed nonnegative')
        value=square.sqrt()
        numerator=math.ceil(float(value.mid())*10**7)+1
        result=fmpq(max(numerator,1),10**7)
        require(value<arb(result),label+': proposed scalar upper bound failed')
        return result

    def fstr(q):
        return str(q)

    def maxq(values):
        return max(values)

    # Small rejection controls exercise the certified LDL test itself. These
    # are not production matrix substitutes and use no disabled assert checks.
    require(positive_definite(acb_mat([[2,acb(0,1)],[acb(0,-1),2]])),
        'LDL self-check rejected a known positive Hermitian matrix')
    require(not positive_definite(acb_mat([[1,2],[2,1]])),
        'LDL self-check accepted a known indefinite matrix')
    require(not positive_definite(acb_mat([[1,1],[1,1]])),
        'LDL self-check accepted a singular matrix as positive definite')

    pvpath=repo/'validation/frequency_domain/methods/surface/raw/pv_mean_M192.npy'
    shellpath=repo/'validation/frequency_domain/results/surface/surface_N128_matrices.npz'
    def real_array(value,shape,label):
        require(isinstance(value,np.ndarray) and value.shape==shape,
            label+': unexpected array shape')
        require(np.issubdtype(value.dtype,np.floating) and not np.iscomplexobj(value),
            label+': expected a real floating array')
        require(bool(np.all(np.isfinite(value))),label+': nonfinite input entry')
        return value

    pvraw=real_array(np.load(pvpath,allow_pickle=False),(20,20),'PV raw matrix')
    pv=real_array(pvraw[:12,:12],(12,12),'PV reduced matrix')
    with np.load(shellpath,allow_pickle=False) as z:
        imag=real_array(z['Im_total'],(12,12),'Im_total')
        q0=real_array(z['Im_branch_0'],(12,12),'Im_branch_0')
        q1=real_array(z['Im_branch_1'],(12,12),'Im_branch_1')
    for name,arr in [('PV reduced matrix',pv),('Im_total',imag),('Im_branch_0',q0),('Im_branch_1',q1)]:
        require(np.array_equal(arr,arr.conj().T),name+': input is not exactly Hermitian')
    source_hashes={str(path.relative_to(repo)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (pvpath,shellpath)}
    expected_hashes={
        'validation/frequency_domain/methods/surface/raw/pv_mean_M192.npy':'ce45000dbd154d556c1cb2a967cc5992149bb2f0b8a98153d654c855b8280e24',
        'validation/frequency_domain/results/surface/surface_N128_matrices.npz':'af6620513ca1e81def55d0ddcc1ddffdaa3d2f4a46d5c2cd3d0acf1f13e8c0bb',
    }
    require(source_hashes==expected_hashes,'Inputs do not match the specified archived cached-kernel identities')
    q=real_array(q0+q1,(12,12),'qT combined matrix')
    require(np.array_equal(pv,pv.T),'PV approximant is not exactly symmetric')
    require(np.array_equal(imag,imag.T),'Total spectral approximant is not exactly symmetric')
    require(np.array_equal(q,q.T),'qT spectral approximant is not exactly symmetric')
    g=from_numpy(pv+1j*imag)
    Q=from_numpy(q)
    require(positive_definite(Q),'Dyadic qT approximant not verified positive definite')
    nq,nc,n=7,12,19
    kstar=arb(fmpq(169,100))
    omega=arb(fmpq(6,5))
    coswave=1-omega*omega/(2*kstar)
    sinwave=(1-coswave*coswave).sqrt()
    gamma=2*kstar*sinwave
    sigma=-kstar*acb(coswave,sinwave)
    sqrt2=arb(2).sqrt()
    bare=acb_mat(nq,nq)
    for i in range(nq):
        bare[i,i]=omega*omega-2*kstar
        if i>0: bare[i,i-1]=kstar
        if i<nq-1: bare[i,i+1]=kstar
    bare[0,0]-=sigma
    bare[nq-1,nq-1]-=sigma
    b=acb_mat(n,1)
    b[0,0]=(gamma/kstar).sqrt()
    E=acb_mat(nc,n)
    for i in range(nc): E[i,nq+i]=1
    Hout=adj(E)*Q*E*(2*kstar)
    baseprofile=[fmpq(1,10),fmpq(1,5),fmpq(1,5),fmpq(1,5),fmpq(1,5),fmpq(1,10)]

    def system(signs):
        U,V,S=acb_mat(nc,nc),acb_mat(nc,nq),acb_mat(nq,nq)
        for s in range(6):
            gs=arb(baseprofile[s]*(1+fmpq(signs[s],20)))
            U[2*s,2*s]=2*gs/3
            U[2*s+1,2*s+1]=4*gs/3
            V[2*s,s]=V[2*s,s+1]=-2*gs/3
            V[2*s+1,s]=-2*gs*sqrt2/3
            V[2*s+1,s+1]=2*gs*sqrt2/3
            S[s,s]+=4*gs/3
            S[s+1,s+1]+=4*gs/3
        a=acb_mat(n,n)
        aa=(bare-S)/kstar
        ab=adj(V)*g
        ba=-V/kstar
        bb=ident(nc)+U*g
        for i in range(nq):
            for j in range(nq): a[i,j]=aa[i,j]
            for j in range(nc): a[i,nq+j]=ab[i,j]
        for i in range(nc):
            for j in range(nq): a[nq+i,j]=ba[i,j]
            for j in range(nc): a[nq+i,nq+j]=bb[i,j]
        D=acb_mat(n,nc)
        for i in range(nq):
            for j in range(nc): D[i,j]=V[j,i]/kstar
        for i in range(nc):
            for j in range(nc): D[nq+i,j]=U[i,j]/kstar
        return a,D,U

    A0,D0,U0=system([0]*6)
    R=A0.inv()
    x0=R*b
    W0=ident(nc)+U0*g
    RW=W0.inv()
    T0=(adj(x0)*Hout*x0)[0,0].real
    require(T0>0,'Nominal transmission not verified positive')
    thresholds=[]
    for count,signs in enumerate(itertools.product((-1,1),repeat=6),1):
        A,D,U=system(signs)
        B=R*(A-A0)
        y=B*x0
        RD=R*D
        vals=dict(r=op_upper(B,'r'),v=vector_upper(y,'v'),
            w=gram_norm_upper(adj(B)*Hout*B,'w'),d=vector_upper(y,'d',Hout),
            wE=op_upper(E*B,'wE'),dE=vector_upper(E*y,'dE'),
            h=op_upper(RD,'h'),hcore=op_upper(RW*U/kstar,'hcore'),
            beta=op_upper(RW*(U-U0)*g,'beta'))
        thresholds.append(dict(signs=list(signs),bounds={key:fstr(value) for key,value in vals.items()}))
        if count%16==0:
            print(json.dumps(dict(stage='finite_vertices',completed=count,elapsed_seconds=time.monotonic()-start)),flush=True)
    constants={key:maxq(fmpq(row['bounds'][key]) for row in thresholds) for key in thresholds[0]['bounds']}
    aa={key:arb(value) for key,value in constants.items()}
    r,v,w,d=(aa[key] for key in ('r','v','w','d'))
    require(r<1 and aa['beta']<1,'Uniform inverse conditions not verified')
    amplitude=T0.sqrt()-d-w*v/(1-r)
    require(amplitude>0,'No positive nominal lower amplitude')
    floor=amplitude*amplitude
    require(floor>arb(fmpq(27,250)),'Did not verify nominal floor > 0.108')
    Fnorm=arb(gram_norm_upper(Hout,'F_norm'))
    Z=arb(vector_upper(E*x0,'z0'))+aa['dE']+aa['wE']*v/(1-r)
    eG,eQ=arb(fmpq(1,50)),arb(fmpq(1,200))
    full_inverse=r+aa['h']*eG
    core_inverse=aa['beta']+aa['hcore']*eG
    require(full_inverse<1 and core_inverse<1,'Hypothetical error-budget inverse condition failed')
    state=aa['h']*eG*Z/(1-full_inverse)
    damp=amplitude-Fnorm*state
    require(damp>0,'Hypothetical error budget removed amplitude margin')
    conditional_floor=damp*damp-2*eQ*(Z+state)**2
    require(conditional_floor>arb(fmpq(1037,10000)),'Did not verify conditional error-budget floor > 0.1037')
    result=dict(status='finite_cached_kernel_certificate_passed',certificate_revision='basis_invariant_n192_sobol23_m256_v1',selection_sha256='d13f65f7fc1aa666fe29bca774584acd88644027204f0729106b2b0c78af14b6',original_reviewed_producer_sha256='1ccabcafc3165fcd790e667c49573075cfe0e6dd35f10694fdb85bff3bca3324',precision_bits=args.precision,python_flint=flint.__version__,
        elapsed_seconds=time.monotonic()-start,checked_positive_LDL_pivots=checked_pivots,
        claim='For the specified exact dyadic cached g0,Q matrices and exact rational mechanical parameters, all six independent complete-motif coupling factors in [19/20,21/20] have T_qT > 27/250.',
        exact_lattice_claim=False,
        input_interpretation='Complex and qT matrices use the exact dyadic values of the loaded NumPy arrays (including the saved-in-memory float sum of the two qT branch arrays). Physical decimal parameters 1.69,1.2,0.1,0.2 and0.05 are exact rationals; sqrt2 and scalar-lead radicals are outward-enclosed.',
        T0_enclosure=T0.str(40),finite_floor_enclosure=floor.str(40),finite_rational_floor='27/250',
        uniform_norm_upper_rationals={key:fstr(value) for key,value in constants.items()},
        hypothetical_kernel_error_conclusion=dict(
            hypotheses='IF exact lattice satisfies ||1.69*(g0_exact-g0_cached)||2<=1/50 and ||1.69*(Q_exact-Q_cached)||2<=1/200 with exact Q>=0, THEN T_qT>1037/10000 throughout same box.',
            errors_measured_or_certified=False,eG='1/50',eQ='1/200',floor_enclosure=conditional_floor.str(40),rational_floor='1037/10000',
            source_norm_upper=Z.str(40),full_inverse_upper=full_inverse.str(40),core_inverse_upper=core_inverse.str(40)),
        method='Arb/acb outward ball arithmetic; interval LDL* positivity proves all chosen operator norm bounds; no approximate flint solve/eigen modes.',
        ldl_controls='Known positive complex-Hermitian matrix accepted; known indefinite and singular matrices rejected using explicit runtime checks.',
        vertices=thresholds,failed_threshold_proposals=failures,
        source_sha256=source_hashes,
        input_validation='Raw PV is a finite real 20x20 array; its selected 12x12 block and each finite real 12x12 spectral array are exactly Hermitian. The source files must match the two specified archived SHA256 identities before arithmetic.',
        remaining=['No rigorous on-shell/PV integration enclosure for exact lattice.','No experimental or model-validity certification.'])
    (out/'finite_kernel_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('vertices',)},indent=2))


if __name__=='__main__':
    main()
