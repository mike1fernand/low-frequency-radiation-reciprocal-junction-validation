"""Exact-rational auxiliary-vertex certificate for all six ±5% errors.

The vertices bound convex auxiliary norms and a linear scalar. No assertion that
transmission extrema occur at vertices is used. There is no adaptive search.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, itertools, json, time

ROOT=Path(__file__).resolve().parent
POINT=ROOT/'verify_point_rational.py'
# Reuse its exact arithmetic/model construction, stopping before output writing.
text=POINT.read_text();assert text.count('receipt={')==1
exec(compile(text.split('receipt={')[0],str(POINT),'exec'),globals())
START_BOX=time.perf_counter()
g0=g.copy();A0=Ac;D0=Dc;R0=R
vhat=[R0[i][0] for i in range(19)]
vhat_eta=vhat[7:]

def mv(A,v):return [sum((a*b for a,b in zip(row,v)),C()) for row in A]
def weighted2(A):return 2*k*sum(quadratic(list(col)) for col in zip(*A))
def cross(x,y):
    return sum((x[i].conj()*Q[i][j]*y[j] for i in range(12) for j in range(12)),C())

# Fraction caps have intentionally visible margins over all exact candidate
# preconditioner quantities. Residual enclosures then cover the true inverse.
caps={
 'r':F(164,1000),'w':F(68,1000),'v':F(1648,10000),
 'ell':F(173,10000),'w_eta':F(649,10000),'d_eta':F(504,10000),
 'h0':F(994,1000),'j0':F(6289,10000),'RD':F(1744,1000),
}
rows=[]
for signs in itertools.product([-1,1],repeat=6):
    g=[gj*(1+F(sign,20)) for gj,sign in zip(g0,signs)]
    Av,Dv=model(S2,S13)
    gap=minus(Av,A0);Bh=mm(R0,gap);Bh_eta=Bh[7:]
    RDh=mm(R0,Dv);RDh_eta=RDh[7:]
    Bhv=mv(Bh,vhat);Bhv_eta=Bhv[7:]

    # Both contact models are linear in sqrt2; the identical lead contribution
    # cancels from the perturbation. Build the exact coefficient of the residual
    # radical interval rather than assuming machine construction is exact.
    A_v00,D_v00=model(F(0),F(0))
    A_v2,D_v2=model(F(1),F(0))
    A_v2=minus(A_v2,A_v00);D_v2=minus(D_v2,D_v00)
    delta_gap=e2*sum(z.l1() for row in minus(A_v2,A2) for z in row)
    delta_d=e2*sum(z.l1() for row in D_v2 for z in row)
    assert delta_gap<F(1,10**60) and delta_d<F(1,10**60)

    # Exact dyadic preconditioner caps.
    rb2=fro2(Bh);wb2=weighted2(Bh_eta)
    vb2=(gamhi/k)*vec2(Bhv)
    web2=fro2(Bh_eta);deb2=(gamhi/k)*vec2(Bhv_eta)
    h02=fro2(RDh_eta);j02=weighted2(RDh_eta);rd2=fro2(RDh)
    ellhat=4*gamhi*abs(cross(vhat_eta,Bhv_eta).r)
    assert rb2<F(16395,100000)**2
    assert wb2<F(6795,100000)**2
    assert vb2<F(1647,10000)**2
    assert web2<F(6485,100000)**2
    assert deb2<F(503,10000)**2
    assert h02<F(9939,10000)**2
    assert j02<F(6288,10000)**2
    assert rd2<F(17438,10000)**2
    assert ellhat<F(17275,1000000)

    # Transfer from proposed inverse/products to the exact physical matrices.
    db=mu*F(16395,100000)+6*delta_gap/(1-eps)
    dd=mu*F(17438,10000)+6*delta_d/(1-eps)
    yerr=c_cap*mu*2
    ynorm=c_cap*(2+mu*2)
    byerr=db*ynorm+F(16395,100000)*yerr
    outputerr=Fcap*byerr
    assert F(16395,100000)+db<caps['r']
    assert F(6795,100000)+Fcap*db<caps['w']
    assert F(1647,10000)+byerr<caps['v']
    assert F(6485,100000)+db<caps['w_eta']
    assert F(503,10000)+byerr<caps['d_eta']
    assert F(9939,10000)+dd<caps['h0']
    assert F(6288,10000)+Fcap*dd<caps['j0']
    assert F(17438,10000)+dd<caps['RD']
    # The scalar is 2Re[(F y0)^dagger F B y0]. Its interval correction
    # bounds both factors and their product, with no directional assumption.
    ell_error=2*(ae*Fcap*F(1647,10000)+F(858,1000)*outputerr+ae*outputerr)
    assert ellhat+ell_error<caps['ell']
    rows.append({'signs':list(signs),'g':[str(x) for x in g],
                 'exact_preconditioner_norm_squared_displays':{
                     'B_frobenius_squared':float(rb2),'FB_frobenius_squared':float(wb2),
                     'B_y0_norm_squared_upper':float(vb2),'P_B_frobenius_squared':float(web2),
                     'P_B_y0_squared_upper':float(deb2),'P_RD_frobenius_squared':float(h02),
                     'F_RD_frobenius_squared':float(j02),'RD_frobenius_squared':float(rd2),
                     'linear_energy_scalar_absolute_upper':float(ellhat)}})

r=caps['r'];w=caps['w'];v=caps['v'];ell=caps['ell']
rem=w*v/(1-r)
Tcache_low=F(4289,5000)**2-ell-2*F(8579,10000)*rem
Lcache=F(8339,10000)
assert Lcache*Lcache<Tcache_low
Zcache=F(569,1000)+caps['d_eta']+caps['w_eta']*v/(1-r)
Huniform=caps['h0']+caps['w_eta']*caps['RD']/(1-r)
Juniform=caps['j0']+w*caps['RD']/(1-r)
Zcap=F(633,1000);Hcap=F(1130,1000);Jcap=F(772,1000)
assert Zcache<Zcap and Huniform<Hcap and Juniform<Jcap

def bound(eg,eq):
    assert Hcap*eg<1
    zz=Zcap/(1-Hcap*eg);dd=Jcap*eg*zz
    assert Lcache>dd
    low=(Lcache-dd)**2-2*eq*zz*zz
    return {'eG':str(eg),'eQ':str(eq),'lower_rational':str(low),'lower_display':float(low)},low
coarse,coarse_lo=bound(F(1,50),F(1,200))
tight,tight_lo=bound(F(1,250),F(1,500))
assert coarse_lo>F(674,1000)
assert tight_lo>F(69,100)
receipt={
 'status':'PASS','date':'2026-11-08','domain':'g_j=(j/7)(1+delta_j), all six delta_j independently in[-1/20,1/20]; eachmotif8linksvariedtogether;omega6/5,remainingparametersfixed',
 'vertex_count':64,'uses_transmission_vertex_extrema_assumption':False,
 'acceptance_arithmetic':'Exact Fraction inequalities with the independently verified point inverse residual; interval radicals width1e-80.',
 'source_kernel_sha256':HASH,'point_verifier_sha256':hashlib.sha256(POINT.read_bytes()).hexdigest(),
 'auxiliary_norm_caps':{x:str(y) for x,y in caps.items()},
 'finite_matrix_energy_floor_rational':str(Tcache_low),'finite_matrix_energy_floor_display':float(Tcache_low),
 'finite_matrix_output_amplitude_lower':str(Lcache),
 'uniform_source_caps':{'eta_c_norm':str(Zcap),'Hnorm':str(Hcap),'Jnorm':str(Jcap)},
 'coarse_kernel_box_lower':coarse,'tight_kernel_box_lower':tight,
 'conclusion':'Under authenticatedscaledkernelerrorcaps eG<=1/250,eQ<=1/500, T_qT>69/100 throughout the entire six-motif ±5% coupling box.',
 'scope':'Candidate-specific continuous coupling-box lowerbound with realphysicallead and completeqTmetric. Not a bandwidth or packet result, no claim that exactnominalpoint rounds to73.6%. Kernelpremiseauthentication belongs to the linked contour-replay chain.',
 'vertex_details':rows,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'seconds':time.perf_counter()-START_BOX,
}
(ROOT/'BOX_RATIONAL_CERTIFICATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='vertex_details'},indent=2))
