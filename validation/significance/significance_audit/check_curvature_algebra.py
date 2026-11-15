#!/usr/bin/env python3
"""Exact spherical-moment algebra for the proposed motif-purity curvature.

This does not prove the separate acoustic threshold expansion or solve a
scattering state. All calculations use standard-library rational arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

D = 8  # n_x,n_y,n_z,M,A,N,P,h; h is an origin shift/one-motif position.
ZERO = (0,)*D


def const(x):
    x = F(x)
    return {ZERO: x} if x else {}


def var(i):
    e = list(ZERO); e[i] = 1
    return {tuple(e): F(1)}


def add(*items):
    out = {}
    for item in items:
        for e, c in item.items():
            out[e] = out.get(e, F(0))+c
    return {e: c for e, c in out.items() if c}


def scale(a, s):
    return {e: c*F(s) for e, c in a.items() if c*F(s)}


def mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x+y for x,y in zip(ea,eb))
            out[e] = out.get(e,F(0))+ca*cb
    return {e: c for e,c in out.items() if c}


def power(a, j):
    out = const(1)
    for _ in range(j): out = mul(out,a)
    return out


def dot(a,b):
    return add(*(mul(x,y) for x,y in zip(a,b)))


def dfact(j):
    out = 1
    while j > 0: out *= j; j -= 2
    return out


def avg(a):
    """Normalized uniform S^2 moment, retaining M,A,N,P,h symbolically."""
    out = {}
    for e,c in a.items():
        if any(j % 2 for j in e[:3]): continue
        factor = F(dfact(e[0]-1)*dfact(e[1]-1)*dfact(e[2]-1),dfact(sum(e[:3])+1))
        ee = (0,0,0)+e[3:]
        out[ee] = out.get(ee,F(0))+c*factor
    return {e:c for e,c in out.items() if c}


def scalar(a):
    assert set(a).issubset({ZERO})
    return a.get(ZERO,F(0))


def show(a):
    names = ['nx','ny','nz','M','A','N','P','h']
    return {','.join(f'{n}^{p}' for n,p in zip(names,e) if p) or '1':str(c)
            for e,c in sorted(a.items())}


def main():
    nx,ny,nz,M,A,N,P,h = [var(i) for i in range(D)]
    n = [nx,ny,nz]
    transverse = add(power(ny,2),power(nz,2))
    first = [mul(M,nx),mul(N,ny),mul(N,nz)]
    second = [add(scale(mul(A,power(nx,2)),F(1,2)),scale(transverse,F(1,8))),
              mul(P,mul(nx,ny)),mul(P,mul(nx,nz))]
    total_correction = add(dot(first,first),scale(second[0],-2))
    L_correction = add(power(dot(n,first),2),scale(mul(nx,dot(n,second)),-2))
    T_correction = add(total_correction,scale(L_correction,-1))
    leading_T = F(2,3); leading_L = F(1,3)
    rT_geom = scale(avg(T_correction),1/leading_T)
    # q_L = x/sqrt(3), while q_T=x.
    rL_geom = scale(avg(L_correction),F(1,3)/leading_L)
    geometry_difference = add(rT_geom,scale(rL_geom,-1))
    expected_geometry = scale(add(scale(P,2),scale(mul(M,N),-2),
                                  scale(power(N,2),-1),const(F(-1,2))),F(1,3))
    assert geometry_difference == expected_geometry
    assert all(e[4] == 0 for e in geometry_difference)  # A cancels exactly.

    bonds = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),
             (1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
    A4 = [[{} for _ in range(3)] for _ in range(3)]
    for b in bonds:
        phase = add(*(scale(n[j],b[j]) for j in range(3)))
        fourth = power(phase,4)
        b2 = sum(x*x for x in b)
        for i in range(3):
            for j in range(3):
                A4[i][j] = add(A4[i][j],scale(fourth,F(-b[i]*b[j],12*b2)))
    axis_entry = A4[0][0]
    A4n = [dot(row,n) for row in A4]
    axial_A4n = mul(nx,A4n[0])
    axial2_h = mul(power(nx,2),dot(n,A4n))
    angular = [scalar(avg(p)) for p in [axis_entry,axial_A4n,axial2_h]]
    assert angular == [F(-3,20),F(-3,28),F(-3,28)]
    # Cluster formulas sum the whole transverse pair; no degenerate basis.
    T_material_num = add(scale(axis_entry,F(-5,2)),scale(axial_A4n,4),scale(axial2_h,F(-3,2)))
    L_material_num = add(scale(axial_A4n,F(1,3)),scale(axial2_h,F(-11,18)))
    rT_mat = scalar(avg(T_material_num))/leading_T
    rL_mat = scalar(avg(L_material_num))/leading_L
    assert (rT_mat,rL_mat) == (F(9,56),F(5,56))
    beta = add(geometry_difference,const(rT_mat-rL_mat))
    expected_beta = add(const(F(-2,21)),scale(add(P,scale(mul(M,N),-1)),F(2,3)),scale(power(N,2),F(-1,3)))
    assert beta == expected_beta

    translated_beta = add(const(F(-2,21)),
        scale(add(P,mul(h,N),scale(mul(add(M,h),N),-1)),F(2,3)),
        scale(power(N,2),F(-1,3)))
    assert translated_beta == beta
    # One motif: 2X=1, M=xi=h and P=xi N=hN.
    one_motif = add(const(F(-2,21)),
        scale(add(mul(h,N),scale(mul(h,N),-1)),F(2,3)),scale(power(N,2),F(-1,3)))
    assert one_motif == add(const(F(-2,21)),scale(power(N,2),F(-1,3)))

    result = {'status':'PASS_EXACT_CONDITIONAL_ALGEBRA',
        'method':'Sparse polynomial expansion and exact normalized spherical moments with Fraction arithmetic',
        'geometry_rT':show(rT_geom),'geometry_rL':show(rL_geom),
        'geometry_difference':show(geometry_difference),
        'quartic_symbol_angular_moments':{'A4_xx':str(angular[0]),'nx_A4n_x':str(angular[1]),'nx2_nA4n':str(angular[2])},
        'material_rT':str(rT_mat),'material_rL':str(rL_mat),'purity_beta':show(beta),
        'origin_translation_invariant':True,'axial_second_moment_cancels':True,
        'one_motif_beta':show(one_motif),'one_motif_sign':'strictly negative for every real N',
        'assumptions':['Actual full eight-link motif Cartesian force basis',
            'Real normalized static source 2 sum_i Xi=1',
            'Equal-axis/equal-face-diagonal H18 tangent bulk with cL/cT=sqrt(3)',
            'Whole-sphere branch power and required omega^3 spectral matrix threshold expansion',
            'No arbitrary finite-aperture detector or independent-link extension'],
        'not_proved':['Acoustic threshold regularity, reviewed separately',
            'Existence of a physically admissible multi-motif profile with positive beta',
            'Finite-frequency operating window or global monotonicity',
            'Historical novelty, electromagnetic or original-DST identification'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=Path(__file__).with_name('CURVATURE_ALGEBRA_CHECK.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
