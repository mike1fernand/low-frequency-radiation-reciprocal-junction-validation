"""Exact integer-coefficient checks of identities in the middle-branch proof.
No floating-point arithmetic or third-party packages are used.
"""
import json
from pathlib import Path

def add(*polys):
    result = {}
    for poly in polys:
        for power, coef in poly.items():
            result[power] = result.get(power, 0) + coef
    return {power: coef for power, coef in result.items() if coef}

def scale(poly, factor):
    return {power: coef*factor for power, coef in poly.items() if coef*factor}

def mul(a, b):
    result = {}
    for u, cu in a.items():
        for v, cv in b.items():
            power = tuple(x+y for x,y in zip(u,v))
            result[power] = result.get(power,0) + cu*cv
    return {power: coef for power, coef in result.items() if coef}

one = {(0,0,0):1}
y,z,delta = ({(1,0,0):1}, {(0,1,0):1}, {(0,0,1):1})
B,P = add(y,z),mul(y,z)
Eyy,Ezz = add(one,mul(delta,y),P), add(one,mul(delta,z),P)
sine_product_squared = mul(add(one,scale(mul(y,y),-1)),add(one,scale(mul(z,z),-1)))
determinant = add(mul(Eyy,Ezz),scale(sine_product_squared,-1))
factored = mul(add(B,delta),add(B,mul(delta,P)))
residual = add(determinant,scale(factored,-1))
assert not residual
# Here y plays c_i and z,delta play the other two cosines.
d_i = add(scale(one,2),mul(y,y),scale(mul(y,add(one,z,delta)),-1))
diagonal_D_over_2 = add(scale(one,3),scale(mul(y,add(one,z,delta)),-1))
rankone_diagonal = add(d_i,one,scale(mul(y,y),-1))
rankone_residual = add(diagonal_D_over_2,scale(rankone_diagonal,-1))
assert not rankone_residual
result = {'determinant_identity_residual':0,'rank_one_diagonal_identity_residual':0,'off_diagonal_identity':'D_ij/2=s_i*s_j by definition','method':'Exact integer polynomial coefficients; inequalities and interlacing checked separately in prose.'}
Path(__file__).with_name('exact_middle_branch_algebra.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
