"""Replay exact algebra bookkeeping for the WR compression-coefficient rewrite.

No scattering, quadrature, symbolic-engine dependency or new scientific cases.
A fresh output path is required; frozen evidence is never overwritten.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json


def check():
    # WR Eq63 monomial, excluding the common4*u_c factor.
    powers={'A':Q(3,2),'Eb':Q(3,2),'Es':Q(-3,2),
            'rho_s':Q(1,2),'rho_b':Q(-1,2),'S':Q(1),'k':Q(2)}
    # Es=rho_s*c_t^2*A; k=omega*sqrt(rho_b/Eb).
    for old,sub in [('Es',{'rho_s':Q(1),'c_t':Q(2),'A':Q(1)}),
                    ('k',{'omega':Q(1),'rho_b':Q(1,2),'Eb':Q(-1,2)})]:
        power=powers.pop(old)
        for key,value in sub.items():
            powers[key]=powers.get(key,Q(0))+power*value
    powers={key:value for key,value in powers.items() if value}
    expected={'Eb':Q(1,2),'rho_s':Q(-1),'rho_b':Q(1,2),
              'S':Q(1),'c_t':Q(-3),'omega':Q(2)}
    assert powers==expected,(powers,expected)
    # Same three rational examples as the original executed check.
    examples=[]
    for alpha in [Q(1,4),Q(1,3),Q(2,5)]:
        sigma=(1-2*alpha)/(2*(1-alpha))
        assert 2*(1+sigma)==(3-4*alpha)/(1-alpha)
        examples.append({'alpha':str(alpha),'sigma':str(sigma)})
    return {'status':'PASS',
      'method':'Exact rational exponent cancellation and3rational material identity examples; no scattering/quadrature run.',
      'identity':'WR Eq63 compression coefficient =4*Z_b*u_c(alpha)*omega^2/(rho_s*c_t^3), Z_b=S*sqrt(E_b*rho_b).',
      'resulting_monomial':{key:str(value) for key,value in powers.items()},
      'material_identity_examples':examples,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'scope':'Replay of the original arithmetic only. SymPy was not installed or used. The three material examples are sample checks; the coefficient cancellation is exact given the displayed modulus relation.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; choose a fresh filename')
    result=check()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write('\n')
    print(json.dumps({'status':result['status'],'output':str(args.output)},indent=2))


if __name__=='__main__':
    main()
