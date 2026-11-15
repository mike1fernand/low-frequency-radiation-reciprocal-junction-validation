"""Compare complete certified enclosures with independent R3 ordinary quadrature.

Containment is a numerical sanity check only. The certified matrix gate does
not use these sampled values or their convergence agreement.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ledger',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    ledger=json.loads(args.ledger.read_text())
    if any(F(x['fallback_parameter_area']) for x in ledger['results'].values()):
        raise ValueError('This check requires complete refined enclosures, not fallback.')
    if not all(x['all_new_jobs_complete'] for x in ledger['results'].values()):
        raise ValueError('Every scheduled job must be complete.')
    source=args.baseline/'matrix/numerical/gl32_128.json'
    manifest=json.loads((args.baseline/'MANIFEST.json').read_text())
    if sha(source)!=manifest['files']['matrix/numerical/gl32_128.json']:
        raise ValueError('R3 diagnostic source differs from sealed manifest.')
    diagnostic=json.loads(source.read_text())
    assert diagnostic['angular_order']==32 and diagnostic['radial_order']==128
    assert diagnostic['energy']=='36/25'
    details={}
    for kind in ('g','Q'):
        intervals=ledger['complete_moments'][kind]
        reference=diagnostic['comparison'][kind]['moments']
        assert len(intervals)==len(reference)==17
        comparisons=[]
        for m,(row,value) in enumerate(zip(intervals,reference)):
            assert len(row)==len(value)==2
            for part,(bounds,number) in enumerate(zip(row,value)):
                lo,hi=map(F,bounds)
                # Recover the exact dyadic value represented by the saved float.
                center=F(*float(number).as_integer_ratio())
                if not lo<=center<=hi:
                    raise ArithmeticError(f'R3 {kind} moment{m} component{part} lies outside full certified enclosure')
                radius=(hi-lo)/2
                offset=abs((lo+hi)/2-center)
                comparisons.append({'moment_index':m,'component':'real' if part==0 else 'imag',
                                    'contained':True,'center_offset_diagnostic':float(offset),
                                    'offset_over_radius_diagnostic':float(offset/radius) if radius else None})
        details[kind]={'component_checks':len(comparisons),'comparisons':comparisons}
    out={'status':'passed','ledger_sha256':sha(args.ledger),
         'R3_ordinary_gl32_128_sha256':sha(source),'results':details,
         'claim':'All17 complex g and17 real Q moment enclosures contain the retained ordinary32/128 quadrature values, compared componentwise.',
         'limits':'Diagnostic agreement only. No new ordinary quadrature was run. It supplies no rigorous error bound and is not used to establish the independent certified operator gate.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','g_components':34,'Q_components':34,'claim_limit':out['limits']}))


if __name__=='__main__':main()
