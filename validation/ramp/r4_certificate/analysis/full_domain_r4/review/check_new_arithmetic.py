"""Independent R4 checks: exact 2x2 Gram tests and compact dyadic arithmetic."""
import copy
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'budget'), str(ROOT/'integrator')]
import budget_ledger as bl
import compact_controller as cc


def main():
    rng = random.Random(41723)
    checks = []
    # Independent 2x2 Hermitian criterion: positive first diagonal and determinant.
    for index in range(30):
        d = [[(F(rng.randrange(-11,12),7), F(rng.randrange(-11,12),5))
              for _ in range(2)] for _ in range(2)]
        upper, _ = bl.exact_operator_upper(d)
        a = sum(z[0]*z[0]+z[1]*z[1] for z in (d[0][0], d[1][0]))
        b = sum(z[0]*z[0]+z[1]*z[1] for z in (d[0][1], d[1][1]))
        cre = sum(d[k][0][0]*d[k][1][0]+d[k][0][1]*d[k][1][1] for k in range(2))
        cim = sum(d[k][0][0]*d[k][1][1]-d[k][0][1]*d[k][1][0] for k in range(2))
        assert upper*upper > a
        assert (upper*upper-a)*(upper*upper-b) > cre*cre+cim*cim
        # Subtract enough to make the top eigenvalue fail, retaining a positive candidate.
        small = upper/F(2)
        h = [[(small*small-a,F(0)),(-cre,-cim)],
             [(-cre,cim),(small*small-b,F(0))]]
        explicit = small*small > a and (small*small-a)*(small*small-b)>cre*cre+cim*cim
        assert bl.exact_positive_definite(h) == explicit
    checks.append('30 complex 2x2 upper bounds and LDL outcomes checked by independent exact determinant criterion')

    for _ in range(200):
        a = cc.normalize_endpoint(rng.randrange(-10**12,10**12),rng.randrange(-220,220))
        b = cc.normalize_endpoint(rng.randrange(-10**12,10**12),rng.randrange(-220,220))
        def frac(end):
            m,e=int(end[0]),int(end[1])
            return F(m)*(F(2)**e)
        assert frac(cc.add_endpoint(a,b)) == frac(a)+frac(b)
    checks.append('200 mixed-sign unequal-exponent compact sums checked independently with Fraction')

    complete = ROOT/'integrator'/'BODY_SAMPLE.json'
    if complete.exists():
        state=json.loads(complete.read_text())
        cc.audit_state(state)
        if not state['pending_ids']:
            independent = [[F(0) for _ in range(4)] for _ in range(17)]
            for stage in state['stages']:
                for m,row in enumerate(stage['increment_endpoints']):
                    for k,(mantissa, exponent) in enumerate(row):
                        independent[m][k] += F(int(mantissa))*F(2)**int(exponent)
            observed = [[F(int(m))*F(2)**int(e) for m,e in row]
                        for row in state['aggregate_endpoints']]
            assert independent == observed
            checks.append('complete physical checkpoint stage aggregate independently summed as exact Fractions')
            # Recompute outer digest so rejection tests geometry, not checksum alone.
            bad = copy.deepcopy(state)
            bad['accepted_ids'].pop()
            bad['state_sha256']=cc.digest({k:v for k,v in bad.items() if k!='state_sha256'})
            try: cc.audit_state(bad)
            except ValueError: pass
            else: raise AssertionError('missing covered node passed after outer checksum recomputation')
            checks.append('missing complete physical leaf rejected after outer checksum recomputation')

    out={'status':'passed','checks':checks,'count':len(checks),
         'limits':'No full-domain integral, source primitive, bridge, or old production simulation was rerun.'}
    dest=Path(__file__).with_name('NEW_ARITHMETIC_CHECKS.json')
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__': main()
