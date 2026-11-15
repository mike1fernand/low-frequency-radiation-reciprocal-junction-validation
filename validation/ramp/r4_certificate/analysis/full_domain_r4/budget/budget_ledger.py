"""Exact matrix error gate for certified *complete-domain* moment rectangles.

Inputs must separately be proved to enclose the complete integrals. Numerical
eigensolvers merely propose rational bounds, accepted only by exact LDL*. This
module does not promote upstream coverage/analyticity claims on its own.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math

import matrix_budget as mb


def add(a, b): return a[0]+b[0], a[1]+b[1]
def sub(a, b): return a[0]-b[0], a[1]-b[1]
def mul(a, b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def conj(a): return a[0], -a[1]


def exact_positive_definite(h):
    """LDL* with exact complex rational entries, no approximate acceptance."""
    n = len(h)
    if any(len(row) != n for row in h):
        raise ValueError('Square matrix required')
    for i in range(n):
        for j in range(n):
            if h[i][j] != conj(h[j][i]):
                raise ValueError('Hermitian matrix required')
    l = [[(F(0), F(0)) for _ in range(n)] for _ in range(n)]
    ds = []
    for j in range(n):
        d = h[j][j][0]
        for k in range(j):
            d -= (l[j][k][0]**2+l[j][k][1]**2)*ds[k]
        if d <= 0:
            return False
        ds.append(d)
        l[j][j] = F(1), F(0)
        for i in range(j+1, n):
            z = h[i][j]
            for k in range(j):
                x = mul(l[i][k], conj(l[j][k]))
                z = sub(z, (x[0]*ds[k], x[1]*ds[k]))
            l[i][j] = z[0]/d, z[1]/d
    return True


def gram(d):
    n = len(d)
    return [[tuple(sum((mul(conj(d[k][i]), d[k][j])[c] for k in range(n)), F(0))
                   for c in range(2)) for j in range(n)] for i in range(n)]


def exact_operator_upper(d, denominator=10**10):
    """Propose using SVD; accept only rational u with u²I-D*D>0."""
    import numpy as np
    if type(denominator) is not int or denominator < 1:
        raise ValueError('Positive integer proposal denominator required')
    n = len(d)
    if n == 0 or any(len(row) != n for row in d):
        raise ValueError('Nonempty square matrix required')
    if any(len(z) != 2 for row in d for z in row):
        raise ValueError('Exact real/imaginary entry pairs required')
    d = [[tuple(mb.rational(x) for x in z) for z in row] for row in d]
    if not any(x for row in d for z in row for x in z):
        return F(0), 0
    g = gram(d)
    approximate = np.array([[complex(float(z[0]), float(z[1])) for z in row] for row in d])
    proposed = float(np.linalg.norm(approximate, 2))
    if not math.isfinite(proposed):
        raise ValueError('Nonfinite bound proposal')
    numerator = max(1, math.ceil(proposed*denominator)+1)
    for trial in range(12):
        upper = F(numerator, denominator)
        test = [[sub((upper*upper if i == j else F(0), F(0)), g[i][j])
                 for j in range(n)] for i in range(n)]
        if exact_positive_definite(test):
            return upper, trial+1
        numerator += max(1, numerator//10**max(1, 10-trial))
    raise ArithmeticError('Unable to certify proposed center operator norm')


def decompose(intervals):
    z = mb.parse_intervals(intervals)
    centers = [tuple((lo+hi)/2 for lo, hi in x) for x in z]
    radii = [tuple((hi-lo)/2 for lo, hi in x) for x in z]
    return centers, radii


def radius_envelope(radii):
    return [[F(0) if mb.slot(i,j)[0] is None else
             mb.sqrt_up(sum(x*x for x in radii[mb.slot(i,j)[0]]))
             for j in range(12)] for i in range(12)]


def assess(intervals, target, kind):
    direct = mb.bound(intervals, target, kind)
    centers, radii = decompose(intervals)
    target = mb.target_matrix(target)
    difference = []
    for i in range(12):
        row = []
        for j in range(12):
            index, sign = mb.slot(i,j)
            value = (F(0), F(0)) if index is None else tuple(sign*x for x in centers[index])
            row.append(sub(value, target[i][j]))
        difference.append(row)
    center_bound, attempts = exact_operator_upper(difference)
    radius_bound = mb.norm_bounds(radius_envelope(radii))
    separated = mb.KSTAR*(center_bound+radius_bound['operator_upper'])
    best = min(direct['scaled_operator_upper'], separated)
    return {
        'kind': kind, 'direct_envelope': direct,
        'moment_centers': centers, 'moment_component_radii': radii,
        'center_operator_upper_exact_ldl': center_bound,
        'center_ldl_proposal_attempts': attempts,
        'scaled_center_operator_upper': mb.KSTAR*center_bound,
        'uncertainty_envelope': radius_bound,
        'scaled_uncertainty_operator_upper': mb.KSTAR*radius_bound['operator_upper'],
        'scaled_separated_operator_upper': separated,
        'scaled_operator_upper': best,
        'threshold': mb.THRESHOLDS[kind],
        'arithmetic_budget_pass': best <= mb.THRESHOLDS[kind],
        'exact_lattice_floor_certified': False,
        'claim': 'Arithmetic only; complete-domain enclosures and PSD premise require separate validation.'
    }


def rank_regions(records):
    """Rank reducible uncertainty, not physical contribution or center shift.

    Each record is one independent integral enclosure with all physical copies
    counted in positive integer multiplicity. Duplicated/correlated copies can
    be summed with that multiplicity without assuming independence.
    """
    output = []
    for record in records:
        weight = record.get('multiplicity', 1)
        if type(weight) is not int or weight < 1:
            raise ValueError('Physical multiplicity must be a positive integer')
        _, radii = decompose(record['intervals'])
        norms = mb.norm_bounds(radius_envelope(radii))
        contributions = []
        rows = mb.multiplicities()['rows']
        for k, rr in enumerate(radii):
            modulus = mb.sqrt_up(rr[0]*rr[0]+rr[1]*rr[1])
            contributions.append(weight*mb.KSTAR*max(row[k] for row in rows)*modulus)
        output.append({'id': record['id'], 'multiplicity': weight,
                       'scaled_uncertainty_operator_upper': weight*mb.KSTAR*norms['operator_upper'],
                       'per_moment_scaled_max_row_contribution': contributions,
                       'priority_note': 'Upper-bound reduction potential, not a true-error estimate.'})
    output.sort(key=lambda r: r['scaled_uncertainty_operator_upper'], reverse=True)
    return output


def recommend_allocations(center_bound, kind, relative_costs=None, reserve=F(1,10)):
    """Sufficient nonuniform global modulus radii from exact row constraints.

    Costs are positive rational heuristics; they affect allocation only. All
    candidate radii are then rescaled to satisfy every matrix row. A reserve
    fraction leaves room for subsequent midpoint drift; it is not its bound.
    """
    center_bound = mb.rational(center_bound)
    reserve = mb.rational(reserve)
    if not 0 <= reserve < 1:
        raise ValueError('Reserve fraction must be in [0,1)')
    remaining = mb.THRESHOLDS[kind]/mb.KSTAR-center_bound
    if remaining <= 0:
        return {'allocation_exists': False, 'reason': 'Current certified center upper bound consumes sufficient budget.'}
    costs = [F(1)]*17 if relative_costs is None else [mb.rational(x) for x in relative_costs]
    if len(costs) != 17 or any(x <= 0 for x in costs):
        raise ValueError('Seventeen positive rational relative costs required')
    rows = mb.multiplicities()['rows']
    # Larger cost permits larger errors, but rare moments can receive more.
    raw = [costs[k]/max(row[k] for row in rows) for k in range(17)]
    denominator = max(sum(row[k]*raw[k] for k in range(17)) for row in rows)
    radii = [(1-reserve)*remaining*x/denominator for x in raw]
    exact_row_use = max(sum(row[k]*radii[k] for k in range(17)) for row in rows)
    if exact_row_use > (1-reserve)*remaining:
        raise ArithmeticError('Allocation row constraint failed')
    return {'allocation_exists': True, 'global_moment_modulus_radii': radii,
            'uniform_physical_chart_density': [r/3 for r in radii],
            'unused_operator_reserve': remaining-exact_row_use,
            'remaining_operator_budget_after_current_center_upper': remaining,
            'reserve_fraction': reserve,
            'meaning': 'Sufficient conditional allocation about the supplied certified center bound; final intervals must still be checked.'}


def allocate_regions(regions, moment_modulus_radii, kind):
    """Turn a global moment allocation into exact region-specific tolerances.

    Positive area*cost weights permit larger tolerance on expensive regions.
    Physical copies enter the denominator, so reuse is never free. These
    conditional allocations are not a substitute for achieved-ball gates.
    """
    if kind not in ('g','Q'):
        raise ValueError('Expected g or Q')
    radii=[mb.rational(x) for x in moment_modulus_radii]
    if len(radii)!=17 or any(x<=0 for x in radii):
        raise ValueError('Seventeen positive global radii required')
    parsed=[]
    for r in regions:
        area=mb.rational(r['area']);cost=mb.rational(r.get('relative_cost','1'))
        copies=r['multiplicity']
        if area<=0 or cost<=0 or type(copies) is not int or copies<1:
            raise ValueError('Positive areas, costs and integral multiplicities required')
        parsed.append((r['id'],area,cost,copies))
    if not parsed or len({r[0] for r in parsed})!=len(parsed):
        raise ValueError('Nonempty region list with unique identities required')
    denominator=sum((area*cost*copies for _,area,cost,copies in parsed),F(0))
    divisor=mb.sqrt_up(F(2)) if kind=='g' else F(1)
    output=[]
    for name,area,cost,copies in parsed:
        local=[x*area*cost/denominator for x in radii]
        output.append({'id':name,'multiplicity':copies,'moment_modulus_radii':local,
                       'component_density_tolerances':[x/area/divisor for x in local]})
    for k in range(17):
        if sum((x['multiplicity']*x['moment_modulus_radii'][k] for x in output),F(0))!=radii[k]:
            raise ArithmeticError('Physical region allocation conservation failed')
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True,
                   help='JSON with complete_moments mapping g and Q to17 complex rectangles')
    p.add_argument('--targets', type=Path, default=Path(__file__).with_name('TARGETS.json'))
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    data = json.loads(args.input.read_text())
    target = json.loads(args.targets.read_text())
    results = {kind: assess(data['complete_moments'][kind], target['targets'][kind], kind)
               for kind in ('g', 'Q')}
    output = {'schema': 'paper1.exact-matrix-budget-r4.v1',
              'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'target_sha256': hashlib.sha256(args.targets.read_bytes()).hexdigest(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'results': results,
              'both_arithmetic_budgets_pass': all(x['arithmetic_budget_pass'] for x in results.values()),
              'exact_lattice_floor_certified': False}
    args.output.write_text(json.dumps(mb.json_exact(output), indent=2)+'\n')
    print(json.dumps({kind: {'scaled_upper': float(x['scaled_operator_upper']),
                            'pass': x['arithmetic_budget_pass']} for kind, x in results.items()}))


if __name__ == '__main__':
    main()
