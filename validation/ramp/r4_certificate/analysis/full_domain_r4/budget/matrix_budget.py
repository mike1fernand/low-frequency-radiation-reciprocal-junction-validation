"""Exact rational moment assembly and *sufficient* operator-error bounds.

Input intervals must separately be proved to contain the complete integrals.
This module checks arithmetic, not that upstream mathematical premise.
"""
from fractions import Fraction as F
from math import isqrt

NAMES = tuple([f'a{i}' for i in range(6)] + [f'c{i}' for i in range(6)] +
              [f'b{i}' for i in range(1, 6)])
KSTAR = F(169, 100)
THRESHOLDS = {'g': F(1, 50), 'Q': F(1, 200)}


def rational(x):
    # JSON floating endpoints would silently turn a claimed decimal enclosure
    # into binary values. Demand exact integers or explicit rational strings.
    if isinstance(x, bool) or not isinstance(x, (str, int, F)):
        raise TypeError('Exact endpoints require int, Fraction, or rational string')
    return F(x)


def sqrt_up(q, bits=180):
    """Smallest dyadic on a 2**-bits grid whose square is >= nonnegative q."""
    q = rational(q)
    if q < 0 or bits < 1:
        raise ValueError('Invalid square-root bound')
    numerator_root, denominator_root = isqrt(q.numerator), isqrt(q.denominator)
    if numerator_root*numerator_root == q.numerator and denominator_root*denominator_root == q.denominator:
        return F(numerator_root,denominator_root)
    scale = 1 << bits
    n = isqrt((q.numerator * scale * scale) // q.denominator)
    if n*n*q.denominator < q.numerator*scale*scale:
        n += 1
    return F(n, scale)


def slot(i, j):
    """Return (moment index, sign), with None for the exact b0 zero."""
    if not (0 <= i < 12 and 0 <= j < 12):
        raise ValueError('Matrix index out of range')
    s, alpha = divmod(i, 2)
    t, beta = divmod(j, 2)
    d = s-t
    if alpha == beta:
        return ((0 if alpha == 0 else 6)+abs(d), 1)
    if d == 0:
        return (None, 0)
    sign = (1 if d > 0 else -1) * (1 if alpha == 0 else -1)
    return (11+abs(d), sign)


def assemble(moments):
    if len(moments) != 17:
        raise ValueError('Exactly 17 ordered moments required')
    return [[0 if slot(i,j)[0] is None else
             slot(i,j)[1]*moments[slot(i,j)[0]]
             for j in range(12)] for i in range(12)]


def multiplicities():
    counts = [0]*17
    rows = [[0]*17 for _ in range(12)]
    for i in range(12):
        for j in range(12):
            k, _ = slot(i,j)
            if k is not None:
                counts[k] += 1
                rows[i][k] += 1
    return {'entries': counts, 'rows': rows}


def parse_intervals(values):
    if len(values) != 17:
        raise ValueError('Exactly 17 ordered complex rectangles required')
    out=[]
    for z in values:
        if len(z) != 2 or any(len(pair) != 2 for pair in z):
            raise ValueError('Each complex interval requires real and imaginary endpoint pairs')
        pair = tuple(tuple(rational(v) for v in p) for p in z)
        if any(lo > hi for lo,hi in pair):
            raise ValueError('Reversed interval endpoints')
        out.append(pair)
    return out


def target_matrix(values):
    if len(values) != 12 or any(len(row) != 12 for row in values):
        raise ValueError('Target must have shape12x12')
    if any(len(z) != 2 for row in values for z in row):
        raise ValueError('Target entry must have exact real and imaginary parts')
    return [[tuple(rational(x) for x in z) for z in row] for row in values]


def envelope(intervals, target):
    """Entrywise modulus bound, compared with *every actual target entry*."""
    intervals = parse_intervals(intervals)
    target = target_matrix(target)
    out=[]
    for i in range(12):
        row=[]
        for j in range(12):
            k, sign = slot(i,j)
            z = ((F(0),F(0)), (F(0),F(0))) if k is None else intervals[k]
            if sign < 0:
                z = tuple((-hi,-lo) for lo,hi in z)
            magnitudes = [max(abs(lo-v),abs(hi-v))
                          for (lo,hi),v in zip(z,target[i][j])]
            row.append(sqrt_up(sum(v*v for v in magnitudes)))
        out.append(row)
    return out


def norm_bounds(E):
    """For |Delta_ij|<=E_ij, bound ||Delta||2 without symmetry assumptions."""
    if len(E) != 12 or any(len(row) != 12 for row in E):
        raise ValueError('Envelope must have shape12x12')
    E=[[rational(v) for v in row] for row in E]
    if any(x < 0 for row in E for x in row):
        raise ValueError('Envelope cannot contain negative entries')
    row=max(sum(row) for row in E)
    col=max(sum(E[i][j] for i in range(12)) for j in range(12))
    induced=sqrt_up(row*col)
    fro=sqrt_up(sum(v*v for r in E for v in r))
    return {'row_sum':row, 'column_sum':col,
            'sqrt_row_column':induced, 'frobenius':fro,
            'operator_upper':min(induced,fro)}


def bound(intervals, target, kind):
    if kind not in THRESHOLDS:
        raise ValueError('Expected g or Q')
    result=norm_bounds(envelope(intervals,target))
    result['scaled_operator_upper']=KSTAR*result['operator_upper']
    result['threshold']=THRESHOLDS[kind]
    result['arithmetic_budget_pass']=result['scaled_operator_upper']<=result['threshold']
    result['certifies_global_integral']=False
    return result


def json_exact(x):
    if isinstance(x,F):
        return str(x)
    if isinstance(x,dict):
        return {k:json_exact(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [json_exact(v) for v in x]
    return x
