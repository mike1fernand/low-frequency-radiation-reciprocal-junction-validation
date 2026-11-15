"""Tensor GL3/GL5 using exact-radical nodes and rigorous Cauchy tails.

The callback must be holomorphic wherever its complex ball evaluation is
finite. In this application that follows from trigonometric entire functions
and a finite interval enclosure of the full inverse of D-E.
"""
from itertools import product
from flint import arb, acb, ctx
from functools import lru_cache


class NeighborhoodFailure(ArithmeticError):
    def __init__(self, axis, reason='no finite analytic neighborhood'):
        self.axis = axis
        self.reason = reason
        super().__init__(reason + ' on axis ' + str(axis))


def rule(order):
    if order == 3:
        r = (arb(3)/5).sqrt()
        return [-r, arb(0), r], [arb(5)/18, arb(4)/9, arb(5)/18]
    if order == 5:
        s = (arb(10)/7).sqrt()
        inner, outer = (5-2*s).sqrt()/3, (5+2*s).sqrt()/3
        wi, wo = (322+13*arb(70).sqrt())/1800, (322-13*arb(70).sqrt())/1800
        return [-outer, -inner, arb(0), inner, outer], [wo, wi, arb(64)/225, wi, wo]
    raise ValueError('only analytically specified GL orders 3 and 5 are supported')


def finite_values(callback, xs):
    vals = list(callback(*xs))
    if not vals or not all(v.is_finite() for v in vals):
        raise ArithmeticError('nonfinite analytic-neighborhood enclosure')
    return vals

@lru_cache(maxsize=16)
def even_residuals(order, precision_bits):
    """Exact Gauss formula enclosed at the calling precision, cached safely."""
    if precision_bits!=ctx.prec:raise ValueError('precision cache mismatch')
    nodes,weights=rule(order)
    return tuple(abs(arb(1)/(k+1)-sum(w*x**k for x,w in zip(nodes,weights))).upper()
                 for k in range(2*order,2*order+21,2))

def tail_factor(order,r):
    """Cauchy factor retaining actual even Gauss moment residuals.

    Odd monomials integrate and quadrature to zero exactly. Terms through
    degree2n-1 agree exactly. After eleven enclosed even residuals, the
    remaining even coefficients have |integral-Q|<=2 by positivity and
    normalization. This formula requires no unproved error-sign property.
    """
    if not (0<=r and r<1):raise ValueError('Cauchy ratio must lie in [0,1)')
    result=arb(0);power=r**(2*order)
    for residual in even_residuals(order,ctx.prec):
        result+=residual*power;power*=r*r
    return result+2*power/(1-r*r)


def gauss_box(callback, bounds, orders=(5,5,5), radius_multipliers=(2,4,8,16), density_tolerance=None):
    """Enclose a complete 3D box; a failed bound returns a split coordinate.

    For order n the normalized quadrature and integral agree through degree
    2n-1. On a containing complex disk with |f|<=M and r=h/rho<1, the error
    is <=M*tail_factor(n,r). Tensor telescoping sums three such bounds.
    The square evaluated here contains the disk. Other coordinates remain
    their entire real intervals, establishing a uniform tensor bound.
    """
    if len(bounds)!=3 or len(orders)!=3:
        raise ValueError('three bounds and three quadrature orders required')
    rules=[rule(n) for n in orders]
    center=[(a+b)/2 for a,b in bounds]
    half=[(b-a)/2 for a,b in bounds]
    if not all(h>0 for h in half):
        raise ValueError('bounds must have positive widths')
    if density_tolerance is not None and not density_tolerance>0:
        raise ValueError('density tolerance must be positive')
    volume=8*half[0]*half[1]*half[2]
    real=[acb(a.union(b)) for a,b in bounds]
    errors=[]; chosen=[]; calls=0; rejected=0
    for axis in range(3):
        choices=[]
        for mult in radius_multipliers:
            if mult<=1: raise ValueError('radius multipliers must exceed one')
            radius=half[axis]*mult
            xs=list(real)
            xs[axis]=acb(center[axis]+arb(0,radius.abs_upper()),arb(0,radius.abs_upper()))
            try:
                calls+=1; vals=finite_values(callback,xs)
                r=half[axis].abs_upper()/radius.lower()
                if not r<1: raise ArithmeticError('noncontracting disk')
                factor=tail_factor(orders[axis],r)
                err=[abs(v).upper()*factor for v in vals]
                choices.append((float(max(e.upper() for e in err)),err,mult))
            except (ArithmeticError,ZeroDivisionError,ValueError):
                rejected+=1
        if not choices: raise NeighborhoodFailure(axis)
        _, err, mult=min(choices,key=lambda v:v[0])
        errors.append(err);chosen.append(mult)
    axis_max=[max(e.upper() for e in ee) for ee in errors]
    if density_tolerance is not None and not sum(axis_max)<density_tolerance:
        raise NeighborhoodFailure(max(range(3),key=lambda j:float(axis_max[j])), 'Cauchy remainder exceeds tolerance')
    accum=[acb(0) for _ in errors[0]]
    for indices in product(*[range(n) for n in orders]):
        xs=[acb(center[j]+half[j]*rules[j][0][indices[j]]) for j in range(3)]
        calls+=1;vals=finite_values(callback,xs)
        if len(vals)!=len(accum):raise ValueError('callback size changed')
        w=arb(1)
        for j in range(3): w*=rules[j][1][indices[j]]
        for m,v in enumerate(vals):accum[m]+=w*v
    remainders=[volume*sum(errors[j][m] for j in range(3)) for m in range(len(accum))]
    result=[volume*x+acb(arb(0,e.upper()),arb(0,e.upper())) for x,e in zip(accum,remainders)]
    return result, {'calls':calls,'rejected_neighborhoods':rejected,
                    'chosen_radius_multipliers':chosen,'orders':list(orders),
                    'maximum_remainder_diagnostic':float(max(e.upper() for e in remainders)),
                    'maximum_density_remainder_diagnostic':float(max(e.upper() for e in remainders)/volume)}
