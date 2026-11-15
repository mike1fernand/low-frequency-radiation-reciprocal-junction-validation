"""Arb source-projected resolvent moments at exact E=36/25.

Moment order: a_0..a_5, c_0..c_5, b_1..b_5.  All values include
positive-octant normalization 1/pi**3 and the cube-ray Jacobian t**2.
For complex t, the source products are analytic continuations of their
real-t values.  Conjugating t or the source on the contour is incorrect.
No angular permutation is applied to the physical source x direction.
"""
from flint import arb, acb, acb_mat, ctx


def finite_matrix(M):
    return all(M[i,j].is_finite() for i in range(3) for j in range(3))


def exact_midpoint_matrix(M):
    """Each real and imaginary midpoint is an exact dyadic Arb number."""
    out=acb_mat([[acb(M[i,j].real.mid(),M[i,j].imag.mid()) for j in range(3)] for i in range(3)])
    if not all(out[i,j].real.is_exact() and out[i,j].imag.is_exact() for i in range(3) for j in range(3)):
        raise ArithmeticError('midpoint preconditioner is not exact')
    return out


def preconditioned_inverse(M):
    """Certify M^-1=(P0*M)^-1*P0 for an exactly represented dyadic P0.

    The midpoint inverse is used only to choose P0. It is never assumed
    exact. Finiteness of the full interval inverse of P0*M proves the claim.
    """
    P0=exact_midpoint_matrix(exact_midpoint_matrix(M).inv())
    result=(P0*M).inv()*P0
    if not finite_matrix(result):raise ArithmeticError('nonfinite preconditioned inverse')
    return result


def certified_inverse(M):
    """Both branches are rigorous; the width threshold affects cost only.

    Point-quadrature balls have radii near rounding precision. Retain their
    cheap direct inverses. For larger neighborhood balls, try both methods
    and keep the finite enclosure with the smaller maximum entry modulus.
    """
    candidates=[]
    try:
        direct=M.inv()
        if finite_matrix(direct):candidates.append(direct)
    except (ArithmeticError,ZeroDivisionError):pass
    threshold=arb(2)**(-(ctx.prec//2))
    wide=any(M[i,j].real.rad()>threshold or M[i,j].imag.rad()>threshold for i in range(3) for j in range(3))
    if candidates and not wide:return candidates[0]
    try:candidates.append(preconditioned_inverse(M))
    except (ArithmeticError,ZeroDivisionError):pass
    if not candidates:raise ArithmeticError('no finite inverse enclosure')
    return min(candidates,key=lambda R:float(max(abs(R[i,j]).upper() for i in range(3) for j in range(3))))


def coordinate_values(t, u, v, chart=0):
    t,u,v=acb(t),acb(u),acb(v)
    if chart == 0: return (t,t*u,t*v)
    if chart == 1: return (t*u,t,t*v)
    if chart == 2: return (t*u,t*v,t)
    raise ValueError('chart must be 0, 1, or 2')


def symbol_matrix(t,u,v,chart=0):
    k=coordinate_values(t,u,v,chart)
    c=[x.cos() for x in k]; s=[x.sin() for x in k]
    return acb_mat([[6-2*c[i]*(1+c[(i+1)%3]+c[(i+2)%3])
                     if i==j else 2*s[i]*s[j]
                     for j in range(3)] for i in range(3)])


def shifted_matrix(t,u,v,chart=0):
    D=symbol_matrix(t,u,v,chart)
    E=arb(36)/25
    for i in range(3): D[i,i]-=E
    return D


def moment_values(t,u,v,chart=0):
    """Return 17 complex balls enclosing every input represented by t,u,v."""
    t=acb(t); kx,ky,kz=coordinate_values(t,u,v,chart)
    R=certified_inverse(shifted_matrix(t,u,v,chart))
    root2=arb(2).sqrt()
    A=2*(ky/2).cos()*(kz/2).cos()
    B=root2*(ky/2).sin()*(kz/2).cos()
    C=root2*(ky/2).cos()*(kz/2).sin()
    aa=A*A*R[0,0]
    cc=B*B*R[1,1]+2*B*C*R[1,2]+C*C*R[2,2]
    bb=A*(B*R[0,1]+C*R[0,2])
    factor=t*t/(arb.pi()**3)
    out=[factor*aa*(d*kx).cos() for d in range(6)]
    out += [factor*cc*(d*kx).cos() for d in range(6)]
    out += [factor*bb*(d*kx).sin() for d in range(1,6)]
    if not all(x.is_finite() for x in out):
        raise ArithmeticError('nonfinite resolvent moment enclosure')
    return out


def ball_serialization(x):
    """Exact binary endpoints plus explicitly diagnostic decimal summaries."""
    def endpoint(a):
        m,e=a.man_exp()
        return {'mantissa':str(m),'exponent2':int(e)}
    def realball(a):
        return {'lower':endpoint(a.lower()), 'upper':endpoint(a.upper()),
                'midpoint_diagnostic':float(a.mid()), 'radius_diagnostic':float(a.rad())}
    if isinstance(x,acb): return {'real':realball(x.real),'imag':realball(x.imag)}
    return realball(x)
