"""Independent symbolic cross-check of the short reference-shell proof."""
from fractions import Fraction
from sympy import Matrix, diag, simplify, symbols, zeros

cosines = symbols('c:3')
sines = symbols('v:3')
bonds = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
symbol = zeros(3)
for bond in bonds:
    axes = [i for i in range(3) if bond[i]]
    if len(axes) == 1:
        cosine = cosines[axes[0]]
    else:
        i, j = axes
        cosine = cosines[i]*cosines[j] - bond[i]*bond[j]*sines[i]*sines[j]
    b = Matrix(bond)
    symbol += 2*(1 - cosine)*b*b.T/b.dot(b)

h = [1 - c for c in cosines]
H = sum(h)
v = Matrix(sines)
claimed = 2*diag(*[H*(1-t) + 2*t*t for t in h]) + 2*v*v.T
identities = {v[i]**2: 1 - cosines[i]**2 for i in range(3)}
assert (symbol-claimed).expand().subs(identities).applyfunc(simplify) == zeros(3)
assert Fraction(9,4) * Fraction(13,16)**2 == Fraction(1521,1024) > Fraction(36,25)
assert min(2*t-t*t/2 for t in [Fraction(12,25),Fraction(9,4)]) == Fraction(528,625)
print('PASS. The nine-bond symbol, radius inequality, and gap constant agree exactly.')
