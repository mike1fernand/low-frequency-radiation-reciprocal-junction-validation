"""Verify that reference-frequency shells lie inside |a k| < 3/2.

The dimensionless equal-stiffness bulk has three axial and six face-diagonal
bond families. At reference frequency 6/5 its resonant eigenvalue is 36/25.
Cubic symmetry permits the ordered first octant. The independent cell tests
below use integer and rational arithmetic, without a floating-point library.
Run with Python 3. The printed cell counts establish finite covering only when
combined with the Fourier symbol and monotonicity arguments in the review.
"""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial


WEIGHT_SCALE = 10**8
RESONANT_EIGENVALUE = Fraction(36, 25)
PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(22, 7)


def endpoint_weight(phase):
    """Integer lower bound for WEIGHT_SCALE * 2(1 - cos(phase))."""
    # For |phase| <= 44/7 the Taylor remainder after the negative term of
    # degree 20 begins positive and alternates with decreasing magnitudes.
    polynomial = 2 * sum(
        (-1)**(n + 1) * phase**(2*n) / factorial(2*n) for n in range(1, 11)
    )
    return max(0, WEIGHT_SCALE * polynomial.numerator // polynomial.denominator)


def interval_weight(start, stop, step, endpoints):
    if start <= 0 <= stop:
        return 0
    left, right = sorted((abs(start), abs(stop)))
    if right * step >= 2 * PI_LOWER:
        return 0  # conservatively allow the minimum near 2 pi
    # In (0,2 pi), the only interior extremum is the maximum at pi.
    return min(endpoints[left], endpoints[right])


def excludes_shell(cell, interval):
    # The sum of the lower weights times unit bond dyads is a matrix lower
    # bound. Use 2*WEIGHT_SCALE times its difference from (36/25)I, so every
    # entry, including the 1/2 of each face-diagonal dyad, is an integer.
    i, j, k = cell
    scale_threshold = 2 * WEIGHT_SCALE * RESONANT_EIGENVALUE
    if scale_threshold.denominator != 1:
        raise AssertionError('Threshold scaling must be exact.')
    a, b, c = [2 * interval[n, n+1] - int(scale_threshold) for n in cell]

    def diagonal_pair(left, right):
        plus = interval[left+right, left+right+2]
        minus = interval[left-right-1, left-right+1]
        return plus + minus, plus - minus

    xy, e = diagonal_pair(i, j)
    xz, f = diagonal_pair(i, k)
    yz, g = diagonal_pair(j, k)
    a += xy + xz
    b += xy + yz
    c += xz + yz
    determinant = a*b*c + 2*e*f*g - a*g*g - b*f*f - c*e*e
    return a > 0 and a*b - e*e > 0 and determinant > 0


def certify(subdivisions=48, radius=Fraction(3, 2)):
    step = PI_UPPER / subdivisions
    endpoints = [endpoint_weight(n * step) for n in range(2*subdivisions + 1)]
    interval = {
        (start, stop): interval_weight(start, stop, step, endpoints)
        for stop in range(-subdivisions, 2*subdivisions + 1)
        for start in range(max(-subdivisions, stop-2), stop+1)
    }
    interior = 0
    excluded = 0
    for cell in combinations_with_replacement(range(subdivisions), 3):
        max_norm_squared = step**2 * sum((n+1)**2 for n in cell)
        if max_norm_squared < radius**2:
            interior += 1
        elif excludes_shell(cell, interval):
            excluded += 1
        else:
            raise ArithmeticError(f'The lower bound did not exclude the shell in cell {cell}.')
    return interior, excluded


if __name__ == '__main__':
    interior, excluded = certify()
    print(f'PASS. {interior + excluded} cells covered. {interior} lie strictly inside')
    print(f'radius 3/2. In all {excluded} others, D - (36/25)I is positive definite.')
