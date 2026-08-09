#!/usr/bin/env python3
"""Independent symbolic checks for the algebra in part2.txt.

Requires SymPy.  The script verifies:
  * the seven-wedge expansion (54), modulo the stated exact term;
  * the associativity equations (E1)--(E8) from the multiplication table (34);
  * the formula (56)--(57) for Sigma and B_raw;
  * the resultant factorizations (45)--(46);
  * the simplification (59) and the three nonzero case values (60)--(62).
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# 1. Verify the edge expansion (54).
# ---------------------------------------------------------------------------
ai, aj, Ai, Aj, Bi, Bj = sp.symbols("ai aj Ai Aj Bi Bj")


def wedge(ui, uj, vi, vj):
    return sp.expand(ui * vj - vi * uj)


edge = sp.expand((ai - aj) ** 3 * (Ai - Aj) * (Bi - Bj))
formula_54 = (
    wedge(Bi, Bj, ai**3 * Ai, aj**3 * Aj)
    + wedge(Ai, Aj, ai**3 * Bi, aj**3 * Bj)
    - wedge(Ai * Bi, Aj * Bj, ai**3, aj**3)
    + 3 * wedge(ai, aj, ai**2 * Ai * Bi, aj**2 * Aj * Bj)
    - 3 * wedge(ai * Bi, aj * Bj, ai**2 * Ai, aj**2 * Aj)
    - 3 * wedge(ai * Ai, aj * Aj, ai**2 * Bi, aj**2 * Bj)
    + 3 * wedge(ai * Ai * Bi, aj * Aj * Bj, ai**2, aj**2)
)
exact_term = ai**3 * Ai * Bi - aj**3 * Aj * Bj
assert sp.expand(edge - formula_54 - exact_term) == 0


# ---------------------------------------------------------------------------
# 2. Build the five-dimensional coordinatewise algebra from table (34).
# ---------------------------------------------------------------------------
ONE, a, b, A, B = range(5)
t, m, n, N, ell, Ell, v, V, k, K = sp.symbols(
    "t m n N ell Ell v V k K", commutative=True
)

Vector = List[sp.Expr]


def zero() -> Vector:
    return [sp.Integer(0)] * 5


def basis(i: int) -> Vector:
    out = zero()
    out[i] = sp.Integer(1)
    return out


def vector(**kwargs) -> Vector:
    names = {"one": ONE, "a": a, "b": b, "A": A, "B": B}
    out = zero()
    for name, value in kwargs.items():
        out[names[name]] = sp.sympify(value)
    return out


products: Dict[Tuple[int, int], Vector] = {}


def set_product(i: int, j: int, value: Vector) -> None:
    products[tuple(sorted((i, j)))] = value


for i, name in enumerate(("one", "a", "b", "A", "B")):
    set_product(ONE, i, vector(**{name: 1}))

set_product(a, A, vector(one=sp.Rational(2, 5), b=t, B=t))
set_product(b, B, vector(one=sp.Rational(2, 5), b=-t, B=-t))
set_product(a, a, vector(b=m / 2, A=n / 2, B=ell))
set_product(A, A, vector(B=m / 2, a=N / 2, b=Ell))
set_product(a, B, vector(a=t, b=N / 4, A=m / 2))
set_product(A, b, vector(A=t, B=n / 4, a=m / 2))
set_product(a, b, vector(a=t, A=ell, B=v / 2))
set_product(A, B, vector(A=t, a=Ell, b=V / 2))
set_product(b, b, vector(a=n / 4, b=-t, A=v / 2, B=k / 2))
set_product(B, B, vector(A=N / 4, B=-t, a=V / 2, b=K / 2))


def add(x: Vector, y: Vector) -> Vector:
    return [sp.expand(xi + yi) for xi, yi in zip(x, y)]


def scale(c, x: Vector) -> Vector:
    return [sp.expand(c * xi) for xi in x]


def multiply(x: Vector, y: Vector) -> Vector:
    out = zero()
    for i, xi in enumerate(x):
        if xi == 0:
            continue
        for j, yj in enumerate(y):
            if yj == 0:
                continue
            out = add(out, scale(xi * yj, products[tuple(sorted((i, j)))]))
    return [sp.expand(z) for z in out]


E = [basis(i) for i in range(5)]


def difference(x: Vector, y: Vector) -> Vector:
    return [sp.factor(z) for z in add(x, scale(-1, y))]


# The exact coefficient identities stated as E1--E8.
E1 = 3 * m * n - 4 * t * v
E2 = N * v + 16 * ell * t + 8 * m * t
E3 = -2 * k * m + 24 * ell * t - n**2
E4 = -k * t + ell * m + 4 * t**2
E5 = N * k + m * n + 8 * t * v
E6 = 20 * Ell * ell + 5 * m**2 + 5 * N * n - 40 * t**2 - 8
E7 = 16 * Ell * ell - 4 * m**2 - N * n + 4 * V * v
E8 = 20 * Ell * ell + 60 * t**2 + 5 * V * v - 8

assoc_1 = difference(multiply(multiply(E[a], E[a]), E[b]), multiply(E[a], multiply(E[a], E[b])))
assert sp.factor(assoc_1[a] - E1 / 8) == 0
assert sp.factor(assoc_1[b] + E2 / 8) == 0
assert sp.factor(assoc_1[B] + E3 / 8) == 0

assoc_2 = difference(multiply(multiply(E[a], E[b]), E[b]), multiply(E[a], multiply(E[b], E[b])))
assert sp.factor(assoc_2[a] - E4 / 2) == 0
assert sp.factor(assoc_2[b] + E5 / 8) == 0

assoc_3 = difference(multiply(multiply(E[a], E[a]), E[A]), multiply(E[a], multiply(E[a], E[A])))
assert sp.factor(assoc_3[a] - E6 / 20) == 0

assoc_4 = difference(multiply(multiply(E[a], E[b]), E[A]), multiply(E[a], multiply(E[b], E[A])))
assert sp.factor(assoc_4[b] - E7 / 16) == 0

assoc_5 = difference(multiply(multiply(E[a], E[b]), E[B]), multiply(E[a], multiply(E[b], E[B])))
assert sp.factor(assoc_5[a] - E8 / 20) == 0


# ---------------------------------------------------------------------------
# 3. Verify Sigma = -(i/80) B_raw from (53)--(57).
# ---------------------------------------------------------------------------
def power(x: Vector, exponent: int) -> Vector:
    out = E[ONE]
    for _ in range(exponent):
        out = multiply(out, x)
    return out


def dot(x: Vector, y: Vector):
    # Sum of coordinates equals five times the coefficient of the all-ones vector.
    return sp.expand(5 * multiply(x, y)[ONE])


def wedge_pair(u: Vector, w: Vector, x: Vector, y: Vector):
    return sp.expand(dot(u, x) * dot(w, y) - dot(u, y) * dot(w, x))


va, vb, vA, vB = E[a], E[b], E[A], E[B]
terms = [
    (1, vB, multiply(power(va, 3), vA)),
    (1, vA, multiply(power(va, 3), vB)),
    (-1, multiply(vA, vB), power(va, 3)),
    (3, va, multiply(multiply(power(va, 2), vA), vB)),
    (-3, multiply(va, vB), multiply(power(va, 2), vA)),
    (-3, multiply(va, vA), multiply(power(va, 2), vB)),
    (3, multiply(multiply(va, vA), vB), power(va, 2)),
]

Sigma = sp.Integer(0)
for coefficient, u0, w0 in terms:
    Sigma += sp.I * sp.Rational(1, 2) * coefficient * wedge_pair(vA, va, u0, w0)
    Sigma += sp.I * coefficient * wedge_pair(vB, vb, u0, w0)
Sigma = sp.expand(Sigma)

B_raw = (
    -960 * K * ell * t
    + 160 * Ell * ell * m
    + 320 * ell * t**2
    - 120 * m**3
    - 75 * m * N * n
    - 800 * m * t**2
    + 100 * m * V * v
    + 192 * m
    - 160 * n * t * V
)
assert sp.factor(Sigma + sp.I * B_raw / 80) == 0


# ---------------------------------------------------------------------------
# 4. Verify the elimination in Case III.
# ---------------------------------------------------------------------------
q, h, sigma = sp.symbols("q h sigma", real=True)
f43 = 3 * h * q**2 + 32 * q**2 - 448 * sigma * q - 256
f44 = -9 * h * q - 8 * sigma * h + 112 * sigma * q**2 - 96 * q

resultant_plus = sp.factor(sp.resultant(f43.subs(sigma, 1), f44.subs(sigma, 1), h))
resultant_minus = sp.factor(sp.resultant(f43.subs(sigma, -1), f44.subs(sigma, -1), h))
assert sp.factor(resultant_plus / (16 * (q - 4) * (q + 2) * (21 * q**2 + 42 * q + 16))) == 1
assert sp.factor(resultant_minus / (-16 * (q + 4) * (q - 2) * (21 * q**2 - 42 * q + 16))) == 1


# ---------------------------------------------------------------------------
# 5. Verify the simplification (59) and the three case evaluations.
# ---------------------------------------------------------------------------
abs_n2, abs_v2 = sp.symbols("abs_n2 abs_v2", real=True)
B_after_58 = (
    -960 * ell * (ell * m + 4 * t**2)
    + 160 * ell**2 * m
    + 320 * ell * t**2
    - 120 * m**3
    - 75 * m * abs_n2
    - 800 * m * t**2
    + 100 * m * abs_v2
    + 192 * m
    - 120 * m * abs_n2
)
B_59 = m * (180 * m**2 + 105 * abs_n2 - 4400 * t**2 - 128) - 3520 * ell * t**2

# Eliminate |n|^2 and |v|^2 using E6 and E8, with ell real.
abs_v2_expr = (8 - 20 * ell**2 - 60 * t**2) / 5
abs_n2_expr = (8 - 20 * ell**2 - 5 * m**2 + 40 * t**2) / 5
assert sp.factor(
    (B_after_58 - B_59).subs({abs_v2: abs_v2_expr, abs_n2: abs_n2_expr})
) == 0

# Case I.
assert sp.factor(B_59.subs({t: 0, abs_n2: 0, ell: 0, m**2: sp.Rational(8, 5)}) - 160 * m) == 0

# Case II.
assert sp.factor(
    B_59.subs({t**2: sp.Rational(1, 20), abs_n2: 0, ell: -m / 2, m**2: 1}) + 80 * m
) == 0

# Case III.
r = sp.sqrt(105)
q2 = (1 - r / 21) ** 2
h_value = (-74 + 14 * r) / 3
t2_value = 7 * (5 * r - 33) / 1440
ell_over_m = -(32 + 3 * h_value) / 64
case_iii_coefficient = sp.simplify(
    180 * q2 * t2_value
    + 105 * h_value * t2_value
    - 4400 * t2_value
    - 128
    - 3520 * ell_over_m * t2_value
)
assert sp.simplify(case_iii_coefficient - 60 * (75 - 7 * r)) == 0

print("All symbolic checks passed.")
