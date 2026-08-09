"""Exact symbolic checks for the d=4 A_4 counterexample proof.

Requires SymPy.  Every computation is over Q with formal conjugate symbols;
the final case computations are over Q(sqrt(105)).

This script checks:
  1. the degree-five wedge expansion (W5);
  2. the associativity equations (E1)--(E8) from the multiplication table;
  3. the formula Sigma = -i B_raw/80 for the transverse functional;
  4. the eliminations and the three nonzero values of B.
"""
from __future__ import annotations

import sympy as sp


def check_wedge_expansion() -> None:
    x, y, A, C, B, D = sp.symbols("x y A C B D")

    def wedge(ui, uj, vi, vj):
        return ui * vj - vi * uj

    lhs = (x - y) ** 3 * (A - C) * (B - D)
    rhs = (
        wedge(B, D, x**3 * A, y**3 * C)
        + wedge(A, C, x**3 * B, y**3 * D)
        - wedge(A * B, C * D, x**3, y**3)
        + 3 * wedge(x, y, x**2 * A * B, y**2 * C * D)
        - 3 * wedge(x * B, y * D, x**2 * A, y**2 * C)
        - 3 * wedge(x * A, y * C, x**2 * B, y**2 * D)
        + 3 * wedge(x * A * B, y * C * D, x**2, y**2)
    )
    # The omitted term is exact: c_i-c_j with c_i=a_i^3 \bar a_i \bar b_i.
    exact = x**3 * A * B - y**3 * C * D
    assert sp.expand(lhs - rhs - exact) == 0


# Basis: 1, a, b, A=bar(a), B=bar(b).
ONE, a, b, A, B = range(5)
t, m = sp.symbols("t m", real=True)
n, N, ell, ELL, v, V, k, K = sp.symbols("n N ell ELL v V k K")


def zvec():
    return [sp.Integer(0)] * 5


def vec(o=0, aa=0, bb=0, AA=0, BB=0):
    return [o, aa, bb, AA, BB]


TABLE = {}


def put(i, j, value):
    TABLE[tuple(sorted((i, j)))] = value


for index in range(5):
    e = zvec()
    e[index] = 1
    put(ONE, index, e)

put(a, A, vec(sp.Rational(2, 5), 0, t, 0, t))
put(b, B, vec(sp.Rational(2, 5), 0, -t, 0, -t))
put(a, a, vec(0, 0, m / 2, n / 2, ell))
put(A, A, vec(0, N / 2, ELL, 0, m / 2))
put(a, B, vec(0, t, N / 4, m / 2, 0))
put(A, b, vec(0, m / 2, 0, t, n / 4))
put(a, b, vec(0, t, 0, ell, v / 2))
put(A, B, vec(0, ELL, V / 2, t, 0))
put(b, b, vec(0, n / 4, -t, v / 2, k / 2))
put(B, B, vec(0, V / 2, K / 2, N / 4, -t))

BASIS = []
for index in range(5):
    e = zvec()
    e[index] = 1
    BASIS.append(e)


def add(x, y):
    return [sp.expand(x[i] + y[i]) for i in range(5)]


def scale(c, x):
    return [sp.expand(c * xi) for xi in x]


def mul(x, y):
    out = zvec()
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            if xi != 0 and yj != 0:
                out = add(out, scale(xi * yj, TABLE[tuple(sorted((i, j)))]))
    return [sp.expand(q) for q in out]


def product(*indices):
    out = BASIS[ONE]
    for index in indices:
        out = mul(out, BASIS[index])
    return out


def associator(x, y, z):
    lhs = mul(mul(BASIS[x], BASIS[y]), BASIS[z])
    rhs = mul(BASIS[x], mul(BASIS[y], BASIS[z]))
    return [sp.factor(lhs[i] - rhs[i]) for i in range(5)]


def check_associativity_equations() -> None:
    aa_b = associator(a, a, b)
    ab_b = associator(a, b, b)
    aa_A = associator(a, a, A)
    ab_A = associator(a, b, A)
    ab_B = associator(a, b, B)

    expected = [
        sp.factor(3 * m * n - 4 * t * v),
        sp.factor(N * v + 16 * ell * t + 8 * m * t),
        sp.factor(-2 * k * m + 24 * ell * t - n**2),
        sp.factor(-k * t + ell * m + 4 * t**2),
        sp.factor(N * k + m * n + 8 * t * v),
        sp.factor(20 * ELL * ell + 5 * N * n + 5 * m**2 - 40 * t**2 - 8),
        sp.factor(16 * ELL * ell - N * n + 4 * V * v - 4 * m**2),
        sp.factor(20 * ELL * ell + 5 * V * v + 60 * t**2 - 8),
    ]
    obtained = [
        sp.factor(8 * aa_b[a]),
        sp.factor(-8 * aa_b[b]),
        sp.factor(8 * aa_b[B]),
        sp.factor(2 * ab_b[a]),
        sp.factor(-8 * ab_b[b]),
        sp.factor(20 * aa_A[a]),
        sp.factor(16 * ab_A[b]),
        sp.factor(20 * ab_B[a]),
    ]
    # Each obtained expression can differ by a harmless nonzero sign.
    for got, want in zip(obtained, expected):
        assert sp.expand(got - want) == 0 or sp.expand(got + want) == 0


def bilinear_dot(x, y):
    # Sum_i x_i y_i = 5 times the coefficient of 1 in the coordinatewise product.
    return sp.expand(5 * mul(x, y)[ONE])


def wedge_pair(u, w, x, y):
    return sp.expand(
        bilinear_dot(u, x) * bilinear_dot(w, y)
        - bilinear_dot(u, y) * bilinear_dot(w, x)
    )


def check_transverse_functional() -> None:
    terms = [
        (1, BASIS[B], product(a, a, a, A)),
        (1, BASIS[A], product(a, a, a, B)),
        (-1, product(A, B), product(a, a, a)),
        (3, BASIS[a], product(a, a, A, B)),
        (-3, product(a, B), product(a, a, A)),
        (-3, product(a, A), product(a, a, B)),
        (3, product(a, A, B), product(a, a)),
    ]

    sigma = 0
    for coefficient, u, w in terms:
        sigma += coefficient * sp.I / 2 * (
            wedge_pair(BASIS[A], BASIS[a], u, w)
            + 2 * wedge_pair(BASIS[B], BASIS[b], u, w)
        )
    sigma = sp.expand(sigma)

    b_raw = (
        -960 * K * ell * t
        + 160 * ELL * ell * m
        + 320 * ell * t**2
        - 120 * m**3
        - 75 * m * N * n
        - 800 * m * t**2
        + 100 * m * V * v
        + 192 * m
        - 160 * n * t * V
    )
    assert sp.expand(sigma + sp.I * b_raw / 80) == 0


def check_case_elimination() -> None:
    q, h = sp.symbols("q h", real=True)
    for sigma in (1, -1):
        f1 = 3 * h * q**2 + 32 * q**2 - 448 * sigma * q - 256
        f2 = -9 * h * q - 8 * sigma * h + 112 * sigma * q**2 - 96 * q
        resultant = sp.factor(sp.resultant(f1, f2, h))
        target = (
            16 * (q - 4) * (q + 2) * (21 * q**2 + 42 * q + 16)
            if sigma == 1
            else -16 * (q - 2) * (q + 4) * (21 * q**2 - 42 * q + 16)
        )
        assert sp.expand(resultant - target) == 0

    root = sp.sqrt(105)
    for sigma in (1, -1):
        q0 = -sigma * (1 - root / 21)
        h0 = (-74 + 14 * root) / 3
        L0 = -q0 * (32 + 3 * h0) / 64  # ell/t
        t2 = sp.simplify(8 / (20 * L0**2 + 5 * q0**2 + 5 * h0 - 40))
        m2 = sp.simplify(q0**2 * t2)
        assert sp.simplify(t2 - 7 * (5 * root - 33) / 1440) == 0
        assert sp.simplify(m2 - (49 * root - 477) / 1080) == 0

        # B/m in Case III.
        coefficient = sp.simplify(
            t2 * (180 * q0**2 + 105 * h0 - 4400)
            - 128
            - 3520 * (L0 / q0) * t2
        )
        assert sp.simplify(coefficient - 60 * (75 - 7 * root)) == 0

    # Cases I and II.
    m_symbol = sp.symbols("m_symbol", nonzero=True, real=True)
    b_case_i = m_symbol * (180 * sp.Rational(8, 5) - 128)
    assert sp.expand(b_case_i - 160 * m_symbol) == 0
    b_case_ii = (
        m_symbol * (180 - 4400 * sp.Rational(1, 20) - 128)
        - 3520 * (-m_symbol / 2) * sp.Rational(1, 20)
    )
    assert sp.expand(b_case_ii + 80 * m_symbol) == 0


if __name__ == "__main__":
    check_wedge_expansion()
    check_associativity_equations()
    check_transverse_functional()
    check_case_elimination()
    print("All exact symbolic checks passed.")
