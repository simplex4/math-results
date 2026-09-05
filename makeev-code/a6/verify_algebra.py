#!/usr/bin/env python3
"""Exact algebra regression checks for the accompanying verification draft.

Requires Python 3.10+ and SymPy. Run: python verify_algebra.py

These checks do not verify the normal-form classification, the spectral
fiber lemma, semialgebraic dimension arguments, or the complete proof.
No numerical tolerances, random trials, or floating-point arithmetic are used.
"""
from __future__ import annotations

import sys
from collections.abc import Sequence
from typing import Any

try:
    import sympy as sp
except ImportError as exc:
    raise SystemExit("SymPy is required. Install it with: python -m pip install sympy") from exc

SQRT2 = sp.sqrt(2)
CHECKS = 0


def zero(expression: Any, name: str) -> None:
    """Require an exact zero; raise an informative exception on failure."""
    global CHECKS
    result = sp.simplify(sp.expand(expression))
    if result != 0:
        raise AssertionError(f"{name}: nonzero remainder {result}")
    CHECKS += 1


def zero_matrix(matrix: sp.MatrixBase, name: str) -> None:
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            zero(matrix[i, j], f"{name}, entry ({i},{j})")


def Q(point: Sequence[Any]) -> Any:
    x, y, u, v = point
    return (x*x - y*y)*u + 2*x*y*v


def H(point: Sequence[Any]) -> Any:
    x, y, u, v = point
    return (x*x + y*y)*(2*x*y*u - (x*x - y*y)*v)


def gradient_Q(point: Sequence[Any]) -> list[Any]:
    x, y, u, v = point
    return [2*x*u + 2*y*v, -2*y*u + 2*x*v, x*x - y*y, 2*x*y]


def edge_X(a: Sequence[Any], b: Sequence[Any]) -> Any:
    x, y, u, v = a
    xx, yy, uu, vv = b
    return -x*yy + y*xx - 2*u*vv + 2*v*uu


def normalized_difference(a: Sequence[Any], b: Sequence[Any]) -> list[Any]:
    return [(aa - bb)/SQRT2 for aa, bb in zip(a, b)]


def conjugate(point: Sequence[Any]) -> list[Any]:
    x, y, u, v = point
    return [x, -y, u, -v]


def check_polynomial_identities() -> None:
    a = sp.symbols("x y u v", real=True)
    b = sp.symbols("X Y U V", real=True)
    x, y, u, v = a
    generator = (-y, x, -2*v, 2*u)
    for polynomial, name in ((Q, "Q"), (H, "H")):
        zero(polynomial([-z for z in a]) + polynomial(a), f"oddness of {name}")
        zero(sum(sp.diff(polynomial(a), z)*dz for z, dz in zip(a, generator)),
             f"infinitesimal invariance of {name}")
    zero(Q(conjugate(a)) - Q(a), "Q under conjugation")
    zero(H(conjugate(a)) + H(a), "H under conjugation")
    zero(edge_X(conjugate(a), conjugate(b)) + edge_X(a, b),
         "annihilator under conjugation")
    zero(Q([aa - bb for aa, bb in zip(a, b)]) - Q(a) + Q(b)
         + sum(ga*bb for ga, bb in zip(gradient_Q(a), b))
         - sum(aa*gb for aa, gb in zip(a, gradient_Q(b))),
         "cubic expansion")
    dx, dy, du, dv = [aa - bb for aa, bb in zip(a, b)]
    displayed = (dx*dx + dy*dy)*(2*dx*dy*du - (dx*dx - dy*dy)*dv)/(4*SQRT2)
    zero(H(normalized_difference(a, b)) - displayed, "normalized quintic edge formula")
    print("PASS: oddness, infinitesimal invariance, conjugation, cubic and edge identities")


def check_edge_group_sums() -> None:
    k, ell = sp.symbols("k ell", positive=True)
    gamma = sp.symbols("gamma", real=True, nonzero=True)
    omega, w = sp.symbols("omega w", real=True)
    inv_sum = 1/k + 1/ell
    u1, u2 = 1/(k*gamma), -1/(ell*gamma)
    cross = 0
    first_outside = 0
    second_outside = 0
    for sign_a in (-1, 1):
        first = [sign_a/sp.sqrt(k), 0, u1, omega]
        first_outside += k/sp.Integer(2)*edge_X(first, [0, 0, 0, w])*H(
            normalized_difference(first, [0, 0, 0, w]))
        for sign_b in (-1, 1):
            second = [0, sign_b/sp.sqrt(ell), u2, omega]
            cross += k*ell/sp.Integer(4)*edge_X(first, second)*H(
                normalized_difference(first, second))
    for sign_b in (-1, 1):
        second = [0, sign_b/sp.sqrt(ell), u2, omega]
        second_outside += ell/sp.Integer(2)*edge_X(second, [0, 0, 0, w])*H(
            normalized_difference(second, [0, 0, 0, w]))
    cross = sp.simplify(cross)
    zero(cross - inv_sum**2/(2*SQRT2*gamma), "cross-edge sum before gamma relation")
    zero(first_outside - w*(omega - w)/(2*SQRT2*k**2*gamma), "K1-to-outside sum")
    zero(second_outside - w*(omega - w)/(2*SQRT2*ell**2*gamma), "K2-to-outside sum")
    # The common moment over outside indices is -1, by zero mean and unit norm.
    total = cross - (1/k**2 + 1/ell**2)/(2*SQRT2*gamma)
    zero(total - 1/(SQRT2*k*ell*gamma), "final obstruction")
    # The displayed cross sum gamma^3/(2 sqrt(2)) uses gamma^2 = 1/k + 1/ell.
    zero((cross - gamma**3/(2*SQRT2))*2*SQRT2*gamma
         - (inv_sum - gamma**2)*(inv_sum + gamma**2),
         "reduction using the gamma-squared relation")
    print("PASS: all three nonzero edge-group sums and the obstruction simplification")


def make_frame(N: int, k: int, ell: int, gamma_sign: int) -> tuple[sp.Matrix, Any]:
    """Construct four exact centered orthonormal columns in normal form.

    The outside entries of v are constant in these test examples.
    Any four centered orthonormal columns extend to a full simplex frame.
    """
    outside = N - k - ell
    if min(k, ell) < 2 or k % 2 or ell % 2 or outside < 1:
        raise ValueError("Even support sizes at least two and an outside index are required.")
    if gamma_sign not in (-1, 1):
        raise ValueError("gamma_sign must be -1 or 1.")
    gamma = gamma_sign*sp.sqrt(sp.Rational(1, k) + sp.Rational(1, ell))
    omega = sp.sqrt(sp.Rational(outside, (k + ell)*N))
    x = [1/sp.sqrt(k)]*(k//2) + [-1/sp.sqrt(k)]*(k//2) + [sp.Integer(0)]*(N - k)
    y = ([sp.Integer(0)]*k + [1/sp.sqrt(ell)]*(ell//2)
         + [-1/sp.sqrt(ell)]*(ell//2) + [sp.Integer(0)]*outside)
    u = [(xx*xx - yy*yy)/gamma for xx, yy in zip(x, y)]
    v = [omega]*(k + ell) + [-(k + ell)*omega/outside]*outside
    return sp.Matrix.hstack(*[sp.Matrix(column) for column in (x, y, u, v)]), gamma


def check_frame(N: int, k: int, ell: int, gamma_sign: int) -> None:
    P, gamma = make_frame(N, k, ell, gamma_sign)
    name = f"N={N}, k1={k}, k2={ell}, sign(gamma)={gamma_sign:+d}"
    zero_matrix(P.T*P - sp.eye(4), f"orthonormal columns: {name}")
    zero_matrix(sp.ones(1, N)*P, f"zero mean: {name}")
    rows = [list(P.row(i)) for i in range(N)]
    G = sp.Matrix([gradient_Q(row) for row in rows])
    S = P.T*G
    zero_matrix(S - S.T, f"symmetric gradient coefficients: {name}")
    zero_matrix(G - P*S, f"gradient matrix equation: {name}")
    Y = sp.zeros(N)
    lambda_H = sp.Integer(0)
    lambda_Q = sp.Integer(0)
    for i in range(N):
        for j in range(i + 1, N):
            difference = normalized_difference(rows[i], rows[j])
            Y[i, j] = Q(difference)
            Y[j, i] = -Y[i, j]
            xx = edge_X(rows[i], rows[j])
            lambda_H += xx*H(difference)
            lambda_Q += xx*Q(difference)
    Pi = sp.eye(N) - sp.ones(N)/N
    zero_matrix(Pi*Y*Pi, f"exact cubic labels: {name}")
    zero(lambda_Q, f"cubic annihilation: {name}")
    zero(lambda_H - 1/(SQRT2*k*ell*gamma), f"quintic obstruction: {name}")
    print(f"PASS: exact normal-form frame ({name})")


def main() -> int:
    print(f"Python {sys.version.split()[0]}; SymPy {sp.__version__}")
    print("Exact algebra regression checks; not a verification of the full proof.")
    check_polynomial_identities()
    check_edge_group_sums()
    for N, k, ell in ((5, 2, 2), (6, 2, 2), (7, 2, 4), (9, 4, 4), (11, 4, 6)):
        for sign in (-1, 1):
            check_frame(N, k, ell, sign)
    print(f"SUCCESS: {CHECKS} exact zero assertions passed.")
    print("The classification and dimension arguments still require independent human verification.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
