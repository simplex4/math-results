#!/usr/bin/env python3
"""Exact algebra regression checks for the accompanying verification draft.

Requires Python 3.10+ and SymPy. Run: python verify_algebra_revised.py

Expanded during the September 4, 2026 review. The original checks are
retained, with additional exact test families. This is a regression suite,
not a formal proof or an exhaustive search over placements.

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
    expression = sp.sympify(expression)
    if expression.has(sp.Float):
        raise AssertionError(f"{name}: inexact floating-point input is not allowed")
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


def make_frame(
    N: int, k: int, ell: int, gamma_sign: int, v_mode: str = "constant"
) -> tuple[sp.Matrix, Any]:
    """Construct four exact centered orthonormal columns in normal form.

    v_mode selects constant, negative_constant, zero, mixed, or negative_mixed.
    The last three modes require at least two outside indices. The zero mode
    tests the omega=0 branch; mixed modes have nonconstant outside entries.
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
    v = sp.Matrix(v)
    if v_mode == "negative_constant":
        v = -v
    elif v_mode in ("zero", "mixed", "negative_mixed"):
        if outside < 2:
            raise ValueError(f"{v_mode} requires at least two outside indices")
        direction = sp.zeros(N, 1)
        direction[k + ell] = 1/SQRT2
        direction[k + ell + 1] = -1/SQRT2
        if v_mode == "zero":
            v = direction
        else:
            cosine = sp.Rational(3 if v_mode == "mixed" else -3, 5)
            v = cosine*v + sp.Rational(4, 5)*direction
    elif v_mode != "constant":
        raise ValueError(f"Unknown v_mode: {v_mode}")
    return sp.Matrix.hstack(*[sp.Matrix(column) for column in (x, y, u, v)]), gamma


def check_frame(N: int, k: int, ell: int, gamma_sign: int, v_mode: str = "constant") -> None:
    P, gamma = make_frame(N, k, ell, gamma_sign, v_mode)
    name = f"N={N}, k1={k}, k2={ell}, sign(gamma)={gamma_sign:+d}, v={v_mode}"
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
    # If k=2, the two first rows differ only in x by sqrt(2), so this
    # really is an exceptional E1 placement in any full-frame extension.
    if k == 2:
        zero_matrix((P.row(0) - P.row(1))/SQRT2 - sp.Matrix([[1, 0, 0, 0]]),
                    f"E1 root: {name}")
    # In this branch the two outside rows differ only in v by sqrt(2).
    if v_mode == "zero":
        zero_matrix((P.row(k + ell) - P.row(k + ell + 1))/SQRT2
                    - sp.Matrix([[0, 0, 0, 1]]), f"E2 root: {name}")
    print(f"PASS: exact normal-form frame ({name})")



def check_quadratic_inverse() -> None:
    """Verify the evaluation inverse on arbitrary symbolic edge labels."""
    for N in (4, 6):
        W = sp.zeros(N)
        labels = iter(sp.symbols(f"a0:{N*(N-1)//2}"))
        for i in range(N):
            for j in range(i + 1, N):
                W[i, j] = W[j, i] = next(labels)
        Pi = sp.eye(N) - sp.ones(N)/N
        B = -Pi*W*Pi
        zero_matrix(B*sp.ones(N, 1), f"quadratic inverse is centered: N={N}")
        for i in range(N):
            for j in range(i + 1, N):
                zero((B[i, i] + B[j, j])/2 - B[i, j] - W[i, j],
                     f"quadratic inverse edge {i},{j}: N={N}")
    print("PASS: quadratic-evaluation inverse on symbolic labels")


def helmert_frame(N: int) -> sp.Matrix:
    """Return a centered orthonormal basis of 1-perp, as N by (N-1) columns."""
    columns = []
    for j in range(1, N):
        scale = sp.sqrt(j*(j + 1))
        columns.append(sp.Matrix([1/scale]*j + [-j/scale] + [0]*(N-j-1)))
    return sp.Matrix.hstack(*columns)


def check_general_frames() -> None:
    """Check identities away from normal-form frames, not just at cubic zeros."""
    for d in (4, 5, 6):
        base = helmert_frame(d + 1)
        rotation = sp.eye(d)
        for a, b in ((0, 1), (0, 2), (1, 3), (0, d-1)):
            step = sp.eye(d)
            step[a, a] = step[b, b] = sp.Rational(3, 5)
            step[a, b] = -sp.Rational(4, 5)
            step[b, a] = sp.Rational(4, 5)
            rotation = rotation*step
        for mode, full in (("base", base), ("rotated", base*rotation)):
            N = full.rows
            name = f"d={d}, {mode}"
            Pi = sp.eye(N) - sp.ones(N)/N
            zero_matrix(full.T*full - sp.eye(d), f"full-frame columns: {name}")
            zero_matrix(full*full.T - Pi, f"full simplex Gram matrix: {name}")
            P = full[:, :4]
            rows = [list(P.row(i)) for i in range(N)]
            G = sp.Matrix([gradient_Q(row) for row in rows])
            Y = sp.zeros(N)
            X = sp.zeros(N)
            for i in range(N):
                for j in range(i + 1, N):
                    Y[i, j] = Q(normalized_difference(rows[i], rows[j]))
                    Y[j, i] = -Y[i, j]
                    X[i, j] = edge_X(rows[i], rows[j])
                    X[j, i] = -X[i, j]
            zero_matrix(sp.ones(1, N)*G, f"zero mean of gradient: {name}")
            zero_matrix(X*sp.ones(N, 1), f"cycle annihilator: {name}")
            zero_matrix(Pi*Y*Pi - (P*G.T - G*P.T)/(2*SQRT2),
                        f"general cubic residual identity: {name}")
            zero(sum(X[i, j]*Y[i, j] for i in range(N) for j in range(i+1, N)),
                 f"cubic annihilation away from assumed normal forms: {name}")
            print(f"PASS: general simplex-frame identities ({name})")


def check_coordinate_change() -> None:
    """Check that the r,s calculation transforms the polynomial and its gradient."""
    r, s0, u, v = sp.symbols("r s u v", real=True)
    point = [(r+s0)/SQRT2, (s0-r)/SQRT2, u, v]
    q_new = 2*r*s0*u + (s0*s0-r*r)*v
    zero(Q(point) - q_new, "passive coordinate expression for Q")
    # old coordinate column = T * new coordinate column
    T = sp.Matrix([[1/SQRT2, 1/SQRT2, 0, 0],
                   [-1/SQRT2, 1/SQRT2, 0, 0],
                   [0, 0, 1, 0], [0, 0, 0, 1]])
    new_gradient = sp.Matrix([sp.diff(q_new, z) for z in (r, s0, u, v)])
    zero_matrix(T.T*sp.Matrix(gradient_Q(point)) - new_gradient,
                "chain rule for transformed gradient")
    print("PASS: passive coordinate change and transformed gradient")


def check_spectral_examples() -> None:
    """Finite exact regression examples; this does NOT verify the fiber bound."""
    t = sp.symbols("t")
    for d in (4, 6):
        a, b = sp.Rational(1, 3), sp.Rational(9, 2)
        S = sp.diag(a, a, b, b, *[4*j+4 for j in range(1, d-3)])
        samples = {"zero": sp.zeros(d, 1),
                   "one-weight": sp.Matrix([sp.Rational(1, 3)] + [0]*(d-1)),
                   "both-weights": sp.Matrix([sp.Rational(1, 4), 0,
                                                sp.Rational(1, 3), 0] + [0]*(d-4)),
                   "all-clusters": sp.Matrix([sp.Rational(j+1, 10*d) for j in range(d)])}
        chi_S = S.charpoly(t).as_expr()
        for mode, c in samples.items():
            if not (c.dot(c) < 1):
                raise AssertionError("Test vector violates the spectral norm hypothesis")
            B = S + c*c.T
            M = t*sp.eye(d) - B
            chi_B = B.charpoly(t).as_expr()
            # Schur complement: this determinant is chi_B*(1+c^T M^{-1}c).
            bordered = M.row_join(c).col_join((-c.T).row_join(sp.ones(1, 1)))
            zero(bordered.det(method="domain-ge") - chi_S,
                 f"secular determinant sign: d={d}, {mode}")
            gcd = sp.Poly(sp.gcd(chi_B, sp.diff(chi_B, t)), t).monic()
            zero(sp.rem((t-a)*(t-b), gcd.as_expr(), t),
                 f"only unchanged double eigenvalues: d={d}, {mode}")
            zero(sp.gcd(gcd.as_expr(), sp.diff(gcd.as_expr(), t)) - 1,
                 f"no higher spectral multiplicities: d={d}, {mode}")
            for eigenvalue in (a, b):
                kernel = (B-eigenvalue*sp.eye(d)).nullspace()
                if len(kernel) == 2:
                    for basis_vector in kernel:
                        zero(c.dot(basis_vector),
                             f"zero component in repeated eigenspace: d={d}, {mode}")
            print(f"PASS: selected rank-one spectral example (d={d}, {mode})")


def main() -> int:
    print(f"Python {sys.version.split()[0]}; SymPy {sp.__version__}")
    print("Exact algebra regression checks; not a verification of the full proof.")
    check_polynomial_identities()
    check_edge_group_sums()
    check_coordinate_change()
    check_quadratic_inverse()
    check_spectral_examples()
    check_general_frames()
    for N, k, ell in ((5, 2, 2), (6, 2, 2), (7, 2, 4), (9, 4, 4), (11, 4, 6)):
        for sign in (-1, 1):
            check_frame(N, k, ell, sign)
    for N, k, ell in ((6, 2, 2), (8, 2, 4), (10, 4, 4)):
        for sign in (-1, 1):
            for mode in ("zero", "mixed", "negative_mixed"):
                check_frame(N, k, ell, sign, mode)
    check_frame(8, 2, 4, 1, "negative_constant")
    print(f"SUCCESS: {CHECKS} exact zero assertions passed.")
    print("The classification and dimension arguments still require independent human verification.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
