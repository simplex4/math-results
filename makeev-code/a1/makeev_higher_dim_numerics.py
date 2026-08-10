#!/usr/bin/env python3
"""Numerical reconnaissance for the same Makeev candidate in d >= 4.

The script searches the Stiefel manifold of complex pairs (a,b) representing
the first four coordinates of a rotated regular simplex.  It is evidence, not
a proof: only the accompanying Magma unit-ideal calculation can certify d=5.

Examples:
    python makeev_higher_dim_numerics.py --d 5 --starts 100
    python makeev_higher_dim_numerics.py --d 8 --starts 50 --seed 8
"""

from __future__ import annotations

import argparse
import numpy as np
from scipy.optimize import least_squares


def unpack(v: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    rows = v.reshape(4, n)
    return rows[0] + 1j * rows[1], rows[2] + 1j * rows[3]


def frame_residual(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Equations saying Re/Im(a,b) are orthonormal in 1^perp."""
    complex_eqs = (
        np.sum(a), np.sum(b), np.sum(a * a), np.sum(b * b),
        np.sum(a * b), np.sum(a * np.conj(b)),
    )
    out: list[float] = []
    for value in complex_eqs:
        out.extend((float(value.real), float(value.imag)))
    out.extend((float(np.vdot(a, a).real - 2),
                float(np.vdot(b, b).real - 2)))
    return np.asarray(out)


def random_frame(rng: np.random.Generator, n: int) -> np.ndarray:
    one = np.ones((n, 1)) / np.sqrt(n)
    raw = rng.normal(size=(n, 4))
    raw -= one @ (one.T @ raw)
    q, _ = np.linalg.qr(raw)
    return np.vstack((q[:, 0], q[:, 1], q[:, 2], q[:, 3])).ravel()


def q_value(da: complex, db: complex) -> float:
    return float(np.real(da * da * np.conj(db)))


def r_value(da: complex, db: complex) -> float:
    return float(np.imag(da * da * np.conj(db)))


def q_cycle_residual(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Independent triangle sums for the Q-label matrix."""
    n = len(a)
    labels = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            labels[i, j] = q_value((a[i] - a[j]) / np.sqrt(2),
                                   (b[i] - b[j]) / np.sqrt(2))
            labels[j, i] = -labels[i, j]
    return np.asarray([
        labels[0, i] + labels[i, j] + labels[j, 0]
        for i in range(1, n) for j in range(i + 1, n)
    ])


def lambda_h(a: np.ndarray, b: np.ndarray) -> float:
    """The exact normalized value Lambda_g(Y_H(g))."""
    total = 0.0
    n = len(a)
    for i in range(n):
        for j in range(i + 1, n):
            xij = (np.imag(a[i] * np.conj(a[j]))
                    + 2 * np.imag(b[i] * np.conj(b[j])))
            da = (a[i] - a[j]) / np.sqrt(2)
            db = (b[i] - b[j]) / np.sqrt(2)
            h = abs(da) ** 2 * r_value(da, db)
            total += xij * h
    return float(total)


def cyclic_root(n: int) -> tuple[np.ndarray, np.ndarray]:
    zeta = np.exp(2j * np.pi / n)
    powers = zeta ** np.arange(n)
    scale = np.sqrt(2 / n)
    return scale * powers, scale * powers**2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d", type=int, required=True)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--max-nfev", type=int, default=10000)
    args = parser.parse_args()
    if args.d < 4:
        parser.error("this construction needs d >= 4")

    n = args.d + 1
    rng = np.random.default_rng(args.seed)

    ac, bc = cyclic_root(n)
    print("cyclic root:",
          f"cycle_norm={np.linalg.norm(q_cycle_residual(ac, bc)):.3e}",
          f"LambdaH={lambda_h(ac, bc):+.12e}")

    def residual(v: np.ndarray) -> np.ndarray:
        a, b = unpack(v, n)
        return np.r_[frame_residual(a, b), q_cycle_residual(a, b)]

    obstruction_values: list[float] = []
    for start in range(args.starts):
        solution = least_squares(
            residual, random_frame(rng, n), method="trf", jac="3-point",
            ftol=1e-13, xtol=1e-13, gtol=1e-13,
            max_nfev=args.max_nfev,
        )
        error = float(np.linalg.norm(solution.fun))
        a, b = unpack(solution.x, n)
        value = lambda_h(a, b)
        if error < 1e-8:
            obstruction_values.append(abs(value))
        rank = np.linalg.matrix_rank(solution.jac, 1e-7)
        print(f"{start:4d} residual={error:.3e}  LambdaH={value:+.12e}  "
              f"rank={rank}  nfev={solution.nfev}", flush=True)

    print(f"converged roots: {len(obstruction_values)}/{args.starts}")
    if obstruction_values:
        print(f"smallest sampled |LambdaH|: {min(obstruction_values):.12e}")


if __name__ == "__main__":
    main()

