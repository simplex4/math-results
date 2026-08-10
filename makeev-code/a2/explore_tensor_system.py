#!/usr/bin/env python3
"""Numerical reconnaissance for the d=5 multiplication-tensor system."""

from __future__ import annotations

import itertools
import numpy as np
from scipy.optimize import least_squares


TRIPLES = list(itertools.combinations_with_replacement(range(5), 3))
INDEX = {t: k for k, t in enumerate(TRIPLES)}


def tidx(i: int, j: int, k: int) -> int:
    return INDEX[tuple(sorted((i, j, k)))]


def tensor(v: np.ndarray) -> np.ndarray:
    out = np.empty((5, 5, 5))
    for i in range(5):
        for j in range(5):
            for k in range(5):
                out[i, j, k] = v[tidx(i, j, k)]
    return out


def complex_mul(a: tuple[float, float], b: tuple[float, float]):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def complex_conj(a: tuple[float, float]):
    return a[0], -a[1]


def obstruction(T: np.ndarray) -> float:
    # Indices x,y,p,q,c are 0,1,2,3,4.
    s = ((T[0, 0, 0] + T[0, 1, 1]) / 2,
         (T[0, 0, 1] + T[1, 1, 1]) / 2)
    m = T[0, 0, 2] - T[1, 1, 2] + 2 * T[0, 1, 3]
    nu = (T[0, 0, 0] - 3 * T[0, 1, 1],
          3 * T[0, 0, 1] - T[1, 1, 1])
    ell = ((T[0, 0, 2] - T[1, 1, 2] - 2 * T[0, 1, 3]) / 2,
           (T[0, 0, 3] - T[1, 1, 3] + 2 * T[0, 1, 2]) / 2)
    lam = ((T[0, 0, 2] + T[1, 1, 2]) / 2,
           -(T[0, 0, 3] + T[1, 1, 3]) / 2)
    vv = (T[0, 2, 2] - T[0, 3, 3] - 2 * T[1, 2, 3],
          T[1, 2, 2] - T[1, 3, 3] + 2 * T[0, 2, 3])
    kk = (T[2, 2, 2] - 3 * T[2, 3, 3],
          3 * T[2, 2, 3] - T[3, 3, 3])
    rho = (T[2, 2, 2] + T[2, 3, 3],
           T[2, 2, 3] + T[3, 3, 3])
    mr = (m, 0.0)

    def prod(*zs):
        ans = (1.0, 0.0)
        for z in zs:
            ans = complex_mul(ans, z)
        return ans

    terms = [
        (-1.5, prod(mr, mr, mr)),
        (0.75, prod(complex_conj(nu), mr, nu)),
        (-12, prod(complex_conj(rho), ell, lam)),
        (-16, prod(complex_conj(lam), lam, mr)),
        (-2, prod(complex_conj(nu), complex_conj(s), ell)),
        (-6, prod(rho, lam, mr)),
        (-8, prod(complex_conj(s), lam, nu)),
        (-8, prod(lam, s, s)),
        (-1, prod(complex_conj(vv), mr, vv)),
        (12, prod(complex_conj(kk), complex_conj(lam), ell)),
        (14, prod(complex_conj(ell), ell, mr)),
        (16, prod(ell, lam, lam)),
        (4, prod(complex_conj(lam), complex_conj(vv), nu)),
        (5, prod(complex_conj(ell), nu, s)),
        (6, prod(complex_conj(lam), complex_conj(rho), mr)),
        (6, prod(complex_conj(vv), ell, s)),
    ]
    return sum(scale * z[0] for scale, z in terms)


def exactness(T: np.ndarray) -> np.ndarray:
    return np.array([
        2*T[0, 1, 2] - T[0, 0, 3] + T[1, 1, 3],
        2*(T[0, 2, 2] - T[0, 3, 3] + 2*T[1, 2, 3])
        - T[0, 0, 0] + 3*T[0, 1, 1],
        2*(T[1, 2, 2] - T[1, 3, 3] - 2*T[0, 2, 3])
        + 3*T[0, 0, 1] - T[1, 1, 1],
        T[0, 0, 4] - T[1, 1, 4],
        2*T[0, 1, 4],
        T[0, 2, 4] + T[1, 3, 4],
        T[1, 2, 4] - T[0, 3, 4],
        2*T[0, 2, 2] + 2*T[0, 3, 3] - T[0, 0, 0] - T[0, 1, 1],
        2*T[1, 2, 2] + 2*T[1, 3, 3] - T[0, 0, 1] - T[1, 1, 1],
    ])


def matrices(T: np.ndarray) -> list[np.ndarray]:
    out = []
    h = 1 / np.sqrt(6)
    for i in range(5):
        M = np.zeros((6, 6))
        M[0, i+1] = M[i+1, 0] = h
        M[1:, 1:] = T[i]
        out.append(M)
    return out


def residual(v: np.ndarray, include_obstruction: bool = True) -> np.ndarray:
    T = tensor(v)
    Ms = matrices(T)
    trace = np.array([np.trace(M) for M in Ms])
    comm = []
    for i in range(5):
        for j in range(i+1, 5):
            C = Ms[i] @ Ms[j] - Ms[j] @ Ms[i]
            comm.extend(C[np.triu_indices(6, 1)])
    ans = np.r_[trace, exactness(T), np.asarray(comm)]
    if include_obstruction:
        ans = np.r_[ans, obstruction(T)]
    return ans


def frame_tensor(rng: np.random.Generator) -> np.ndarray:
    one = np.ones((6, 1)) / np.sqrt(6)
    raw = rng.normal(size=(6, 5))
    raw -= one @ (one.T @ raw)
    q, _ = np.linalg.qr(raw)
    return np.array([np.sum(q[:, i] * q[:, j] * q[:, k]) for i, j, k in TRIPLES])


def main() -> None:
    rng = np.random.default_rng(20260810)
    for k in range(30):
        x0 = frame_tensor(rng) if k < 10 else rng.uniform(-0.5, 0.5, len(TRIPLES))
        sol = least_squares(residual, x0, jac="3-point", ftol=1e-13,
                            xtol=1e-13, gtol=1e-13, max_nfev=10000)
        print(k, np.linalg.norm(sol.fun), obstruction(tensor(sol.x)), sol.nfev,
              np.linalg.matrix_rank(sol.jac, 1e-8), flush=True)


if __name__ == "__main__":
    main()
