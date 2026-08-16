#!/usr/bin/env python3
"""Build exact multiplication-tensor certificates in every dimension d >= 4.

For a fixed dimension, the script eliminates the trace and cubic-exactness
linear equations, forms the quadratic commutator equations for the simplex
multiplication tensor, appends the cubic transversality obstruction, and emits
a self-contained Magma program.  It uses only the Python standard library.

The default is d=5 for backward compatibility.  In that dimension the Magma
output is byte-for-byte identical to the certificate that was verified in
14 seconds.  A polynomial is represented by a sparse dictionary whose
monomial is a sorted tuple of indices of the free tensor variables.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import reduce
from itertools import combinations, combinations_with_replacement
from math import gcd
import argparse


@dataclass(frozen=True)
class Poly:
    terms: dict[tuple[int, ...], Q]

    def __post_init__(self):
        object.__setattr__(self, "terms", {m: Q(c) for m, c in self.terms.items() if c})

    @staticmethod
    def c(c=0):
        c = Q(c)
        return Poly({(): c} if c else {})

    @staticmethod
    def x(i: int):
        return Poly({(i,): Q(1)})

    def __add__(self, other):
        other = as_poly(other); d = self.terms.copy()
        for m, c in other.terms.items():
            d[m] = d.get(m, Q(0)) + c
        return Poly(d)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-as_poly(other))

    def __rsub__(self, other):
        return as_poly(other) - self

    def __mul__(self, other):
        other = as_poly(other); d: dict[tuple[int, ...], Q] = {}
        for m, c in self.terms.items():
            for n, e in other.terms.items():
                k = tuple(sorted(m + n)); d[k] = d.get(k, Q(0)) + c * e
        return Poly(d)

    __rmul__ = __mul__

    def __truediv__(self, c):
        c = Q(c)
        return Poly({m: a/c for m, a in self.terms.items()})


def as_poly(x) -> Poly:
    return x if isinstance(x, Poly) else Poly.c(x)


def tensor_coordinates(dimension: int):
    triples = list(combinations_with_replacement(range(dimension), 3))
    return triples, {triple: i for i, triple in enumerate(triples)}


def linear_rows(dimension: int, tindex: dict[tuple[int, int, int], int]) -> list[list[Q]]:
    """Trace and exactness equations in symmetric-tensor coordinates."""
    rows: list[list[Q]] = []
    n_tensor = (dimension * (dimension + 1) * (dimension + 2)) // 6

    def tidx(i, j, k):
        return tindex[tuple(sorted((i, j, k)))]

    def row(terms):
        r = [Q(0)] * n_tensor
        for c, ijk in terms:
            r[tidx(*ijk)] += Q(c)
        return r

    # tr(M_i)=sum_j T_ijj=0.
    for i in range(dimension):
        rows.append(row([(1, (i, j, j)) for j in range(dimension)]))

    # The five equations common to every d >= 4: Im(m)=0,
    # 4u=conj(n), and 2d_0=s (the symbol d_0 is a moment, not dimension).
    rows.extend([
        row([(2,(0,1,2)),(-1,(0,0,3)),(1,(1,1,3))]),
        row([(2,(0,2,2)),(-2,(0,3,3)),(4,(1,2,3)),(-1,(0,0,0)),(3,(0,1,1))]),
        row([(2,(1,2,2)),(-2,(1,3,3)),(-4,(0,2,3)),(3,(0,0,1)),(-1,(1,1,1))]),
    ])

    # For each extra real coordinate c_r, exactness also says
    # sum c_r a^2=sum c_r a*conj(b)=0 (four real equations).
    for r in range(4, dimension):
        rows.extend([
            row([(1,(0,0,r)),(-1,(1,1,r))]),
            row([(2,(0,1,r))]),
            row([(1,(0,2,r)),(1,(1,3,r))]),
            row([(1,(1,2,r)),(-1,(0,3,r))]),
        ])

    rows.extend([
        row([(2,(0,2,2)),(2,(0,3,3)),(-1,(0,0,0)),(-1,(0,1,1))]),
        row([(2,(1,2,2)),(2,(1,3,3)),(-1,(0,0,1)),(-1,(1,1,1))]),
    ])
    assert len(rows) == 5 * dimension - 11
    return rows


def rref(rows: list[list[Q]]):
    a = [r[:] for r in rows]; pivots = []; r = 0
    for col in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][col]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        d = a[r][col]; a[r] = [x/d for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][col]:
                d = a[i][col]; a[i] = [x-d*y for x, y in zip(a[i], a[r])]
        pivots.append(col); r += 1
        if r == len(a):
            break
    return a, pivots


def tensor_substitution(dimension: int, triples, tindex):
    rr, pivots = rref(linear_rows(dimension, tindex))
    assert len(pivots) == 5 * dimension - 11
    free = [i for i in range(len(triples)) if i not in pivots]
    pos = {old: new for new, old in enumerate(free)}
    t = [None] * len(triples)
    for old in free:
        t[old] = Poly.x(pos[old])
    for row, pivot in zip(rr, pivots):
        t[pivot] = sum(((-row[j]) * t[j] for j in free if row[j]), Poly.c())
    assert all(x is not None for x in t)
    return t, free


def cadd(a, b): return a[0]+b[0], a[1]+b[1]
def cscale(c, a): return c*a[0], c*a[1]
def cmul(a, b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def cconj(a): return a[0], -a[1]
def cprod(*args): return reduce(cmul, args, (Poly.c(1), Poly.c(0)))


def normalize(p: Poly) -> Poly:
    if not p.terms:
        return p
    den_lcm = 1
    for c in p.terms.values():
        den_lcm = den_lcm * c.denominator // gcd(den_lcm, c.denominator)
    ints = [int(c * den_lcm) for c in p.terms.values()]
    g = reduce(gcd, (abs(x) for x in ints if x))
    q = Poly({m: c * Q(den_lcm, g) for m, c in p.terms.items()})
    lead = min(q.terms)
    return -q if q.terms[lead] < 0 else q


def build_system(dimension: int = 5):
    if dimension < 4:
        raise ValueError("dimension must be at least 4")
    triples, tindex = tensor_coordinates(dimension)
    flat, free = tensor_substitution(dimension, triples, tindex)

    def T(i, j, k):
        return flat[tindex[tuple(sorted((i, j, k)))]]

    equations = []
    for i, j in combinations(range(dimension), 2):
        for k, ell in combinations(range(dimension), 2):
            p = sum(
                T(i,k,m)*T(j,m,ell)-T(j,k,m)*T(i,m,ell)
                for m in range(dimension)
            )
            p += Q(
                (k == i) * (j == ell) - (k == j) * (i == ell),
                dimension + 1,
            )
            equations.append(normalize(p))
    equations = list({tuple(sorted(p.terms.items())): p for p in equations if p.terms}.values())

    s = ((T(0,0,0)+T(0,1,1))/2, (T(0,0,1)+T(1,1,1))/2)
    m = (T(0,0,2)-T(1,1,2)+2*T(0,1,3), Poly.c(0))
    nu = (T(0,0,0)-3*T(0,1,1), 3*T(0,0,1)-T(1,1,1))
    ell = ((T(0,0,2)-T(1,1,2)-2*T(0,1,3))/2,
           (T(0,0,3)-T(1,1,3)+2*T(0,1,2))/2)
    lam = ((T(0,0,2)+T(1,1,2))/2, -(T(0,0,3)+T(1,1,3))/2)
    vv = (T(0,2,2)-T(0,3,3)-2*T(1,2,3),
          T(1,2,2)-T(1,3,3)+2*T(0,2,3))
    kk = (T(2,2,2)-3*T(2,3,3), 3*T(2,2,3)-T(3,3,3))
    rho = (T(2,2,2)+T(2,3,3), T(2,2,3)+T(3,3,3))
    terms = [
        cscale(Q(-3,2),cprod(m,m,m)),
        cscale(Q(3,4),cprod(cconj(nu),m,nu)),
        cscale(-12,cprod(cconj(rho),ell,lam)),
        cscale(-16,cprod(cconj(lam),lam,m)),
        cscale(-2,cprod(cconj(nu),cconj(s),ell)),
        cscale(-6,cprod(rho,lam,m)),
        cscale(-8,cprod(cconj(s),lam,nu)),
        cscale(-8,cprod(lam,s,s)),
        cscale(-1,cprod(cconj(vv),m,vv)),
        cscale(12,cprod(cconj(kk),cconj(lam),ell)),
        cscale(14,cprod(cconj(ell),ell,m)),
        cscale(16,cprod(ell,lam,lam)),
        cscale(4,cprod(cconj(lam),cconj(vv),nu)),
        cscale(5,cprod(cconj(ell),nu,s)),
        cscale(6,cprod(cconj(lam),cconj(rho),m)),
        cscale(6,cprod(cconj(vv),ell,s)),
    ]
    obstruction = normalize(sum((z[0] for z in terms), Poly.c(0)))
    return equations, obstruction, flat, free, triples


def emit_cpp(path: str, dimension: int = 5) -> None:
    eqs, obs, flat, free, _ = build_system(dimension); polys = eqs + [obs]
    offsets = [0]; terms = []
    for p in polys:
        for monomial, coefficient in sorted(p.terms.items()):
            assert coefficient.denominator == 1 and -32768 <= coefficient.numerator <= 32767
            v = list(monomial) + [255] * (3-len(monomial))
            terms.append((coefficient.numerator, len(monomial), *v))
        offsets.append(len(terms))
    with open(path, "w", encoding="utf-8") as out:
        out.write("// Generated exactly by tensor_polynomials.py.\n")
        out.write("struct Term { short c; unsigned char degree,v0,v1,v2; };\n")
        out.write(f"constexpr int TP_DIM={dimension}, TP_NV={len(free)}, "
                  f"TP_NE={len(polys)}, TP_NT={len(terms)};\n")
        out.write("constexpr int TP_OFF[TP_NE+1]={" + ",".join(map(str, offsets)) + "};\n")
        out.write("constexpr Term TP_TERMS[TP_NT]={\n")
        for c, d, a, b, e in terms:
            out.write(f"  {{{c},{d},{a},{b},{e}}},\n")
        out.write("};\n")
        lin_offsets = [0]; lin_terms = []
        for p in flat:
            for monomial, coefficient in sorted(p.terms.items()):
                assert len(monomial) == 1
                lin_terms.append((coefficient.numerator, coefficient.denominator, monomial[0]))
            lin_offsets.append(len(lin_terms))
        out.write("struct LinTerm { short num,den; unsigned char var; };\n")
        out.write(f"constexpr int TP_NL={len(lin_terms)};\n")
        out.write(f"constexpr int TP_LOFF[{len(flat)+1}]={{" +
                  ",".join(map(str,lin_offsets)) + "};\n")
        out.write("constexpr LinTerm TP_LIN[TP_NL]={\n")
        for n,d,v in lin_terms:
            out.write(f"  {{{n},{d},{v}}},\n")
        out.write("};\n")


def emit_magma(path: str, dimension: int = 5) -> None:
    eqs, obs, _, free, triples = build_system(dimension); polys = eqs + [obs]

    def tensor_name(triple):
        if dimension <= 10:
            return "t" + "".join(map(str, triple))
        return "t_" + "_".join(map(str, triple))

    names = [tensor_name(triples[i]) for i in free]

    def term_text(monomial, coefficient):
        factors = [names[i] for i in monomial]
        if not factors:
            return str(coefficient)
        product = "*".join(factors)
        if coefficient == 1:
            return product
        if coefficient == -1:
            return "-" + product
        return f"({coefficient})*{product}"

    def poly_text(p):
        pieces = [term_text(m, c) for m, c in sorted(p.terms.items())]
        return " + ".join(pieces).replace("+ -", "- ") if pieces else "0"

    with open(path, "w", encoding="utf-8") as out:
        out.write("/*\n  Generated exactly by tensor_polynomials.py.\n"
                  f"  This is the {len(free)}-variable multiplication-tensor test: {len(eqs)}\n"
                  "  quadratic commutator equations and one cubic obstruction.\n"
                  "  The success line proves nonvanishing over C, hence over R.\n*/\n")
        out.write(f"Q := Rationals();\nP := PolynomialRing(Q, {len(free)}, \"grevlex\");\n")
        out.write("AssignNames(~P, [" + ",".join(f'\"{n}\"' for n in names) + "]);\n")
        out.write(f"z := [P.i : i in [1..{len(free)}]];\n")
        for i, name in enumerate(names, 1):
            out.write(f"{name} := z[{i}];\n")
        out.write("eqs := [\n")
        for i, p in enumerate(polys):
            comma = "," if i + 1 < len(polys) else ""
            out.write("  " + poly_text(p) + comma + "\n")
        out.write("];\n")
        out.write(f"assert #eqs eq {len(polys)};\nI := ideal<P | eqs>;\n")
        out.write("SetGBGlobalModular(true);\ntime proper := IsProper(I);\n")
        out.write("assert not proper;\n")
        out.write('print "Certificate verified: tensor ideal is the unit ideal.";\n')


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dimension", type=int, default=5)
    ap.add_argument("--emit-cpp")
    ap.add_argument("--emit-magma")
    args = ap.parse_args()
    if args.dimension < 4:
        ap.error("--dimension must be at least 4")
    if args.emit_cpp:
        emit_cpp(args.emit_cpp, args.dimension)
        raise SystemExit
    if args.emit_magma:
        emit_magma(args.emit_magma, args.dimension)
        raise SystemExit
    eqs, obs, flat, free, triples = build_system(args.dimension)
    print("dimension:", args.dimension)
    print("linear equations:", 5 * args.dimension - 11)
    print("free variables:", len(free), [triples[i] for i in free])
    print("commutator equations:", len(eqs), "terms:", sum(len(p.terms) for p in eqs))
    print("obstruction terms:", len(obs.terms), "degree:", max(map(len, obs.terms)))
    print("largest commutator:", max(map(lambda p: len(p.terms), eqs)))
