# A lower-degree certificate for the dimension-five Makeev candidate

## Status

The multiplication-tensor reduction below is exact.  It replaces the former
25-variable system, whose last equation had degree seven, by either of the
following equivalent computations:

- an exact Magma ideal with 21 variables, 55 quadratic equations, and one
  cubic equation;
- a real interval branch-and-bound computation on the box `[-1,1]^21`.

Neither exhaustive computation was completed in the environment in which
these files were prepared.  Thus the dimension-five theorem is still
conditional.  The conclusive outputs are stated below.

## 1. The multiplication tensor

Let

\[
 e_0=6^{-1/2}{\bf1},e_1=x,e_2=y,e_3=p,e_4=q,e_5=c
\]

be an orthonormal basis of \(\mathbb R^6\), where
\(a=x+iy\), \(b=p+iq\), and the last five vectors lie in
\({\bf1}^{\perp}\).  For \(1\leq i,j,k\leq5\), put

\[
 T_{ijk}=\sum_{\nu=0}^{5}e_i(\nu)e_j(\nu)e_k(\nu).
\]

The tensor is fully symmetric.  Coordinatewise multiplication by \(e_i\)
has, in the displayed basis, the real symmetric matrix

\[
 M_i=
 \begin{pmatrix}
 0&6^{-1/2}\delta_{ij}\\
 6^{-1/2}\delta_{ik}&T_{ijk}
 \end{pmatrix}.
\]

Here the lower-right block is indexed by \(j,k=1,\ldots,5\).  Coordinatewise
multiplication is associative, so

\[
 [M_i,M_j]=0 \qquad(1\leq i<j\leq5).
\]

Also

\[
 \operatorname{tr}M_i=\sum_{j=1}^{5}T_{ijj}=0. \tag{1}
\]

These conditions characterize the simplex frames; they are not merely
necessary.  Indeed, suppose a real symmetric tensor satisfies (1) and the
commutator equations.  The commuting real symmetric matrices \(M_i\) are
simultaneously orthogonally diagonalizable.  Associativity gives

\[
 M_iM_j=\frac{\delta_{ij}}6I+\sum_kT_{ijk}M_k.
\]

Taking traces and using (1) yields

\[
 \operatorname{tr}(M_iM_j)=\delta_{ij}.
\]

Consequently their five eigenvalue vectors are an orthonormal basis of
\({\bf1}^{\perp}\).  Finally,

\[
 T_{ijk}=\operatorname{tr}(M_iM_jM_k)
 =\sum_\nu e_i(\nu)e_j(\nu)e_k(\nu),
\]

which recovers the original frame.  This proves the converse.

Only the lower-right commutator entries need to be imposed.  In zero-based
indices \(i,j,k,l,m\in\{0,\ldots,4\}\), they are

\[
 \sum_m\bigl(T_{ikm}T_{jml}-T_{jkm}T_{iml}\bigr)
 +\frac{\delta_{ki}\delta_{jl}-\delta_{kj}\delta_{il}}6=0. \tag{2}
\]

After duplicates are removed, (2) gives 55 nonzero quadratics.

## 2. Exactness is linear in the tensor

Use zero-based indices \(0,1,2,3,4=x,y,p,q,c\).  The nine real exactness
conditions

\[
 \operatorname{Im}m=0,\qquad4u=\bar n,\qquad r=w=0,\qquad2d=s
\]

become the following linear equations:

\[
\begin{aligned}
 2T_{012}-T_{003}+T_{113}&=0,\\
 2(T_{022}-T_{033}+2T_{123})-T_{000}+3T_{011}&=0,\\
 2(T_{122}-T_{133}-2T_{023})+3T_{001}-T_{111}&=0,\\
 T_{004}-T_{114}=T_{014}&=0,\\
 T_{024}+T_{134}=T_{124}-T_{034}&=0,\\
 2T_{022}+2T_{033}-T_{000}-T_{011}&=0,\\
 2T_{122}+2T_{133}-T_{001}-T_{111}&=0.
\end{aligned} \tag{3}
\]

Equations (1) and (3) are 14 independent linear equations in the 35 entries
of a symmetric five-dimensional cubic tensor.  Exact rational row reduction
leaves 21 variables.  `tensor_polynomials.py` performs that reduction and
generates both certificate systems.

## 3. The obstruction

The moment variables in the earlier reduction are linear functions of
\(T\).  For example,

\[
\begin{aligned}
 s&=\tfrac12(T_{000}+T_{011})
   +\tfrac i2(T_{001}+T_{111}),\\
 m&=T_{002}-T_{112}+2T_{013},\\
 n&=(T_{000}-3T_{011})+i(3T_{001}-T_{111}),\\
 \ell&=\tfrac12(T_{002}-T_{112}-2T_{013})
   +\tfrac i2(T_{003}-T_{113}+2T_{012}).
\end{aligned}
\]

The generator contains the analogous formulas for \(\lambda,v,k,\rho\) and
substitutes them into the previously derived obstruction.  The result is a
cubic polynomial \(E(T)\), and direct expansion gives

\[
 E(T)=8\sqrt2\,\Lambda_g(Y_H(g)) \tag{4}
\]

on the exactness locus.  This scaling was also checked numerically on
independently computed cubic zeros.

Thus the required assertion is precisely that (2) and \(E=0\) have no real
solution after (1) and (3) are eliminated.

## 4. Exact Magma route

Run

```sh
magma makeev_d5_tensor_certificate.m
```

The conclusive output is

```text
Certificate verified: tensor ideal is the unit ideal.
```

This test is stronger than necessary because it excludes complex as well as
real solutions.  It has lower degree and fewer variables than the earlier
certificate.  If Magma reports a proper ideal, that is inconclusive for the
real problem; use the interval route.

## 5. Certified real interval route

Generate the checked polynomial table and compile:

```sh
python3 tensor_polynomials.py --emit-cpp tensor_polynomials.inc
g++ -O3 -DNDEBUG -std=c++20 -frounding-math -ffp-contract=off \
    makeev_d5_interval.cpp -o makeev_d5_interval
```

Do not use `-ffast-math`.

For a short smoke test:

```sh
./makeev_d5_interval --max-boxes 100000 --max-depth 250
```

Exit code 2 and `INCOMPLETE` mean only that the budget was exhausted.  They
are not a mathematical result.  An unlimited single-process run is

```sh
./makeev_d5_interval --max-boxes 0 --max-depth 400
```

The final conclusive line is

```text
Certificate verified: no real tensor zero exists.
```

The search can be divided into independent power-of-two shards.  For example,
the supplied driver runs 256 shards on all available cores:

```sh
python3 run_makeev_interval_shards.py \
    --binary ./makeev_d5_interval --shards 256 --max-boxes 0 --max-depth 400
```

Every shard must exit zero.  Logs are placed in `interval-logs/`.

The verifier uses the following certified exclusions:

1. all 35 tensor entries satisfy \(|T_{ijk}|\leq1\), so every possible
   solution is contained in the initial box;
2. polynomial ranges are evaluated with outward-rounded IEEE-754 intervals;
3. like Taylor monomials are collected before interval evaluation, retaining
   cancellations;
4. a bounded left-null multiplier of the overdetermined Jacobian gives an
   inconsistency test;
5. a rank-20 interval Gauss--Seidel step contracts normal directions while
   retaining the one-dimensional circle orbit.

The code returns success only after its assigned box or shard is empty.

## 6. Reconnaissance performed

- On actual numerical cubic zeros, the nine tensor equations agreed with the
  triangle-cycle equations to numerical precision, and (4) had ratio
  `11.313708498... = 8 sqrt(2)`.
- Thirty numerical least-squares searches for a common tensor zero did not
  find one; the smallest residual norm was approximately `0.08565`.
- A limited outward-rounded run processed 20,000 boxes (9,753 range
  exclusions and 235 overdetermined-consistency exclusions), but an
  exhaustive run was not completed here.

Only one of the two conclusive certificate outputs above upgrades the
dimension-five argument to an unconditional proof.
