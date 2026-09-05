# Review of `makeev.tex`

**Review date:** September 4, 2026  
**Materials reviewed:** the complete uploaded `makeev.tex` and `verify_algebra.py`; the relevant passages of the cited Kuperberg paper, Chen preprint, predecessor repository manuscript, and semialgebraic background sources.  
**Deliverables:** `makeev_revised.tex`, its compiled PDF, `verify_algebra_revised.py`, this report, the original and expanded algebra-audit logs, and a unified source diff.

## Assessment

I did not find a fatal mathematical error in the proposed all-dimensional argument. In particular, the spectral-fiber argument and the exceptional normal-form argument withstand the checks described below. I found one literally false normalization sentence, several compressed or potentially misleading explanations, a missing explicit normalization of the polar, and minor editorial issues. The revised source corrects or clarifies these without changing the main theorem, its degree bound, or the existential nature of the perturbation.

This is an AI-assisted mathematical review, not independent human verification, peer review, or proof-assistant certification. The passing symbolic checks do not establish the universal quantifiers in the normal-form proposition or the semialgebraic dimension argument. My assessment of those steps rests on reading and rederiving their arguments, not on extrapolating from the computations.

References to line numbers below refer to the **original uploaded source**, not the revised source. The theorem and equation numbers have been preserved.

## 1. Corrections and clarifications made

### 1.1 The unit-vector sentence is literally false

Original line 115 says that all function evaluations in the paper are at unit vectors. However, the proof evaluates the cubic and its gradient at the simplex vertices

\[
p_i=A(e_i-N^{-1}\mathbf1),\qquad \|p_i\|^2=1-N^{-1}=\frac d{d+1},
\]

which are not unit vectors. The root-label evaluations, by contrast, are at

\[
(p_i-p_j)/\sqrt2=A\alpha_{ij},
\]

which are unit vectors. This is a writing error, not a missing factor in the calculations. I changed the sentence to refer specifically to **root-label evaluations** and explicitly distinguished the auxiliary vertex evaluations. I checked rather than altered the factors \(1/(2\sqrt2)\) in the cubic formulas and \(1/(4\sqrt2)\) in the quintic edge formula.

### 1.2 The conjectured body and its scale are now explicit

The original definition of \(C_d\) is correct, but “appropriately scaled polar” leaves the scale unstated. For the simplex used in the paper, the exact identity is

\[
C_d=\frac1{\sqrt2}(\Delta-\Delta)^\circ,
\qquad \Delta=\operatorname{conv}\{e_i-N^{-1}\mathbf1\}.
\]

The polar is taken in \(V=\mathbf1^\perp\). This follows directly because the defining root inequalities are equivalent to \(|y_i-y_j|\le1/\sqrt2\). The revised text defines the polar, gives this scale, and explicitly distinguishes the polar of the difference body from the difference body itself. It also defines the convention for “universal cover,” allowing translations, rotations, and reflections, consistently with the quantifier over \(O(d)\).

The regular-simplex conjecture in Kuperberg’s Conjecture 1 concerns this dual difference body. Chen’s positive maximal-facet-count theorem concerns other polytopes and a different conjecture of Makeev. I added a sentence to prevent the two conjectures from being conflated.

### 1.3 The spectral multiplicity and residue arguments are expanded

Lemma 4.4 is not invalidated by repeated eigenvalues, but the reason deserves more explicit treatment. The revision explains why every repeated eigenvalue of \(B\) must be an unchanged double eigenvalue of \(S\). It also states precisely that the determinant identity is a rational-function identity, and that its left side is determined by the prescribed spectrum of \(S\), not by knowing the matrix \(S\) itself.

Zero components of \(c\), including \(c=0\), are now explicitly included. For fixed \(B\) and fixed \(\lambda\), the argument gives at most \(2^d\) possible vectors \(c\). This is an expansion of the existing proof, not a new assumption.

### 1.4 “Stratum” does not assert a smooth manifold

The signed orbit-type pieces can be singular or disconnected. The original proof only uses semialgebraic dimension, so this is not a gap. The revised text explicitly says that “stratum” means a piece of a finite semialgebraic partition and that no smoothness, connectedness, or constant individual-root orbit dimension is assumed. Representatives are selected using an explicitly specified lexicographic ordering.

I also made explicit that the exceptional set is semialgebraic and nonempty, and explained why invariant odd functions vanish on its special root directions. The latter is exactly why those directions cannot be handled by arbitrary interpolation.

### 1.5 The normal-form coefficient deductions are displayed

In Proposition 6.1, the symmetric matrices underlying the reduced scalar equations were implicit. I inserted their full forms. In Case 1, the matrix in the ordered coordinates \(x,y,u,v\) is

\[
\begin{pmatrix}
2U&2W&0&2\tau\\
2W&-2a&-2b&0\\
0&-2b&\gamma&0\\
2\tau&0&0&0
\end{pmatrix}.
\]

In the temporary \(r,s,u,v\) coordinates of Case 2, it is

\[
\begin{pmatrix}
-2a&0&2\eta&-2c\\
0&2e&0&2f\\
2\eta&0&0&0\\
-2c&2f&0&\gamma
\end{pmatrix}.
\]

These make the signs and factors of two in the product equations directly checkable.

### 1.6 Case 2 uses a coordinate change, not an additional symmetry

The substitution \(r=(x-y)/\sqrt2\), \(s=(x+y)/\sqrt2\) by itself is not being claimed as a symmetry of the original pair \(Q,H\). It is a passive coordinate change, under which the polynomial is rewritten as

\[
Q=2rsu+(s^2-r^2)v.
\]

If the frame matrix becomes \(PT\), the gradient matrix becomes \(GT\), and the coefficient matrix becomes \(T^TST\). This preserves symmetry of the coefficient matrix without asserting invariance of the original polynomial formula. The actual final normalization of the placement is the allowed circle rotation \(h_{\pi/4}\), which also rotates the second weight plane. I inserted this distinction explicitly.

### 1.7 The outside set is necessarily nonempty

The normal form has \(v_i=\omega\) on \(K_1\cup K_2\). If this union were all indices, the zero-mean condition would force \(v=0\), contradicting its unit norm. Therefore the outside set is nonempty, and \(k_1+k_2\le N-1\). I stated this consequence and its short proof.

I also clarified that the normal-form proposition gives necessary conditions; it does not claim that every frame satisfying the displayed conditions is exceptional. Some of the original script’s test frames satisfy the normal-form equations without being asserted to be exceptional, and that is legitimate for testing the evaluation formula.

### 1.8 Minor editorial changes

I replaced ambiguous \(\Delta x^2\) notation by \((\Delta x)^2\), defined the adjoint and symmetric-operator notation, corrected “corollarys” to “corollaries,” used the published Arabic section numbers for Kuperberg, and qualified the abstract’s reference to the predecessor construction as “previously proposed.” The verification disclaimer remains. A revision note documents the nature of the changes.

## 2. Audit of the all-dimensional dimension argument

### 2.1 Quadratic evaluations and signed orbit interpolation

Lemma 2.3 has the right normalization. For a centered symmetric matrix \(B\), the edge values form the zero-diagonal symmetric matrix

\[
W_{ij}=(B_{ii}+B_{jj})/2-B_{ij}.
\]

Then \(\Pi W\Pi=-B\). Thus evaluations on the \(n=d(d+1)/2\) positive roots determine a symmetric operator on \(V\). The equal-dimension argument is valid. The expanded verifier additionally tests the explicit inverse on arbitrary symbolic edge labels for two sizes.

The orbit-separating map \(\Psi\) is adequate, including when one or both complex coordinates vanish. An orbit contains its negative exactly on the two pure weight planes; a nonzero fixed component prevents this. Consequently, outside the exceptional set, the signed representatives give \(2m\) distinct orbit invariants.

Lagrange interpolation on these \(2m\) points has degree at most \(2m-1\). Composing with the degree-at-most-three map \(\Psi\), then taking the odd part, yields the stated bound

\[
3(2m-1)\le6n-3=3d(d+1)-3.
\]

Taking the odd part does not spoil circle invariance, because negation commutes with the action. I found no degree-bound or interpolation-surjectivity error.

### 2.2 Spectral fibers: why no additional continuous freedom survives

Fix \(B=S+cc^T\), with \(S=A^TD_\lambda A\). Each of the two double eigenspaces of \(S\) has a nonzero vector perpendicular to \(c\), so both \(\lambda_1\) and \(\lambda_2\) remain eigenvalues of \(B\). There are only finitely many choices of the pair for fixed \(B\).

For each such pair, the ordered eigenvalue bounds place every eigenvalue of \(B\) in the interval \([s,s+1)\) associated with its eigenvalue of \(S\). The intervals are disjoint because of the gaps in \(D_\lambda\) and the bound \(\|c\|<1\). A simple cluster has one eigenvalue and a double cluster has two, at least one unchanged. A repeated value of \(B\) therefore must be the unchanged center of a double cluster.

The secular identity

\[
\frac{\det(tI-S)}{\det(tI-B)}
=1+\sum_{\beta\in\operatorname{spec}B}
\frac{\|P_\beta c\|^2}{t-\beta}
\]

then determines all squared projection lengths. At a double eigenvalue of \(B\), the same double factor is present in the numerator, so there is no pole and the corresponding component of \(c\) vanishes. This is the essential step: it rules out rotating a nonzero component of \(c\) continuously in a two-dimensional eigenspace. On the remaining one-dimensional eigenspaces, only sign choices remain.

Thus \(c\) has finitely many possibilities, \(S=B-cc^T\) is then fixed, and the remaining placement freedom is a coset of

\[
O(2)\times O(2)\times\{\pm1\}^{d-4}.
\]

Its dimension is two. The many extra fixed-coordinate directions in high dimension contribute only discrete signs because their auxiliary eigenvalues are all distinct. This is why the fiber bound does not grow with \(d\). The proof covers \(c=0\) as well as nonzero \(c\).

### 2.3 The two dimension inequalities are used in the correct directions

On a signed piece \(Z\), let \(m\) be the number of classes and \(\ell=\dim L_Z\). The auxiliary spectral map has domain dimension \(\dim Z+2+\ell\), image dimension at most \(m\), and fibers of dimension at most two. Hence

\[
\dim Z+\ell\le m.
\]

The actual incidence set uses unrestricted \(c\in L_Z\), not \(\|c\|<1\). Its dimension is

\[
\dim\mathcal I_Z=\dim Z+\ell+p-m\le p,
\qquad p=\dim\mathcal P_d.
\]

This alone would not prove avoidance. The strict improvement comes from the circle contained in every nonempty fiber over a polynomial. That circle acts freely on the *placement* because \(A\) is surjective and the action on \(E_1\) is faithful. Individual roots may have stabilizers or may lie in \(E_0\); that does not invalidate freeness on placements.

The lower fiber-dimension inequality gives image dimension at most \(p-1\). The finite union over signed pieces is still smaller-dimensional and semialgebraic, so its complement is dense. I found neither a reversed inequality nor an illicit restriction on the linear functional in this step. Hardt-type semialgebraic triviality is sufficient; a manifold or transversality assumption on the pieces is unnecessary.

## 3. Audit of exceptional cubic solutions

### 3.1 The reduction to four coordinate columns is valid

Although the simplex frame has \(d\) columns, \(Q\) depends only on its first four coordinates. Those columns are centered and orthonormal. The gradient columns also have zero mean: mixed products sum to zero and \(\sum x_i^2=\sum y_i^2=1\). Therefore the cubic residual is exactly

\[
\mathcal R_Q(A)=\frac1{2\sqrt2}(PG^T-GP^T).
\]

Exactness implies \(G=PS\) with \(S\) symmetric, without dropping any contribution from the remaining coordinates. The global cubic annihilation identity likewise follows from the full tight-frame identity and the infinitesimal circle invariance of \(Q\).

### 3.2 Case 1: a root in the first weight plane

Aligning the root makes \(x=(e_0-e_1)/\sqrt2\). All other coordinate columns have equal entries at these indices. The first coefficient matrix in Section 1.5 gives precisely the three product equations used in the original proof.

The contradiction argument for \(\tau\ne0\) is valid. It forces \(W=0\), supplies another index with \(y_i=\tau\), and comparison of the square and product equations gives \(b=\tau\), \(a=0\). All entries of \(y\) then lie in \(\{0,\tau,2\tau\}\), contradicting zero mean.

Once \(\tau=0\), the special-index equation gives \(\gamma U=1/2\), so division by \(\gamma\) and the deduction \(b=0\) are justified. The disjoint-support, constant-magnitude, and balanced-sign conclusions follow from the displayed product equations, zero mean, and unit norm. I found no omitted zero-denominator case.

### 3.3 Case 2: a root in the second weight plane

The first reduction gives \(X^2=Y^2\). Simultaneous complex conjugation legitimately permits \(X=Y\): it preserves cubic exactness and reverses both factors of the quintic pairing, leaving the pairing unchanged.

In the temporary \(r,s,u,v\) coordinates, the product relations and their symmetric coefficients are correct. Under \(\eta\ne0\), the proof handles both \(c\ne0\) and \(c=0\), and both branches force all entries of \(s\) into \(\{0,\eta,2\eta\}\). This is a valid zero-mean contradiction. After \(\eta=0\), the disjoint-support equations force \(c=f=0\), and then \(\gamma\ne0\). The remaining norm and sign conclusions follow.

Finally, the *allowed* rotation \(h_{\pi/4}\) produces \(z'_1=r+is\), \(z'_2=-v+iu\). This correctly yields the required normal form, with the new imaginary column vanishing on both supports, hence \(\omega=0\). The revised discussion separates this rotation from the earlier passive coordinate substitution.

### 3.4 The quintic pairing has the claimed nonzero value

I checked the three nonzero edge-group contributions and their signs. The sum on cross-support edges is \(\gamma^3/(2\sqrt2)\). The two groups of support-to-outside edges contribute

\[
\frac{k_1^{-2}+k_2^{-2}}{2\sqrt2\,\gamma}
\sum_{j\in O}v_j(\omega-v_j).
\]

Zero mean and unit norm of \(v\) give the outside moment \(-1\), including when \(\omega=0\) and when outside entries vary. Using \(\gamma^2=k_1^{-1}+k_2^{-1}\) gives

\[
\Lambda_A(H)=\frac1{\sqrt2\,k_1k_2\gamma}\ne0.
\]

The sign depends on \(\gamma\); the argument needs nonvanishing, not a globally positive sign. Both signs are tested in the original and expanded scripts.

## 4. Completion and geometric consequence

The compactness argument in Proposition 7.1 is sound. A sequence of exceptional zeros for nonzero \(\varepsilon_j\to0\) would, by cubic annihilation, have \(\Lambda_{A_j}(H)=0\). A subsequential limit would be an exceptional cubic zero with vanishing quintic pairing, contradicting Corollary 6.3.

For fixed admissible \(\varepsilon\), compactness gives a positive residual gap on the exceptional set. The uniform finite-dimensional bound on perturbation residuals makes this property open in coefficient space. Density of nonexceptional-avoiding polynomials then supplies a perturbation in any requested sufficiently small ball. This does not require the nonexceptional bad set to be closed, nor a uniform residual gap as \(\varepsilon\to0\).

The constant-width construction in Lemma 8.1 has the correct baseline \(1/2\). The homogeneous extension has a positive-semidefinite Hessian off the origin for small \(\tau\), and the separate line-through-the-origin argument handles the point where the Hessian is not defined. The supporting-gradient argument establishes that the constructed convex set has exactly the desired support function.

If this body were contained in a translated image of \(C_d\), the two opposite support inequalities would have equal sums, namely one. Both must be equalities, so \(\tau f\) agrees with a linear functional on all rotated roots. This contradicts the functional theorem. Constant width one also gives diameter exactly one, as proved in the manuscript. No appeal to a completion theorem for arbitrary diameter-one sets is needed for this counterexample direction.

## 5. Computational and document checks

The original script ran successfully under **Python 3.13.5 / SymPy 1.14.0**, with **1,322 exact zero assertions**. Its own limitations are accurately stated: it does not certify the classification or dimension argument.

The expanded self-contained script retains the original checks and adds passive-coordinate-change identities, symbolic quadratic-evaluation inverses, six deterministic full simplex-frame tests in dimensions 4, 5, and 6, eight selected spectral examples in dimensions 4 and 6, and extra normal-form frames with zero and positive/negative nonzero \(\omega\), including nonconstant outside entries. There are 29 normal-form frame tests in total. Selected tests explicitly check the presence of an \(E_1\) or \(E_2\) root when the construction guarantees one. It completed with **4,748 exact zero assertions**, with no floating-point arithmetic or random sampling.

These finite examples are regression tests, not an exhaustive normal-form classification, a proof of the spectral fiber bound, or a certification of the theorem. The determinant tests in particular check signs and selected multiplicity cases; the general finiteness-of-fibers argument is the mathematical one in Section 2.2 above.

The revised LaTeX compiled in three passes. The final compilation has no warnings, undefined references, or overfull/underfull box messages. The resulting **16-page PDF** was rendered and all pages inspected for layout problems. The source diff records the changes; the original uploads were not overwritten.

## 6. External-source checks and scope

Kuperberg, *Circumscribing constant-width bodies with polytopes*, New York Journal of Mathematics 5 (1999), 91–100: Conjecture 1 and Proposition 1 support the conjecture identification and functional formulation; Section 5 contains the low-dimensional quadratic-evaluation discussion cited by the manuscript. I checked the published PDF, including its displayed statement and normalization.

Chen, *Note on Lebesgue’s universal cover problem*, arXiv:2607.27227v2 (August 6, 2026): the cited preprint exists; its introduction distinguishes the regular-simplex question from the maximal-facet-count conjecture, and Appendix A, Proposition 13 gives the attributed above-threshold dimension-count argument. Its Lemma 8 has a normalization mismatch in its proof as written: the displayed baseline one yields width two, whereas the statement says width one. This does **not** carry into the uploaded manuscript, whose self-contained Lemma 8.1 correctly uses baseline one-half; I retained that correct proof.

The predecessor repository manuscript at `simplex4/math-results`, file `makeev-code/a4/makeev_d4_d5.tex`, exists and contains the stated cubic–quintic seed and dimension-four/five claims. Checking the repository text and its attribution is not an independent verification of its certificates, which were not needed or run here. The current manuscript reproduces the algebra it uses and does not logically depend on those certificates.

The cited Bochnak–Coste–Roy, Hardt, and Hardt–Lambrechts–Turchin–Volić sources provide appropriate semialgebraic background. The needed dimension inequalities are also explicitly derived from semialgebraic triviality in the manuscript. No outside result was silently substituted for a missing step in the proposed proof.

## 7. What remains unprovided, and what that means

The revision still gives no explicit coefficients for the auxiliary perturbation \(R\), no numerical value of \(\varepsilon_0\), and no proof that the unmodified quintic \(Q+\varepsilon H\) works in every dimension. These are accurately disclosed limitations of the construction, not contradictions of its stated existence theorem. I did not attempt a comprehensive novelty determination or claim independent human verification.

**Bottom line:** the edits repair wording and make several delicate deductions more explicit; I found no reason in this review to retract or weaken the all-dimensional existence theorem. Independent expert review should concentrate on Lemma 4.4, the stratum/incidence count in Lemma 4.5 and Proposition 4.6, and the two-case proof of Proposition 6.1, rather than treating the passing algebra script as a substitute for those arguments.
