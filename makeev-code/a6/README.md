# Odd polynomial counterexamples for the Ad root system

## Status and scope

This is an AI-generated verification draft prepared on September 4, 2026,
from a research conversation with ChatGPT. It presents a proposed proof
for independent human verification and has not been independently
verified or peer reviewed.

The manuscript proposes a counterexample in every dimension d >= 4,
using unit-normalized Ad roots. The polynomial has degree at most
3d(d+1)-3. The final small invariant perturbation is existential: the
manuscript does not supply its coefficients or prove that the original
cubic–quintic seed works unchanged in all dimensions.

## Files

- `makeev_all_dimensions.tex`: standalone LaTeX source; bibliography embedded.
- `makeev_all_dimensions.pdf`: compiled manuscript.
- `verify_algebra.py`: optional exact symbolic regression checks.
- `verification_output.txt`: output from the supplied script in this environment.
- `README.md`: this file.

## Building the paper

Use a reasonably complete TeX Live or MiKTeX installation. No external
figures, bibliography files, custom fonts, shell escape, or network access
are needed. The source uses the AMS article class and standard packages.

From the directory containing the source:

```sh
pdflatex -interaction=nonstopmode -halt-on-error makeev_all_dimensions.tex
pdflatex -interaction=nonstopmode -halt-on-error makeev_all_dimensions.tex
```

An additional pass may be needed after changing headings, references, or
the table of contents. Alternatively:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error makeev_all_dimensions.tex
```

The supplied PDF was compiled with pdfTeX 1.40.26. The final compiler log
contained no unresolved-reference, overfull-box, underfull-box, or package
warnings. Rendered pages were also checked for layout and clipping.

## Running the optional algebra checks

The script needs Python 3.10 or later and SymPy. It was run with Python
3.13.5 and SymPy 1.14.0. To run it:

```sh
python verify_algebra.py
```

If SymPy is not installed, it can be installed into the chosen Python
environment with `python -m pip install sympy`.

The script uses exact symbolic arithmetic, with no random trials or
floating-point tolerances. The supplied run passed 1,322 exact zero
assertions, including matrix entries. It checks the polynomial identities,
edge-group sums, and ten explicitly specified normal-form frames.

These are regression checks, not a proof certificate. In particular, they
do not verify the spectral fiber bound, semialgebraic dimension arguments,
exhaustiveness of the necessary normal form, or the complete theorem.
No script output is used as a premise in the written proof.

## Suggested verification order

The proof has two independent branches. The generic branch uses the
quadratic evaluation isomorphism (Lemma 2.3), spectral fiber bound
(Lemma 4.4), signed-stratum inequality (Lemma 4.5), and incidence-dimension
argument (Proposition 4.6). The exceptional branch uses the symmetric
gradient equation (Lemma 5.3), necessary normal form (Proposition 6.1),
and nonzero edge sum (Lemma 6.2). Section 7 combines them by compactness
and density. Appendix A records the main normalization and parameter
issues worth checking.

In particular, the auxiliary vector in the spectral map is restricted to
a unit ball, but the linear-functional parameter in the agreement
incidence set is unrestricted. In the exceptional calculation, roots
have length one, so the quintic edge evaluation has the factor
1/(4 sqrt(2)).

## Attribution

Section 1.1 and the bibliography distinguish the inherited seed and
symmetry strategy from the additional arguments proposed here. They
cite Kuperberg's functional reformulation, the earlier AI-generated
four- and five-dimensional manuscript linked by VibeMathed, Chen's
paper and its appendix on finite-dimensional dimension counting,
and the sources for semialgebraic triviality and dimension theory.
The elementary graph and matrix identities used in the proof are
proved in the manuscript. The earlier dimension-specific
computer-algebra certificates are not used as proof premises.
