# Makeev dimension-5 extension

The natural dimension-5 candidate is the same polynomial used in dimension 4,
with the fifth real coordinate ignored:

\[
f_\varepsilon(z_1,z_2,x)=\operatorname{Re}(z_1^2\overline z_2)
+\varepsilon |z_1|^2\operatorname{Im}(z_1^2\overline z_2).
\]

Current status: the proof has been reduced to an exact rational unit-ideal
calculation, but that calculation has **not** been run in this environment.
Do not call the dimension-5 theorem unconditional until Magma prints the
success line below.

Files:

- `makeev_d5_note.tex` — derivation, compactness argument, and exact algebraic
  reduction.
- `makeev_d5_certificate.m` — smaller dimension-5-specific Gröbner-basis test.
- `makeev_general_certificate.m` — quadratic/cubic test parameterized by `d`;
  change the first `d := 5;` assignment to test a higher fixed dimension.
- `makeev_higher_dim_numerics.py` — fast nonrigorous reconnaissance.

Run the dimension-5 certificate with:

```sh
magma makeev_d5_certificate.m
```

The conclusive output is:

```text
Certificate verified: the saturated ideal is the unit ideal.
```

If the computation exhausts memory or runs for too long, try the general
quadratic/cubic formulation:

```sh
magma makeev_general_certificate.m
```

That formulation introduces more variables but lowers the maximum degree from
7 to 3.  A non-unit ideal does not by itself produce a real counterexample to
transversality; it may consist only of complex solutions.  In that event,
compute a real-root classification or send back the Gröbner basis/dimension
output for the next reduction.

For numerical searches (NumPy and SciPy required):

```sh
python makeev_higher_dim_numerics.py --d 5 --starts 100
python makeev_higher_dim_numerics.py --d 6 --starts 100 --seed 6
```

Numerical nonvanishing is evidence only.  The exact certificate is what turns
the compactness argument into a proof.

