# Makeev counterexamples in dimensions 4 and 5

This bundle contains the unified paper and the exact algebraic certificate
sources described in it.

## Files

- `makeev_d4_d5.tex` — unified paper, including the complete direct
  dimension-4 transversality appendix and the full Python generator listing.
- `tensor_polynomials.py` — standard-library generator for every fixed
  dimension `d >= 4`.
- `makeev_d4_tensor_certificate.m` — generated 11-variable Magma test for
  dimension 4.
- `makeev_d5_tensor_certificate.m` — generated 21-variable Magma test for
  dimension 5; this is byte-for-byte the file already verified in Magma.

The older `makeev_general_certificate.m` is deliberately excluded: its
coordinate ideal is too slow in dimensions 4 and 5 and is not used by the
paper.

## Reproduce the generated certificates

```sh
python3 tensor_polynomials.py --dimension 4 \
  --emit-magma makeev_d4_tensor_certificate.new.m
python3 tensor_polynomials.py --dimension 5 \
  --emit-magma makeev_d5_tensor_certificate.new.m

cmp makeev_d4_tensor_certificate.m makeev_d4_tensor_certificate.new.m
cmp makeev_d5_tensor_certificate.m makeev_d5_tensor_certificate.new.m
```

The generator should report the following summaries when run without an
output option:

| dimension | linear equations | free variables | commutator equations | obstruction terms |
|---:|---:|---:|---:|---:|
| 4 | 9 | 11 | 21 | 73 |
| 5 | 14 | 21 | 55 | 163 |

## Run the exact Magma checks

```sh
magma makeev_d4_tensor_certificate.m
magma makeev_d5_tensor_certificate.m
```

The conclusive line is

```text
Certificate verified: tensor ideal is the unit ideal.
```

It is printed only after Magma's `IsProper` returns `false` and the assertion
that the ideal is the unit ideal passes. The dimension-5 file produced this
line in approximately 14 seconds in the verified run. The paper's
dimension-4 result is also proved directly in its first appendix, so it does
not depend on the separate dimension-4 Magma run.

## Check exact source identities

```sh
sha256sum tensor_polynomials.py \
  makeev_d4_tensor_certificate.m \
  makeev_d5_tensor_certificate.m
```

Expected digests:

```text
fe71e120af7f7c38806701e07042f88bff461c5a19738653d78a0fbe9a7677ff  tensor_polynomials.py
4ae9030dd0b57c9541c7e6861d9d87a956015e009bd5220c2903eb222018d8d7  makeev_d4_tensor_certificate.m
403d831600c2ac4a6d60ba768506ff031eb1d4817eadde04d80e40e3aa14d2e6  makeev_d5_tensor_certificate.m
```

## Generate a higher-dimensional test

For example:

```sh
python3 tensor_polynomials.py --dimension 6 \
  --emit-magma makeev_d6_tensor_certificate.m
```

The formulation is exact for every fixed `d >= 4`, but no claim is made that
Magma will finish quickly as `d` grows.

## Compile the paper

```sh
pdflatex makeev_d4_d5.tex
pdflatex makeev_d4_d5.tex
```

The second pass resolves cross-references.
