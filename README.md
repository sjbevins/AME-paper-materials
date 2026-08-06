# Symmetry-guided construction and exact certification of AME states

Reproducibility package for the paper

> **Symmetry-guided construction and exact certification of absolutely
> maximally entangled states**
> Samuel Bevins and Yunus Bidav
> (arXiv identifier to be added on announcement)

The paper constructs Hermitian self-dual MDS codes with parameters
[12,6,7]_25, [18,9,10]_121, and [18,9,10]_169 and, from them, certifies
the five absolutely maximally entangled states

| state | equivalent pure quantum code |
|---|---|
| AME(12,5)  | [[12,0,7]]_5   |
| AME(17,11) | [[17,0,9]]_11  |
| AME(18,11) | [[18,0,10]]_11 |
| AME(17,13) | [[17,0,9]]_13  |
| AME(18,13) | [[18,0,10]]_13 |

Every claim is re-established by stand-alone exact verifiers included
here; no output of the original search programs needs to be trusted.

## Quick start: verify everything

Requires Python 3.8+ (standard library only for the core checks).

```bash
cd verification
bash verify_all.sh
```

This runs, comparing each program's output against its packaged
expected output:

1. `check_consistency.py` — ties the packaged data files to the
   matrices printed in the paper and embedded in the verifiers, and
   recomputes all balanced-cut certificates from the packaged arrays
   (this package's integrity check; there is no file-hash manifest).
2. `verify_ame_12_5.py` — AME(12,5), the [12,6,7]_25 Hermitian
   self-dual MDS certificate, and the [[11,1,6]]_5 corollary.
3. `verify_ame_18_11.py` — AME(18,11) via the [18,9,10]_121 code.
4. `verify_ame_18_13.py` — AME(18,13), all eighteen AME(17,13)
   deletions, and the expanded-basis and symplectic cross-checks.
5. `verify_17_11_numba.py` — all eighteen AME(17,11) deletions
   (skipped automatically if NumPy/Numba are not installed).
6. The corrected abelian-template census (skipped automatically if a
   C++17 compiler with Boost's header-only multiprecision library is
   not available).

Each verifier is deliberately self-contained: it embeds its own field
arithmetic and its own copy of the construction data, so agreement
between the verifiers, the packaged data, and the paper is evidence
rather than shared code.

## Layout

- `paper/` — complete LaTeX source of the manuscript (REVTeX 4.2).
  Build with `pdflatex main && bibtex main && pdflatex main && pdflatex
  main`.  `figures/gen_figures.py` regenerates the four data-driven
  figures from `data/`; all other figures are hand-written TikZ.
- `data/` — the construction record (`constructions.json`), the five
  graph matrices in plain text and NumPy formats, and the archived
  Huber–Wyderka AME-table snapshot used in the comparison section.
- `verification/` — the stand-alone verifiers, the consistency check,
  expected outputs, and `verify_all.sh`.
- `evidence/` — partial hit lists and original scripts from the
  discovery runs; provenance only, not used as mathematical
  certificates.

## Requirements

| task | needs |
|---|---|
| core verification | Python 3.8+ (standard library only) |
| q=11 all-deletion census | NumPy, Numba |
| template census | C++17 compiler, Boost (header-only `cpp_int`) |
| build the paper | TeX Live or MiKTeX with REVTeX 4.2, TikZ, booktabs |
| regenerate figures | Python 3.8+ (standard library only) |

## License

MIT — see `LICENSE`.  If you use this package, please cite the paper
(see `CITATION.cff`).
