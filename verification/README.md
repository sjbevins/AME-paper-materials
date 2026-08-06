# Portable exact verification

The proof layer is intentionally separate from the discovery programs.

Core standard-library checks:

```bash
python3 verify_ame_12_5.py
python3 verify_ame_18_11.py
python3 verify_ame_18_13.py
```

The first script proves AME(12,5), the Hermitian self-dual MDS certificate, and the structural [[11,1,6]]_5 corollary. The second proves AME(18,11). The third proves both AME(18,13) and all eighteen AME(17,13) deletions.

For the q=11 odd-party census, use either:

```bash
python3 derive_17_11.py            # NumPy; dual rank implementations for the primary deletion
python3 verify_17_11_numba.py      # NumPy + Numba; all 18 deletions
```

The corrected template classifier requires a C++17 compiler and Boost's header-only multiprecision library (only `boost/multiprecision/cpp_int.hpp` is included; any maintained Boost release provides it):

```bash
g++ -O3 -std=c++17 template_screen_generic.cpp -o template_screen_generic
./template_screen_generic 11
./template_screen_generic 13
```

`check_consistency.py` (standard library only) ties the packaged data to the manuscript: it confirms that the plain-text and NumPy matrices agree, that they equal the matrices embedded in the verifiers, that `constructions.json` reproduces the embedded kernels, that the 17-vertex arrays are principal deletions of the 18-vertex ones, and it recomputes all 462 + 2 x 24,310 balanced-cut certificates from the packaged text files. This replaces a file-hash manifest: agreement of independent recomputations is the integrity check.

`verify_all.sh` runs the consistency check, the core exact scripts, the fast q=11 all-deletion census when NumPy/Numba are available, and the corrected census when a suitable C++ compiler is available. It compares regenerated output against the included expected files. Discovery hit lists are in `../evidence/`; they are provenance only and are not used as mathematical certificates.
