"""Packaged-data consistency check for the five-AME-results package.

Standard library only.  Package integrity is established mathematically,
not by file hashes: every packaged copy of every matrix is compared
against the other copies and against the matrices embedded in the
stand-alone verifiers, and the balanced-cut certificates are recomputed
directly from the packaged plain-text arrays.

Checks:

  [D1] For each of the five graph arrays, the plain-text matrix in
       data/, the NumPy array in data/, and the NumPy array beside the
       verifiers agree entry by entry (the .npy files are parsed with a
       minimal reader; NumPy itself is not required).
  [D2] Every array is symmetric with entries in the correct prime field.
  [D3] The 12- and 18-vertex arrays equal the matrices embedded in
       verify_ame_12_5.py and verify_ame_18_11.py (extracted from their
       sources as literals, without executing them).  Following the
       manuscript convention, the comparison is off-diagonal: the
       packaged arrays retain the diagonal quadratic phase terms from
       local-Clifford reduction, the embedded copies may zero them, and
       no cut block uses the diagonal.
  [D4] data/constructions.json agrees with the verifiers: A_5 equals the
       embedded systematic block, and the q=11 and q=13 kernels
       reproduce the embedded group-circulant data through
       A[x,y] = a(y-x) on Z_3 x Z_3 in the recorded coordinate order.
  [D5] Each 17-vertex array is a principal single-vertex deletion of its
       18-vertex parent.
  [D6] Recomputed certificates from the packaged text arrays alone:
       all 462 balanced cuts of the 12-vertex array are nonsingular over
       F_5 with determinant-residue counts (89,154,131,88), and all
       24,310 balanced cuts of each 18-vertex array are nonsingular over
       F_11 and F_13 respectively.

Expected output ends with ALL CONSISTENCY CHECKS PASSED.
"""
import ast
import itertools
import json
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

FAIL = []


def check(label, ok):
    print(("PASS" if ok else "FAIL") + "  " + label)
    if not ok:
        FAIL.append(label)


# ------------------------------------------------------------ file readers
def read_txt_matrix(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append([int(t) for t in line.split()])
    return rows


def read_npy_int_matrix(path):
    """Minimal .npy reader for little-endian integer 2-D arrays."""
    with open(path, "rb") as f:
        magic = f.read(6)
        assert magic == b"\x93NUMPY", "not a .npy file: %s" % path
        major, _minor = f.read(1)[0], f.read(1)[0]
        if major == 1:
            (hlen,) = struct.unpack("<H", f.read(2))
        else:
            (hlen,) = struct.unpack("<I", f.read(4))
        header = ast.literal_eval(f.read(hlen).decode("latin1"))
        descr, fortran, shape = (header["descr"], header["fortran_order"],
                                 header["shape"])
        assert not fortran and len(shape) == 2
        fmt = {"<i8": "q", "<i4": "i", "<i2": "h", "|i1": "b"}[descr]
        n = shape[0] * shape[1]
        flat = struct.unpack("<%d%s" % (n, fmt), f.read())
    return [list(flat[i * shape[1]:(i + 1) * shape[1]])
            for i in range(shape[0])]


def source_literals(path, names):
    """Extract top-level literal assignments from a script, unexecuted."""
    tree = ast.parse(open(path, encoding="utf8").read())
    out = {}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign)
                and isinstance(node.targets[0], ast.Name)):
            name = node.targets[0].id
            if name in names and name not in out:
                try:
                    out[name] = ast.literal_eval(node.value)
                except ValueError:
                    pass
    missing = set(names) - set(out)
    assert not missing, "literals not found in %s: %s" % (path, missing)
    return out


# ------------------------------------------------------- exact linear algebra
def det_mod(rows, p, inv):
    """Determinant of a square integer matrix modulo the prime p."""
    m = [row[:] for row in rows]
    n = len(m)
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if m[r][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            det = -det
        det = det * m[c][c] % p
        ipiv = inv[m[c][c] % p]
        for r in range(c + 1, n):
            f = m[r][c] * ipiv % p
            if f:
                m[r] = [(m[r][k] - f * m[c][k]) % p for k in range(n)]
    return det % p


def balanced_cut_census(theta, p):
    """Nonsingularity of Theta[S,S^c] for all balanced cuts containing
    coordinate 0; returns (singular_count, residue_histogram)."""
    n = len(theta)
    k = n // 2
    inv = [0] + [pow(x, p - 2, p) for x in range(1, p)]
    hist = [0] * p
    singular = 0
    coords = list(range(n))
    for s in itertools.combinations(range(1, n), k - 1):
        S = (0,) + s
        inS = set(S)
        Sc = [j for j in coords if j not in inS]
        block = [[theta[i][j] for j in Sc] for i in S]
        d = det_mod(block, p, inv)
        hist[d] += 1
        if d == 0:
            singular += 1
    return singular, hist


# ------------------------------------------------------------------- load all
NAMES = {
    "theta_ame125": 5,
    "theta_ame1811": 11,
    "theta_ame1813": 13,
    "theta_ame1711": 11,
    "theta_ame1713": 13,
}

txt, npy_data, npy_verif = {}, {}, {}
for name in NAMES:
    txt[name] = read_txt_matrix(os.path.join(DATA, name + ".txt"))
    npy_data[name] = read_npy_int_matrix(os.path.join(DATA, name + ".npy"))
    npy_verif[name] = read_npy_int_matrix(os.path.join(HERE, name + ".npy"))

# [D1] three copies of each array agree
for name in sorted(NAMES):
    check("[D1] %s: text, data/.npy, verification/.npy identical" % name,
          txt[name] == npy_data[name] == npy_verif[name])

# [D2] symmetry and field range
for name, p in sorted(NAMES.items()):
    m = txt[name]
    sym = all(m[i][j] == m[j][i] for i in range(len(m))
              for j in range(len(m)))
    rng = all(0 <= x < p for row in m for x in row)
    check("[D2] %s: symmetric with entries in F_%d" % (name, p), sym and rng)

# [D3] packaged arrays equal the matrices embedded in the verifiers
v125 = source_literals(os.path.join(HERE, "verify_ame_12_5.py"),
                       {"THETA", "A"})
v1811 = source_literals(os.path.join(HERE, "verify_ame_18_11.py"),
                        {"THETA", "A"})
v1813 = source_literals(os.path.join(HERE, "verify_ame_18_13.py"),
                        {"A_GEN"})
def off_diagonal_equal(m1, m2, p):
    n = len(m1)
    return len(m2) == n and all(
        m1[i][j] % p == m2[i][j] % p
        for i in range(n) for j in range(n) if i != j)


check("[D3] theta_ame125 equals THETA embedded in verify_ame_12_5.py "
      "(off-diagonal)",
      off_diagonal_equal(txt["theta_ame125"], v125["THETA"], 5))
check("[D3] theta_ame1811 equals THETA embedded in verify_ame_18_11.py "
      "(off-diagonal)",
      off_diagonal_equal(txt["theta_ame1811"], v1811["THETA"], 11))

# [D4] constructions.json agrees with the embedded search results
cons = json.load(open(os.path.join(DATA, "constructions.json")))
order = [tuple(g) for g in cons["coordinate_order_Z3xZ3"]]
check("[D4] constructions.json A_5 equals A embedded in verify_ame_12_5.py",
      cons["AME_12_5"]["A_5_integer_encoding"] == v125["A"])

kernel11 = cons["AME_18_11_and_17_11"]["kernel_integer_encoding"]
A11 = [[kernel11[(y[0] - x[0]) % 3][(y[1] - x[1]) % 3]
        for y in order] for x in order]
check("[D4] q=11 kernel + group rule reproduce A embedded in "
      "verify_ame_18_11.py", A11 == v1811["A"])

kernel13 = cons["AME_18_13_and_17_13"]["kernel_integer_encoding"]
agen13 = [[c0 + 13 * c1 for (c0, c1) in row] for row in v1813["A_GEN"]]
check("[D4] q=13 kernel equals A_GEN embedded in verify_ame_18_13.py",
      kernel13 == agen13)

# [D5] the 17-vertex arrays are principal deletions of their parents
for child, parent in (("theta_ame1711", "theta_ame1811"),
                      ("theta_ame1713", "theta_ame1813")):
    big, small = txt[parent], txt[child]
    matches = []
    for v in range(len(big)):
        keep = [i for i in range(len(big)) if i != v]
        if [[big[i][j] for j in keep] for i in keep] == small:
            matches.append(v)
    check("[D5] %s is the deletion of vertex %s of %s"
          % (child, matches if matches else "<none>", parent),
          len(matches) >= 1)

# [D6] balanced-cut certificates recomputed from the packaged text arrays
sing, hist = balanced_cut_census(txt["theta_ame125"], 5)
check("[D6] theta_ame125: all 462 balanced cuts nonsingular over F_5",
      sing == 0 and sum(hist) == 462)
check("[D6] theta_ame125: determinant-residue counts equal %s"
      % (cons["AME_12_5"]["determinant_residue_checksum"],),
      hist[1:] == cons["AME_12_5"]["determinant_residue_checksum"])

for name, p in (("theta_ame1811", 11), ("theta_ame1813", 13)):
    sing, hist = balanced_cut_census(txt[name], p)
    check("[D6] %s: all 24,310 balanced cuts nonsingular over F_%d"
          % (name, p), sing == 0 and sum(hist) == 24310)

print()
if FAIL:
    print("CONSISTENCY CHECK FAILURES: %d" % len(FAIL))
    raise SystemExit(1)
print("ALL CONSISTENCY CHECKS PASSED")
