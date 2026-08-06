"""Stand-alone verifier for the AME(12,5) paper.

Uses ONLY the Python standard library and exact integer arithmetic.
Checks every item needed for Theorem 1 and Corollary 7:

  [T1] Theta is symmetric with entries in F_5, and ALL 462 balanced-cut
       blocks Theta[S,S^c] (representative: the side containing
       coordinate 1, rows/columns sorted) are nonsingular over F_5.
       Det-residue checksum must be (89,154,131,88) for residues 1..4.
  [C1] A Abar^T = -I_6 over F_25 (36 Hermitian inner products), hence
       C = rowspace[I|A] is Hermitian self-dual.
  [C2] A is superregular: all 923 square minors are nonzero, with
       counts (36,225,400,225,36,1) by size; hence C is MDS with d = 7.
  [C3] Corollary 7 ingredients: with P the puncturing and B the
       shortening of C at coordinate 1, P is [11,6], B is [11,5],
       P is contained in B^perp_h and dim B^perp_h = 6, so B^perp_h = P.
       (MDS propagation then gives d(P)=6, d(B)=7 exactly, and the
       Hermitian construction yields the pure [[11,1,6]]_5 code.)

No floating point, no GPU, no third-party package, no enumeration of
25^6 codewords is used. Expected output ends with ALL CHECKS PASSED.
"""
import itertools

P5 = 5

THETA = [
    [0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 3],
    [0, 0, 0, 3, 2, 4, 2, 1, 4, 4, 4, 0],
    [0, 0, 0, 2, 1, 0, 0, 4, 1, 1, 0, 2],
    [0, 3, 2, 0, 2, 2, 3, 2, 4, 4, 3, 2],
    [0, 2, 1, 2, 0, 0, 3, 4, 1, 3, 3, 3],
    [0, 4, 0, 2, 0, 0, 3, 4, 4, 4, 3, 1],
    [4, 2, 0, 3, 3, 3, 0, 3, 0, 2, 3, 1],
    [4, 1, 4, 2, 4, 4, 3, 0, 1, 4, 0, 4],
    [4, 4, 1, 4, 1, 4, 0, 1, 0, 4, 4, 4],
    [4, 4, 1, 4, 3, 4, 2, 4, 4, 0, 3, 1],
    [4, 4, 0, 3, 3, 3, 3, 0, 4, 3, 0, 2],
    [3, 0, 2, 2, 3, 1, 1, 4, 4, 1, 2, 0],
]

A = [
    [1, 1, 1, 1, 1, 2],
    [23, 1, 2, 3, 5, 3],
    [19, 7, 1, 22, 17, 12],
    [9, 15, 11, 2, 18, 15],
    [1, 4, 22, 9, 2, 1],
    [6, 7, 18, 15, 10, 2],
]


# ---------------- F_25 = F_5[alpha], alpha^2 = 4*alpha + 3 -------------
def fmul(x, y):
    a, b = x % 5, x // 5
    c, d = y % 5, y // 5
    # (a + b alpha)(c + d alpha) = ac + (ad+bc) alpha + bd alpha^2
    lo = (a * c + 3 * b * d) % 5
    hi = (a * d + b * c + 4 * b * d) % 5
    return lo + 5 * hi


def fadd(x, y):
    return (x % 5 + y % 5) % 5 + 5 * ((x // 5 + y // 5) % 5)


def fneg(x):
    return (-x % 5) % 5 + 5 * ((-(x // 5)) % 5)


def fpow(x, e):
    r = 1
    while e:
        if e & 1:
            r = fmul(r, x)
        x = fmul(x, x)
        e >>= 1
    return r


def finv(x):
    assert x != 0
    return fpow(x, 23)          # x^(25-2)


def conj(x):
    return fpow(x, 5)           # Frobenius


def det25(M):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = fneg(d)
        d = fmul(d, M[c][c])
        iv = finv(M[c][c])
        for r in range(c + 1, n):
            if M[r][c]:
                f = fmul(M[r][c], iv)
                for cc in range(c, n):
                    M[r][cc] = fadd(M[r][cc], fneg(fmul(f, M[c][cc])))
    return d


def det5(M):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % 5), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = (-d) % 5
        d = d * M[c][c] % 5
        iv = pow(M[c][c], 3, 5)
        for r in range(c + 1, n):
            f = M[r][c] * iv % 5
            if f:
                for cc in range(c, n):
                    M[r][cc] = (M[r][cc] - f * M[c][cc]) % 5
    return d


def rank25(rows):
    M = [row[:] for row in rows]
    n, m = len(M), len(M[0])
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        iv = finv(M[r][c])
        M[r] = [fmul(iv, x) for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [fadd(M[i][cc], fneg(fmul(f, M[r][cc])))
                        for cc in range(m)]
        r += 1
        if r == n:
            break
    return r


def main():
    # [T1] graph-state certificate
    assert all(THETA[i][j] == THETA[j][i] for i in range(12)
               for j in range(12)), "Theta not symmetric"
    assert all(0 <= x < 5 for row in THETA for x in row)
    hist = [0] * 5
    for S in itertools.combinations(range(1, 12), 5):
        rows = (0,) + S
        cols = [j for j in range(12) if j not in rows]
        hist[det5([[THETA[i][j] for j in cols] for i in rows])] += 1
    assert hist[0] == 0, "singular balanced cut found"
    assert hist[1:] == [89, 154, 131, 88], f"checksum mismatch: {hist[1:]}"
    print("[T1] 462 balanced cuts nonsingular over F_5; "
          "det checksum (89,154,131,88) OK")

    # [C1] Hermitian unitarity / self-duality
    negI = [[fneg(1) if i == j else 0 for j in range(6)] for i in range(6)]
    prod = [[0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(6):
            s = 0
            for t in range(6):
                s = fadd(s, fmul(A[i][t], conj(A[j][t])))
            prod[i][j] = s
    assert prod == negI, "A Abar^T != -I"
    print("[C1] A Abar^T = -I_6 over F_25 (36 products) OK "
          "-> C = rowspace[I|A] is Hermitian self-dual")

    # [C2] superregularity
    counts = []
    for s in range(1, 7):
        c = 0
        for rs in itertools.combinations(range(6), s):
            for cs in itertools.combinations(range(6), s):
                assert det25([[A[r][col] for col in cs] for r in rs]) != 0, \
                    f"zero minor at rows {rs} cols {cs}"
                c += 1
        counts.append(c)
    assert counts == [36, 225, 400, 225, 36, 1], counts
    print("[C2] all 923 minors of A nonzero, counts "
          "(36,225,400,225,36,1) OK -> C is MDS, d = 7")

    # [C3] Corollary 7 ingredients
    G = [[(1 if i == j else 0) for j in range(6)] + A[i] for i in range(6)]
    Prows = [row[1:] for row in G]                 # puncture coordinate 1
    Brows = [row[1:] for row in G[1:]]             # shorten coordinate 1
    assert rank25(Prows) == 6 and rank25(Brows) == 5
    for p in Prows:
        for b in Brows:
            s = 0
            for t in range(11):
                s = fadd(s, fmul(p[t], conj(b[t])))
            assert s == 0, "P not contained in B^perp_h"
    # dim B^perp_h = 11 - dim B = 6 = dim P, hence equality
    print("[C3] P = [11,6], B = [11,5], P subset of B^perp_h, "
          "dim B^perp_h = 6 = dim P -> B^perp_h = P; MDS propagation "
          "gives d(P)=6, d(B)=7 -> pure [[11,1,6]]_5 exists")

    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
