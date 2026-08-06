"""Stand-alone verifier for AME(18,11).  Python standard library only,
exact integer arithmetic throughout.

Field: F_121 = F_11(alpha), alpha^2 = -1; encoding x = a + 11*b for
a + b*alpha; conjugation (Frobenius x -> x^11) is a + b*alpha -> a - b*alpha.

Checks:
  [C1] A Abar^T = -I_9 over F_121 (81 Hermitian products)  =>  the code
       C = rowspace[I_9 | A] in F_121^18 is Hermitian self-dual.
  [C2] every square minor of A is nonzero (48,619 minors)  =>  C is MDS,
       [18,9,10]_121  =>  pure [[18,0,10]]_11 stabilizer code = AME(18,11).
  [T1] the derived graph matrix Theta (18x18 symmetric over F_11) has all
       24,310 balanced 9|9 cut blocks nonsingular over F_11 -- the
       graph-state certificate of the SAME state, checked independently.

Expected output ends with ALL CHECKS PASSED.
"""
import itertools

P = 11

A = [[92, 118, 119, 43, 88, 43, 119, 118, 92],
     [119, 92, 118, 43, 43, 88, 92, 119, 118],
     [118, 119, 92, 88, 43, 43, 118, 92, 119],
     [119, 118, 92, 92, 118, 119, 43, 88, 43],
     [92, 119, 118, 119, 92, 118, 43, 43, 88],
     [118, 92, 119, 118, 119, 92, 88, 43, 43],
     [43, 88, 43, 119, 118, 92, 92, 118, 119],
     [43, 43, 88, 92, 119, 118, 119, 92, 118],
     [88, 43, 43, 118, 92, 119, 118, 119, 92]]

THETA = [
    [10, 9, 9, 0, 5, 8, 0, 8, 5, 2, 7, 7, 0, 9, 0, 7, 7, 2],
    [9, 10, 9, 8, 0, 5, 5, 0, 8, 7, 2, 7, 0, 0, 9, 2, 7, 7],
    [9, 9, 10, 5, 8, 0, 8, 5, 0, 7, 7, 2, 9, 0, 0, 7, 2, 7],
    [0, 8, 5, 10, 9, 9, 0, 5, 8, 7, 7, 2, 2, 7, 7, 0, 9, 0],
    [5, 0, 8, 9, 10, 9, 8, 0, 5, 2, 7, 7, 7, 2, 7, 0, 0, 9],
    [8, 5, 0, 9, 9, 10, 5, 8, 0, 7, 2, 7, 7, 7, 2, 9, 0, 0],
    [0, 5, 8, 0, 8, 5, 10, 9, 9, 0, 9, 0, 7, 7, 2, 2, 7, 7],
    [8, 0, 5, 5, 0, 8, 9, 10, 9, 0, 0, 9, 2, 7, 7, 7, 2, 7],
    [5, 8, 0, 8, 5, 0, 9, 9, 10, 9, 0, 0, 7, 2, 7, 7, 7, 2],
    [2, 7, 7, 7, 2, 7, 0, 0, 9, 1, 2, 2, 0, 6, 3, 0, 3, 6],
    [7, 2, 7, 7, 7, 2, 9, 0, 0, 2, 1, 2, 3, 0, 6, 6, 0, 3],
    [7, 7, 2, 2, 7, 7, 0, 9, 0, 2, 2, 1, 6, 3, 0, 3, 6, 0],
    [0, 0, 9, 2, 7, 7, 7, 2, 7, 0, 3, 6, 1, 2, 2, 0, 6, 3],
    [9, 0, 0, 7, 2, 7, 7, 7, 2, 6, 0, 3, 2, 1, 2, 3, 0, 6],
    [0, 9, 0, 7, 7, 2, 2, 7, 7, 3, 6, 0, 2, 2, 1, 6, 3, 0],
    [7, 2, 7, 0, 0, 9, 2, 7, 7, 0, 6, 3, 0, 3, 6, 1, 2, 2],
    [7, 7, 2, 9, 0, 0, 7, 2, 7, 3, 0, 6, 6, 0, 3, 2, 1, 2],
    [2, 7, 7, 0, 9, 0, 7, 7, 2, 6, 3, 0, 3, 6, 0, 2, 2, 1],
]


def fmul(x, y):
    a, b = x % 11, x // 11
    c, d = y % 11, y // 11
    return (a * c - b * d) % 11 + 11 * ((a * d + b * c) % 11)


def fadd(x, y):
    return (x % 11 + y % 11) % 11 + 11 * ((x // 11 + y // 11) % 11)


def fneg(x):
    return (-x % 11) % 11 + 11 * ((-(x // 11)) % 11)


def conj(x):
    return x % 11 + 11 * ((-(x // 11)) % 11)


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
    return fpow(x, 119)


def det121(M):
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


def det11(M):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % 11), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = (-d) % 11
        d = d * M[c][c] % 11
        iv = pow(M[c][c], 9, 11)
        for r in range(c + 1, n):
            f = M[r][c] * iv % 11
            if f:
                for cc in range(c, n):
                    M[r][cc] = (M[r][cc] - f * M[c][cc]) % 11
    return d


def main():
    # sanity of the field implementation
    assert fmul(11, 11) == 10, "alpha^2 != -1"
    for x in (1, 5, 11, 37, 92, 118):
        assert fmul(x, finv(x)) == 1
        assert fmul(x, conj(x)) < 11, "norm not in F_11"
    # [C1]
    for i in range(9):
        for j in range(9):
            s = 0
            for t in range(9):
                s = fadd(s, fmul(A[i][t], conj(A[j][t])))
            want = fneg(1) if i == j else 0
            assert s == want, f"A Abar^T mismatch at ({i},{j})"
    print("[C1] A Abar^T = -I_9 over F_121 -> C Hermitian self-dual")
    # [C2]
    cnt = 0
    for size in range(1, 10):
        for rs in itertools.combinations(range(9), size):
            for cs in itertools.combinations(range(9), size):
                assert det121([[A[r][c] for c in cs] for r in rs]) != 0
                cnt += 1
    assert cnt == 48619
    print("[C2] all 48,619 minors of A nonzero -> [18,9,10]_121 MDS "
          "-> pure [[18,0,10]]_11")
    # [T1]
    assert all(THETA[i][j] == THETA[j][i] for i in range(18)
               for j in range(18))
    tot = 0
    for S in itertools.combinations(range(1, 18), 8):
        rows = (0,) + S
        cols = [j for j in range(18) if j not in rows]
        assert det11([[THETA[i][j] for j in cols] for i in rows]) != 0, \
            f"singular cut {rows}"
        tot += 1
    assert tot == 24310
    print("[T1] all 24,310 balanced cuts of Theta nonsingular over F_11 "
          "-> the graph state is AME(18,11)")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
