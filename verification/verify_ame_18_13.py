"""Stand-alone verifier for AME(18,13) and AME(17,13).

Standard library only, exact integer arithmetic, no imports from the rest
of the repository: this file builds its own F_169, its own linear algebra,
and re-establishes every claim of the paper from the nine embedded field
elements of the group-algebra generator.  It is deliberately independent
of the search engine and of verify_general.py so that agreement between
them is evidence, not shared code.

    python verify_ame_18_13.py     ->  ALL CHECKS PASSED

Conventions.  F_169 = F_13(beta) with beta^2 + beta + 2 = 0; an element
c0 + c1*beta is encoded as the integer c0 + 13*c1.  Conjugation is the
Frobenius x -> x^13.  A code C <= F_169^18 is Hermitian self-dual when
C = C^{perp_h} for <u,v>_h = sum_j u_j conj(v_j).
"""
from itertools import combinations

P = 13
Q = 169
N = 18
K = 9

# ---------------------------------------------------------------- F_169
# beta^2 = -beta - 2 = 12*beta + 11  (mod 13)
def _mul_raw(x, y):
    a0, a1 = x % P, x // P
    b0, b1 = y % P, y // P
    lo = a0 * b0
    mid = a0 * b1 + a1 * b0
    hi = a1 * b1                       # coefficient of beta^2
    c0 = (lo + 11 * hi) % P
    c1 = (mid + 12 * hi) % P
    return c0 + P * c1


def _add_raw(x, y):
    return ((x % P + y % P) % P) + P * ((x // P + y // P) % P)


MUL = [[_mul_raw(x, y) for y in range(Q)] for x in range(Q)]
ADD = [[_add_raw(x, y) for y in range(Q)] for x in range(Q)]
NEG = [((-(x % P)) % P) + P * ((-(x // P)) % P) for x in range(Q)]


def mul(x, y):
    return MUL[x][y]


def add(x, y):
    return ADD[x][y]


def sub(x, y):
    return ADD[x][NEG[y]]


def power(x, e):
    r, b = 1, x
    while e:
        if e & 1:
            r = MUL[r][b]
        b = MUL[b][b]
        e >>= 1
    return r


CONJ = [power(x, P) for x in range(Q)]
INV = [0] * Q
for _x in range(1, Q):
    INV[_x] = power(_x, Q - 2)

# field self-checks
assert all(MUL[x][INV[x]] == 1 for x in range(1, Q))
assert all(CONJ[CONJ[x]] == x for x in range(Q))              # Frobenius^2 = id
assert sum(1 for x in range(Q) if CONJ[x] == x) == P           # fixed field
# closed form of the conjugate: (c0 + c1 b) -> (c0 - c1) - c1 b
for _x in range(Q):
    _c0, _c1 = _x % P, _x // P
    assert CONJ[_x] == (_c0 - _c1) % P + P * ((-_c1) % P)

# ------------------------------------------------- the group-algebra data
# a : Z_3 x Z_3 -> F_169, indexed a[u][v], coordinate index u + 3v
A_GEN = [[(6, 12), (9, 9), (1, 7)],
         [(1, 10), (11, 7), (0, 5)],
         [(8, 12), (6, 9), (2, 7)]]
a = {}
for u in range(3):
    for v in range(3):
        c0, c1 = A_GEN[u][v]
        a[(u, v)] = c0 + P * c1


def coord(t):
    return (t % 3, t // 3)


A = [[a[((coord(y)[0] - coord(x)[0]) % 3,
         (coord(y)[1] - coord(x)[1]) % 3)] for y in range(K)]
     for x in range(K)]

# ------------------------------------------------------- linear algebra
def det169(M):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c]:
                piv = r
                break
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = NEG[d]
        d = MUL[d][M[c][c]]
        iv = INV[M[c][c]]
        for r in range(c + 1, n):
            if M[r][c]:
                f = MUL[M[r][c]][iv]
                Mr, Mc = M[r], M[c]
                for cc in range(c, n):
                    Mr[cc] = sub(Mr[cc], MUL[f][Mc[cc]])
    return d


def rank_p(M):
    """Rank over F_13 (plain integers mod 13)."""
    M = [[x % P for x in row] for row in M]
    if not M:
        return 0
    rows, cols = len(M), len(M[0])
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if M[rr][c]:
                piv = rr
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        iv = pow(M[r][c], P - 2, P)
        M[r] = [(iv * x) % P for x in M[r]]
        for rr in range(rows):
            if rr != r and M[rr][c]:
                f = M[rr][c]
                M[rr] = [(M[rr][j] - f * M[r][j]) % P for j in range(cols)]
        r += 1
    return r


# ------------------------------------------------------------- checks
def check_group_algebra():
    """The twisted convolution identity, equivalent to A Abar^T = -I."""
    bad = 0
    for g in [(u, v) for u in range(3) for v in range(3)]:
        s = 0
        for h in [(u, v) for u in range(3) for v in range(3)]:
            hg = ((h[0] + g[0]) % 3, (h[1] + g[1]) % 3)
            s = add(s, MUL[a[hg]][CONJ[a[h]]])
        want = NEG[1] if g == (0, 0) else 0
        if s != want:
            bad += 1
    assert bad == 0
    print("check 1  twisted convolution  sum_h a(h+g) conj(a(h)) "
          "= -delta_{g,0}   (9 sums)  OK")


def check_unitarity():
    cnt = 0
    for i in range(K):
        for j in range(K):
            s = 0
            for t in range(K):
                s = add(s, MUL[A[i][t]][CONJ[A[j][t]]])
            assert s == (NEG[1] if i == j else 0)
            cnt += 1
    print(f"check 2  A Abar^T = -I_9 over F_169                    "
          f"({cnt} products)  OK")


def check_mds():
    total = 0
    for size in range(1, K + 1):
        for rs in combinations(range(K), size):
            for cs in combinations(range(K), size):
                if det169([[A[r][c] for c in cs] for r in rs]) == 0:
                    raise AssertionError(f"zero minor {rs} {cs}")
                total += 1
    assert total == 48619
    print(f"check 3  every square submatrix of A nonsingular       "
          f"({total} minors)  OK")
    return total


def expand():
    """F_13-basis of C in symplectic coordinates over the basis (1, beta)."""
    G = [[1 if i == j else 0 for j in range(K)] + A[i] for i in range(K)]
    beta = P                                   # 0 + 1*beta
    rows = []
    for row in G:
        rows.append(row)
        rows.append([MUL[beta][x] for x in row])
    cnt = 0
    for u in rows:
        for v in rows:
            s = 0
            for j in range(N):
                s = add(s, MUL[u[j]][CONJ[v[j]]])
            assert s == 0
            cnt += 1
    print(f"check 4  expanded F_13-basis Hermitian self-orthogonal "
          f"({cnt} products)  OK")
    M = []
    for u in rows:
        M.append([x % P for x in u] + [x // P for x in u])
    cnt = 0
    for i in range(N):
        for j in range(N):
            s = 0
            for t in range(N):
                s += M[i][t] * M[j][N + t] - M[j][t] * M[i][N + t]
            assert s % P == 0
            cnt += 1
    print(f"check 5  symplectic self-orthogonality over F_13       "
          f"({cnt} products)  OK")
    return M


def to_graph(M):
    M = [row[:] for row in M]
    for _ in range(4 * N + 16):
        if rank_p([row[:N] for row in M]) == N:
            break
        base = rank_p([row[:N] for row in M])
        moved = False
        for j in range(N):
            M2 = [row[:] for row in M]
            for row in M2:
                row[j], row[N + j] = (-row[N + j]) % P, row[j]
            if rank_p([row[:N] for row in M2]) > base:
                M, moved = M2, True
                break
        if not moved:
            for row in M:
                row[0], row[N] = (-row[N]) % P, row[0]
    assert rank_p([row[:N] for row in M]) == N
    r = 0
    for c in range(N):
        piv = next(rr for rr in range(r, N) if M[rr][c] % P)
        M[r], M[piv] = M[piv], M[r]
        iv = pow(M[r][c] % P, P - 2, P)
        M[r] = [(iv * x) % P for x in M[r]]
        for rr in range(N):
            if rr != r and M[rr][c] % P:
                f = M[rr][c]
                M[rr] = [(M[rr][j] - f * M[r][j]) % P
                         for j in range(2 * N)]
        r += 1
    for i in range(N):
        for j in range(N):
            assert M[i][j] % P == (1 if i == j else 0)
    T = [[M[i][N + j] % P for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            assert T[i][j] == T[j][i]
    print("check 6  graph form [I | Theta] reached, Theta symmetric  OK")
    return T


def check_cuts(T):
    bad = tot = 0
    for S in combinations(range(N), K):
        if 0 not in S:
            continue
        Sc = [x for x in range(N) if x not in S]
        tot += 1
        if rank_p([[T[i][j] for j in Sc] for i in S]) < K:
            bad += 1
    assert bad == 0 and tot == 24310
    print(f"check 7  every balanced cut of Theta nonsingular       "
          f"({tot} blocks)  OK  ==> AME(18,13)")
    return tot


def check_deletions(T):
    npass = tot = 0
    for v in range(N):
        keep = [i for i in range(N) if i != v]
        T17 = [[T[i][j] for j in keep] for i in keep]
        bad = 0
        for S in combinations(range(17), 8):
            Sc = [x for x in range(17) if x not in S]
            tot += 1
            if rank_p([[T17[i][j] for j in Sc] for i in S]) < 8:
                bad += 1
        npass += (bad == 0)
    assert npass == N
    print(f"check 8  all {N} single-vertex deletions AME           "
          f"({tot} rank tests)  OK  ==> AME(17,13)")


if __name__ == "__main__":
    print("Stand-alone verification of AME(18,13) and AME(17,13)")
    print("F_169 = F_13(beta), beta^2 + beta + 2 = 0; "
          "element c0 + c1*beta encoded as c0 + 13*c1\n")
    check_group_algebra()
    check_unitarity()
    check_mds()
    M = expand()
    T = to_graph(M)
    check_cuts(T)
    check_deletions(T)
    print("\nTheta =")
    for row in T:
        print("  " + " ".join(f"{x:2d}" for x in row))
    print("\nALL CHECKS PASSED")
