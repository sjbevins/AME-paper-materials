"""Derive and certify AME(17,11) from the verified AME(18,11) state.

Projection lemma: measuring one qudit of an AME(2m, q) state in the
computational basis leaves an AME(2m-1, q) state on the rest -- for any
(m-1)-subset T of the survivors, T plus the measured party had m parties
and was maximally mixed in the original state, so conditioned on the
(uniformly random) outcome the marginal on T is maximally mixed.  For a
graph state, projecting vertex v onto |0> deletes vertex v: the state on
the remaining vertices is the graph state of Theta with row/column v
removed.

Certificate for odd n = 17: |Theta17> is AME(17,11) iff for EVERY
8-subset S of the 17 vertices the 8 x 9 block Theta17[S, S^c] has full
row rank 8 over F_11 (C(17,8) = 24,310 checks; the complementary 9x8
blocks are transposes, and smaller cuts follow by monotonicity).

Two independent rank implementations (own integer elimination vs numpy
pivoted elimination) must agree on every subset.  As a bonus census the
numpy path then checks all 18 possible vertex deletions.
"""
from __future__ import annotations
import itertools
import sys
from pathlib import Path

import numpy as np

P = 11
N18, N17, K = 18, 17, 8

THETA18 = np.load(Path(__file__).with_name('theta_ame1811.npy')).astype(np.int64) % P

INVP = [0] + [pow(x, P - 2, P) for x in range(1, P)]


def rank_own(M):
    """Pure-python integer Gaussian elimination rank over F_11."""
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0])
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] % P), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        iv = INVP[M[r][c] % P]
        M[r] = [(x * iv) % P for x in M[r]]
        for i in range(r + 1, rows):
            f = M[i][c] % P
            if f:
                M[i] = [(M[i][t] - f * M[r][t]) % P for t in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def rank_np(M):
    M = M.copy() % P
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c] % P), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * INVP[M[r, c] % P]) % P
        for i in range(r + 1, rows):
            if M[i, c] % P:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
        if r == rows:
            break
    return r


def certify(theta17, dual=True, tag="v"):
    n = theta17.shape[0]
    assert n == N17
    bad = 0
    tot = 0
    for S in itertools.combinations(range(N17), K):
        Sc = [x for x in range(N17) if x not in S]
        B = theta17[np.ix_(list(S), Sc)]
        r1 = rank_np(B)
        if dual:
            r2 = rank_own([[int(B[i, j]) for j in range(B.shape[1])]
                           for i in range(B.shape[0])])
            assert r1 == r2, f"[{tag}] implementations disagree at {S}"
        if r1 < K:
            bad += 1
        tot += 1
    return bad, tot


if __name__ == "__main__":
    # primary certificate: delete vertex 17 (the last), dual-implementation
    T17 = np.delete(np.delete(THETA18, 17, axis=0), 17, axis=1)
    bad, tot = certify(T17, dual=True, tag="delete-17")
    print(f"primary (delete vertex 18): deficient 8-cuts = {bad}/{tot} "
          f"(dual-implementation agreement on all subsets)")
    assert bad == 0 and tot == 24310
    np.save(Path(__file__).with_name('theta_ame1711_regenerated.npy'), T17)
    print("AME(17,11) CERTIFIED: all 24,310 cut blocks have full rank 8 "
          "over F_11")
    print("Theta17 ="); print(T17.tolist())
    # bonus census: every single-vertex deletion also yields AME(17,11)
    ok = 0
    for v in range(N18):
        Tv = np.delete(np.delete(THETA18, v, axis=0), v, axis=1)
        b, _ = certify(Tv, dual=False, tag=f"delete-{v}")
        ok += (b == 0)
    print(f"bonus census: {ok}/18 single-vertex deletions give AME(17,11) "
          f"(projection lemma predicts 18/18)")
    assert ok == 18
    print("ALL CHECKS PASSED")
