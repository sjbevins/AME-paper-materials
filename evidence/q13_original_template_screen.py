"""Complete classification of untwisted abelian permutation templates at a
given cell (n, q), by character combinatorics alone (no linear algebra).

For an abelian group H acting on n coordinates with orbits O_1,...,O_r of
stabilisers K_1,...,K_r, the permutation module decomposes with

    mult(chi) = #{ i : chi|_{K_i} = 1 },

the Hermitian form pairs V_chi with V_{chi^{-q}}, and chi is self-paired
iff chi^{q+1} = 1.  A template admits an invariant Hermitian self-dual
code only if every self-paired character has EVEN multiplicity (a
nondegenerate Hermitian space of odd dimension has no Lagrangian).

Family size of a surviving template:

    prod over self-paired chi   L(mult(chi))
  * prod over paired {chi,chi'} S(mult(chi))

    L(2e) = prod_{i=1..e} (q^{2i-1} + 1)      (Lagrangians, dim 2e)
    S(d)  = sum_{t=0..d} gauss_binom(q^2,d,t) (U of any dim, plus its
                                               annihilator in the partner)

Characters must take values in F_{q^2}, i.e. exp(H) | q^2-1, and |H| must
be coprime to char F_{q^2} or the module is not semisimple.

Control: the row for Z3xZ3 with orbits (9,9) at (18,13) must reproduce the
family size 12,252,982,784 that the search engine computed independently.
"""
from __future__ import annotations

import itertools
from fractions import Fraction


def abelian_groups(m):
    """All abelian groups of order m, as tuples of invariant factors."""
    def parts(k):
        if k == 0:
            yield ()
            return
        for first in range(k, 0, -1):
            for rest in parts(k - first):
                if not rest or first >= rest[0]:
                    yield (first,) + rest
    facs = {}
    x = m
    d = 2
    while d * d <= x:
        while x % d == 0:
            facs[d] = facs.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        facs[x] = facs.get(x, 0) + 1
    per_prime = []
    for p, e in sorted(facs.items()):
        opts = []
        for lam in parts(e):
            opts.append(tuple(p ** a for a in sorted(lam)))
        per_prime.append(opts)
    out = []
    for combo in itertools.product(*per_prime) if per_prime else [()]:
        cyc = []
        maxlen = max((len(c) for c in combo), default=0)
        for i in range(maxlen):
            f = 1
            for c in combo:
                cc = list(c) + [1] * (maxlen - len(c))
                f *= cc[i]
            cyc.append(f)
        out.append(tuple(sorted(cyc)))
    return sorted(set(out)) if out else [()]


def elements(inv):
    return list(itertools.product(*[range(d) for d in inv]))


def subgroups(inv):
    """All subgroups of the abelian group with invariant factors inv."""
    els = elements(inv)
    n = len(inv)

    def add(a, b):
        return tuple((a[i] + b[i]) % inv[i] for i in range(n))

    def span(gens):
        S = {tuple([0] * n)}
        frontier = [tuple([0] * n)]
        while frontier:
            x = frontier.pop()
            for g in gens:
                y = add(x, g)
                if y not in S:
                    S.add(y)
                    frontier.append(y)
        return frozenset(S)

    subs = set()
    for r in range(0, 4):
        for gens in itertools.combinations(els, r):
            subs.add(span(gens))
    return sorted(subs, key=lambda S: (len(S), sorted(S)))


def char_trivial_on(inv, c, K):
    """Is chi_c trivial on subgroup K?  chi_c(h) = prod zeta_{d_i}^{c_i h_i}."""
    for h in K:
        s = Fraction(0)
        for i, d in enumerate(inv):
            s += Fraction(c[i] * h[i], d)
        if s.denominator != 1:
            return False
    return True


def gauss_binom(Q, d, t):
    num = den = 1
    for i in range(t):
        num *= Q ** (d - i) - 1
        den *= Q ** (i + 1) - 1
    return num // den


def lagrangians(q, d):
    assert d % 2 == 0
    out = 1
    for i in range(1, d // 2 + 1):
        out *= q ** (2 * i - 1) + 1
    return out


def paired_options(Q, d):
    return sum(gauss_binom(Q, d, t) for t in range(d + 1))


def orbit_multisets(sizes, n):
    """All multisets of orbit sizes (from the available list) summing to n."""
    sizes = sorted(set(sizes), reverse=True)

    def rec(i, rem, cur):
        if rem == 0:
            yield tuple(cur)
            return
        if i == len(sizes):
            return
        s = sizes[i]
        k = 0
        while k * s <= rem:
            yield from rec(i + 1, rem - k * s, cur + [s] * k)
            k += 1
    return list(rec(0, n, []))


def screen(n, q, verbose=True, max_report=None):
    Q = q * q
    p = None
    x = Q
    d = 2
    while d * d <= x:
        if x % d == 0:
            p = d
            break
        d += 1
    if p is None:
        p = Q if all(Q % t for t in range(2, int(Q ** 0.5) + 1)) else None
    # char of F_{q^2} = char of F_q
    ch = 2
    while q % ch:
        ch += 1
    rows = []
    for m in range(1, n + 1):
        if m % ch == 0:
            continue                                  # not semisimple
        for inv in abelian_groups(m):
            if not inv:
                continue
            if (Q - 1) % max(inv) != 0:
                continue                              # chars not in F_{q^2}
            subs = subgroups(inv)
            by_size = {}
            for K in subs:
                by_size.setdefault(m // len(K), []).append(K)
            trivial_cache = {}
            for K in subs:
                for c in elements(inv):
                    trivial_cache[(K, c)] = char_trivial_on(inv, c, K)
            for shape in orbit_multisets(list(by_size), n):
                # multiset of stabilisers per orbit size (order is irrelevant,
                # so enumerate combinations_with_replacement, not the product)
                counts = {}
                for s in shape:
                    counts[s] = counts.get(s, 0) + 1
                per_size = [list(itertools.combinations_with_replacement(
                    by_size[s], counts[s])) for s in sorted(counts)]
                for groups in itertools.product(*per_size):
                    pick = [K for grp in groups for K in grp]
                    mult = {}
                    for c in elements(inv):
                        mult[c] = sum(1 for K in pick
                                      if trivial_cache[(K, c)])
                    ok = True
                    selfp, pairs = [], []
                    done = set()
                    for c in elements(inv):
                        if c in done:
                            continue
                        cc = tuple((-q * c[i]) % inv[i]
                                   for i in range(len(inv)))
                        if mult[c] != mult[cc]:
                            ok = False
                            break
                        if cc == c:
                            if mult[c] % 2:
                                ok = False
                                break
                            if mult[c]:
                                selfp.append(mult[c])
                            done.add(c)
                        else:
                            if mult[c]:
                                pairs.append(mult[c])
                            done.add(c)
                            done.add(cc)
                    if not ok:
                        continue
                    fam = 1
                    for dd in selfp:
                        fam *= lagrangians(q, dd)
                    for dd in pairs:
                        fam *= paired_options(Q, dd)
                    rows.append((fam, inv, shape, tuple(sorted(selfp)),
                                 tuple(sorted(pairs))))
    rows.sort()
    if verbose:
        print(f"parity-viable untwisted abelian templates at "
              f"(n,q) = ({n},{q}):", flush=True)
        print(f"{'family':>26}  {'H':<14} {'orbits':<18} "
              f"self-paired dims / paired dims", flush=True)
        for fam, inv, shape, sp, pr in (rows[:max_report] if max_report
                                        else rows):
            hname = "x".join(f"Z{d}" for d in inv)
            print(f"{fam:>26,}  {hname:<14} {str(shape):<18} "
                  f"{list(sp)} / {list(pr)}", flush=True)
        print(f"total viable: {len(rows)}", flush=True)
    return rows


def profiles(rows):
    """Collapse the raw template list to distinct block profiles."""
    agg = {}
    for fam, inv, shape, sp, pr in rows:
        key = (fam, inv, shape, sp, pr)
        agg[key] = agg.get(key, 0) + 1
    out = sorted(agg.items())
    return out


if __name__ == "__main__":
    rows = screen(18, 13, verbose=False)
    hit = [r for r in rows if r[1] == (3, 3) and r[2] == (9, 9)]
    assert hit, "the Z3xZ3 (9,9) template is missing from the screen"
    assert hit[0][0] == 12_252_982_784, f"family {hit[0][0]} != engine value"
    prof = profiles(rows)
    print(f"parity-viable untwisted abelian permutation templates at "
          f"(n,q) = (18,13): {len(rows)} templates in {len(prof)} distinct "
          f"block profiles\n", flush=True)
    print(f"{'family':>26}  {'H':<10} {'orbits':<16} {'self-p':<10} "
          f"{'paired':<16} {'#':>5}", flush=True)
    for (fam, inv, shape, sp, pr), cnt in prof[:14]:
        hname = "x".join(f"Z{d}" for d in inv)
        print(f"{fam:>26,}  {hname:<10} {str(shape):<16} {str(list(sp)):<10} "
              f"{str(list(pr)):<16} {cnt:>5}", flush=True)
    print(" ...", flush=True)
    for (fam, inv, shape, sp, pr), cnt in prof:
        if inv == (3, 3) and shape == (9, 9):
            hname = "x".join(f"Z{d}" for d in inv)
            rank = [i for i, x in enumerate(prof) if x[0][0] == fam][0]
            print(f"{fam:>26,}  {hname:<10} {str(shape):<16} "
                  f"{str(list(sp)):<10} {str(list(pr)):<16} {cnt:>5}"
                  f"   <-- searched (rank {rank + 1} by family size)",
                  flush=True)
    print("\nCONTROL: Z3xZ3 (9,9) family = 12,252,982,784 "
          "matches the search engine", flush=True)
