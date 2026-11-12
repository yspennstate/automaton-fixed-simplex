"""Exact checks for the automaton A_k: the fixed points (Theorem 1.1), the mean offspring matrix M of the pivot
tree of the merged automaton B, and the identity of Remark 2.1, for offspring laws with finite support.

Run: python check_automaton.py
"""
from fractions import Fraction as F
from itertools import product


def A(children, k):
    """Output of A_k; state 0 is a, states 1..k are b_1..b_k."""
    odd = [c for c in children if c != 0]
    return odd[0] if len(odd) == 1 else 0


def psi(chi, x, k):
    """Automaton distributional map by enumerating every ordered child configuration."""
    out = [F(0)] * (k + 1)
    for n, w in chi.items():
        for conf in product(range(k + 1), repeat=n):
            pr = w
            for c in conf:
                pr *= x[c]
            out[A(conf, k)] += pr
    return out


def G1(chi, s):
    return sum(w * n * s ** (n - 1) for n, w in chi.items() if n >= 1)


def G2(chi, s):
    return sum(w * n * (n - 1) * s ** (n - 2) for n, w in chi.items() if n >= 2)


# offspring laws with a rational s_*: G'(s_*) = 1
LAWS = {
    "binary": ({2: F(1)}, F(1, 2)),
    "unary-binary r=2/5": ({1: F(2, 5), 2: F(3, 5)}, F(1, 2)),
    "leaf (1/4)d0+(3/4)d2": ({0: F(1, 4), 2: F(3, 4)}, F(2, 3)),
    "leaf (1/4)d0+(3/4)d3": ({0: F(1, 4), 3: F(3, 4)}, F(2, 3)),
    "leaf (5/32)d0+(27/32)d4": ({0: F(5, 32), 4: F(27, 32)}, F(2, 3)),
}

n_ok = 0
for name, (chi, s_star) in LAWS.items():
    assert sum(chi.values()) == 1
    assert G1(chi, s_star) == 1, name
    for k in (2, 3):
        # every law with x(a) = s_* is fixed
        for mu in ([F(1)] + [F(0)] * (k - 1), [F(1, k)] * k, [F(j + 1, k * (k + 1) // 2) for j in range(k)]):
            x = [s_star] + [(1 - s_star) * m for m in mu]
            assert psi(chi, x, k) == x, (name, k, mu)
            n_ok += 1
        # delta_a is fixed; laws with x(a) != s_* are not
        assert psi(chi, [F(1)] + [F(0)] * k, k) == [F(1)] + [F(0)] * k
        for s in (F(1, 7), F(1, 3), F(5, 6), F(0)):
            if s == s_star:
                continue
            x = [s] + [(1 - s) / k] * k
            assert psi(chi, x, k) != x, (name, k, s)
            n_ok += 1
        # the closed form of the map
        x = [F(2, 5)] + [F(3, 5) / k] * k
        y = psi(chi, x, k)
        assert y[1] == x[1] * G1(chi, x[0]) and y[0] == 1 - (1 - x[0]) * G1(chi, x[0])
        n_ok += 1

    # pivot tree of B at (s_*, p): count pivotal children by flipping, exactly
    s, p = s_star, 1 - s_star
    B = lambda conf: 1 if sum(conf) == 1 else 0  # 1 = b, 0 = a
    M = {(u, v): F(0) for u in (1, 0) for v in (1, 0)}
    mass = {1: F(0), 0: F(0)}
    for n, w in chi.items():
        for conf in product((0, 1), repeat=n):
            pr = w
            for c in conf:
                pr *= p if c else s
            col = B(conf)
            mass[col] += pr
            for i, c in enumerate(conf):
                flipped = conf[:i] + (1 - c,) + conf[i + 1:]
                if B(flipped) != col:
                    M[(col, c)] += pr
    assert mass[1] == p and mass[0] == s
    M = {key: val / mass[key[0]] for key, val in M.items()}
    g2 = G2(chi, s)
    want = {(1, 1): F(1), (1, 0): s * g2, (0, 1): p * p * g2 / s, (0, 0): F(1)}
    assert M == want, (name, M, want)
    # spectral radius 1 + p G''(s): (t-1)^2 = M_ba M_ab = (p G''(s))^2
    assert M[(1, 0)] * M[(0, 1)] == (p * g2) ** 2
    print(f"{name}: s_* = {s_star}, M = [[1, {M[(1, 0)]}], [{M[(0, 1)]}, 1]], spectral radius 1 + {p * g2} > 1")
    n_ok += 1

# Remark 2.1: chi(n) = 1/(n(n-1)), automaton '1 iff at most one child is 0'; the sum to degree K is q - q^K/K
for q in (F(0), F(1, 4), F(1, 2), F(9, 10), F(1)):
    K = 60
    tot = sum((q ** n + n * (1 - q) * q ** (n - 1)) / (n * (n - 1)) for n in range(2, K + 1))
    assert tot == q - q ** K / K, q
    n_ok += 1
print(f"all {n_ok} checks passed")
