"""
magic_proofs_check.py — magic-results.md 의 명제들을 작은 n에서 assert (회귀 검증).

검증:
  보조정리 1  : 평탄상태 M₂=0 ⟺ 받침이 아핀부분공간.
  따름정리 1  : 그래프상태 M₂=0 ⟺ f 아핀;  비선형 f → >0;  Shor modexp → >0.
  명제 2      : Grover(M=1) 닫힌형 = sre2.
  명제 3      : Grover 정점 magic → 3 (큰 n).

Reproduction:
  python -u -m experiments.magic_proofs_check
"""
from __future__ import annotations

import math
import time

import numpy as np

from magic import sre2, sre2_bruteforce


def flat(support, nq):
    psi = np.zeros(2 ** nq)
    psi[list(support)] = 1.0
    return psi / math.sqrt(len(support))


def subspace(basis, nq, offset=0):
    pts = {offset}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return sorted(pts)


def graph_state(values, n_in, n_out):
    psi = np.zeros(2 ** (n_in + n_out))
    for x in range(2 ** n_in):
        psi[(x << n_out) | values[x]] = 1.0
    return psi / math.sqrt(2 ** n_in)


def sre2_grover_closed(n, a, b):
    N = 2 ** n
    S = 1.0
    S += (N - 1) * (a * a - b * b) ** 4
    S += (N - 1) * (b * b * (N - 2) + 2 * a * b) ** 4
    S += (N - 1) * (N / 2 - 1) * (2 * b * (a - b)) ** 4
    return -math.log2(S / N)


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(7)
    checks = 0

    # 보조정리 1: 아핀부분공간 → 0
    nq = 5
    for basis, off in [([1, 2, 4], 0), ([3, 12], 5), ([1, 2, 4, 8, 16], 0), ([7], 0)]:
        assert abs(sre2(flat(subspace(basis, nq, off), nq))) < 1e-9
        checks += 1
    # 비아핀 받침 → >0  (닫혀있지 않은 집합)
    for S in [[0, 1, 2, 4], [0, 1, 2, 3, 5], [1, 2, 4, 8, 7, 11]]:
        assert sre2(flat(S, nq)) > 1e-6
        checks += 1

    # 보조정리 1 ↔ brute force (n≤4)
    for _ in range(5):
        k = rng.integers(2, 8)
        S = sorted(rng.choice(16, size=k, replace=False).tolist())
        assert abs(sre2(flat(S, 4)) - sre2_bruteforce(flat(S, 4))) < 1e-9
        checks += 1

    # 따름정리 1: Simon 선형오라클 → 0
    for n in [2, 3, 4]:
        for s in [1, 3]:
            p = (s & -s).bit_length() - 1
            vals = [x ^ (s if (x >> p) & 1 else 0) for x in range(2 ** n)]
            assert abs(sre2(graph_state(vals, n, n))) < 1e-9
            checks += 1
    # 비선형 f → >0
    base = [x ^ (3 if (x >> 0) & 1 else 0) for x in range(16)]
    nl = list(base)
    for x in range(16):
        if (x >> 0) & 1 and (x >> 1) & 1:
            nl[x] ^= (1 << 2)
    assert sre2(graph_state(nl, 4, 4)) > 1e-6
    checks += 1
    # Shor modexp → >0
    vals = [pow(7, x, 15) for x in range(2 ** 5)]
    assert sre2(graph_state(vals, 5, 4)) > 1e-6
    checks += 1

    # 명제 2: Grover 닫힌형 = sre2
    for n in [3, 4, 6, 8]:
        N = 2 ** n
        for _ in range(3):
            ang = rng.uniform(0.1, math.pi / 2 - 0.1)
            a, b = math.sin(ang), math.cos(ang) / math.sqrt(N - 1)
            psi = np.full(N, b)
            psi[0] = a
            assert abs(sre2(psi) - sre2_grover_closed(n, a, b)) < 1e-9
            checks += 1

    # 명제 3: 정점 → 3
    angs = np.linspace(0.30, 0.60, 1500) * (math.pi / 2)
    peaks = {}
    for n in [10, 20, 30]:
        N = 2 ** n
        a = np.sin(angs)
        b = np.cos(angs) / math.sqrt(N - 1)
        peaks[n] = max(sre2_grover_closed(n, ai, bi) for ai, bi in zip(a, b))
    assert peaks[30] > peaks[20] > peaks[10]
    assert peaks[30] < 3.0 and abs(peaks[30] - 3.0) < 0.01
    checks += 3

    # 명제 2′: 아핀 W → Grover 정점 유한(≤3+ε), flat_W 안정자; 비아핀 W → flat_W magic>0
    def two_amp(n, W, ang):
        N = 2 ** n
        a, b = math.sin(ang) / math.sqrt(len(W)), math.cos(ang) / math.sqrt(N - len(W))
        psi = np.full(N, b)
        psi[list(W)] = a
        return psi

    def _sub(basis):
        pts = {0}
        for bb in basis:
            pts |= {p ^ bb for p in pts}
        return sorted(pts)

    nn = 8
    angs = np.linspace(0.05, math.pi / 2 - 0.05, 120)
    for basis in [[1], [1, 2], [1, 2, 4]]:
        W = _sub(basis)
        assert abs(sre2(flat(W, nn))) < 1e-9                       # 아핀 → flat 안정자
        assert max(sre2(two_amp(nn, W, a)) for a in angs) < 3.0 + 1e-6   # 정점 유한 ≤3
        checks += 2
    Wnl = sorted(rng.choice(2 ** nn, size=8, replace=False).tolist())
    assert sre2(flat(Wnl, nn)) > 1e-6                              # 비아핀 → magic>0
    checks += 1

    elapsed = time.time() - t0
    print(f"ALL PROOFS CHECKED — {checks} assertions passed.")
    print(f"  보조정리 1, 따름정리 1, 명제 2: OK")
    print(f"  명제 3 정점:  n=10 → {peaks[10]:.4f},  n=20 → {peaks[20]:.4f},  "
          f"n=30 → {peaks[30]:.4f}  (→ 3.000)")
    print(f"# Elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
