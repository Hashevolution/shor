"""
Hybrid (C) + Regev b-trick 만 측정 — 큰 N 에서 빠른 검증.

(C) only 와 b-trick only 는 N=437 결과로 충분히 검증됨. 큰 N 에서는 hybrid 만
측정 — 정리 5 의 E[K] ≤ 1.07 (d=4) 예측이 큰 N 에서도 hold 하는지.

실행:
    python -m experiments.hybrid_large_n
"""

from __future__ import annotations
import math
import random
import sys

import numpy as np

from classical import classical_order
from multi_base import (
    convergent_denominators, divisors, minimize_order,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


def lambda_semiprime(N: int) -> int:
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            return math.lcm(p - 1, N // p - 1)
    raise ValueError


def hybrid_factor(
    N: int, setup, runs: list, Q: int,
) -> tuple[int | None, int]:
    """Hybrid 알고리즘 1 회 실행. (factor, K_used) 반환."""
    from multi_base import MultiBaseState, factor_from_exponent
    state = MultiBaseState()
    for K, run in enumerate(runs, start=1):
        for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
            cands = set(convergent_denominators(ki, Q, N - 1))
            if state.L > 1:
                cands.update(divisors(state.L))
            valid = [d for d in cands if d > 0 and pow(ai, d, N) == 1]
            if valid:
                r = minimize_order(ai, N, min(valid))
                if r > 0 and r == classical_order(ai, N):
                    state.update(ai, r)
                    # b-trick
                    b_pow = pow(bi, r, N)
                    if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                        for delta in (-1, 1):
                            g = math.gcd((b_pow + delta) % N, N)
                            if 1 < g < N:
                                return g, K
        # L 이 짝수면 factor_from_exponent
        if state.L > 1:
            rng = random.Random(0)
            res = factor_from_exponent(N, state.L, rng, max_attempts=5)
            if res and 1 < res.factor < N:
                return res.factor, K
    return None, len(runs)


def measure_hybrid(
    N: int, d: int, n_trials: int = 30, max_runs: int = 10,
    noise_kwargs: dict | None = None, corrupt_prob: float = 0.0,
    seed: int = 0,
) -> dict:
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    noise_kwargs = noise_kwargs or {}

    Ks = []
    succ = 0
    for t in range(n_trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)
        setup = regev_setup_bases(N, d, rng_py)
        runs = [
            simulate_regev_run(
                setup.a, N, Q, rng_np,
                corrupt_prob=corrupt_prob, noise_kwargs=noise_kwargs,
            )
            for _ in range(max_runs)
        ]
        factor, K = hybrid_factor(N, setup, runs, Q)
        if factor:
            succ += 1
            Ks.append(K)
        else:
            Ks.append(max_runs)

    return {
        "N": N, "d": d, "trials": n_trials, "noise": noise_kwargs,
        "corrupt": corrupt_prob,
        "mean_K": sum(Ks) / len(Ks),
        "max_K": max(Ks),
        "success": succ,
    }


def main(argv):
    Ns = [int(x) for x in argv[1:]] if len(argv) > 1 else [437, 1147, 2491, 4087]
    print(f"# Hybrid (C)+b-trick — 큰 N 검증")
    print(f"# d=4, n_trials=30, max_runs=10, noise-free\n")
    print(f"  {'N':>5} {'d':>2}  {'mean K':>8}  {'max K':>6}  {'success':>9}")
    for N in Ns:
        r = measure_hybrid(N, d=4, n_trials=30, max_runs=10)
        print(f"  {r['N']:>5} {r['d']:>2}  {r['mean_K']:>8.2f}  "
              f"{r['max_K']:>6}  {r['success']:>3}/{r['trials']}", flush=True)

    print(f"\n# 노이즈 robustness at N=4087")
    for label, kw in [
        ("noise-free", {}),
        ("depol p=0.3", {"depolarizing": 0.3}),
        ("phase σ=1.0", {"phase_sigma": 1.0}),
    ]:
        r = measure_hybrid(4087, d=4, n_trials=20, max_runs=10, noise_kwargs=kw)
        print(f"  {label:<14}: mean K = {r['mean_K']:.2f}, "
              f"success {r['success']}/{r['trials']}", flush=True)


if __name__ == "__main__":
    main(sys.argv)
