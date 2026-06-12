"""
SR 증폭 메커니즘 탐색.

가설:
- M1: SR 가 per-coordinate 효과면, larger d 가 효과 증폭.
- M2: 노이즈 의 k-broadening 을 인공적으로 흉내내는 near-k enumeration.

M1 실험: d=2, 4, 8, 16 × phase σ scan.
M2 실험: noise-free 에서 (k, k±1, k±2) 후보 vs 표준 (k 만).

실행:
    python -m experiments.sr_amplify
"""

from __future__ import annotations
import math
import random
import sys

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    factor_from_exponent,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


def measure_hybrid_extended(
    N: int, d: int, noise_kwargs: dict, trials: int = 200, max_runs: int = 20,
    seed: int = 0, k_neighbors: int = 0,
) -> float:
    """Hybrid + 선택적 k-neighbor enumeration.

    k_neighbors > 0 면: 각 측정 k 에 대해 (k-k_n, ..., k, ..., k+k_n) 모두 시도.
    """
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    Ks = []
    for t in range(trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)
        setup = regev_setup_bases(N, d, rng_py)
        state = MultiBaseState()
        K_found = max_runs
        for K in range(1, max_runs + 1):
            run = simulate_regev_run(
                setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs,
            )
            for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
                # k-neighbor enumeration
                ks_to_try = [ki]
                if k_neighbors > 0:
                    for delta in range(1, k_neighbors + 1):
                        ks_to_try.extend([(ki + delta) % Q, (ki - delta) % Q])

                for k_try in ks_to_try:
                    cands = set(convergent_denominators(k_try, Q, N - 1))
                    if state.L > 1:
                        cands.update(divisors(state.L))
                    valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
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
                                        K_found = K
                                        break
                                if K_found < max_runs:
                                    break
                            break
                if K_found < max_runs:
                    break
            if K_found < max_runs:
                break
            if state.L > 1:
                rng_f = random.Random(t)
                res = factor_from_exponent(N, state.L, rng_f, max_attempts=5)
                if res and 1 < res.factor < N:
                    K_found = K
                    break
        Ks.append(K_found)
    return sum(Ks) / len(Ks)


def m1_larger_d(N: int = 437, trials: int = 200):
    """M1: larger d 가 SR 증폭하는가."""
    print(f"\n# M1: d 변화에 따른 SR 효과 (N={N}, {trials} trials)")
    print(f"  {'d':>3}  {'σ=0':>8}  {'σ=0.05':>8}  {'σ=0.10':>8}  "
          f"{'σ=0.20':>8}  {'baseline-σ_min':>15}  {'SR %':>8}")
    for d in [2, 4, 8]:
        results = {}
        for sigma in [0.0, 0.05, 0.10, 0.20]:
            K = measure_hybrid_extended(
                N, d, {"phase_sigma": sigma} if sigma > 0 else {},
                trials=trials,
            )
            results[sigma] = K
        baseline = results[0.0]
        min_K = min(results.values())
        improvement = (baseline - min_K) / baseline * 100 if baseline > 0 else 0
        print(f"  {d:>3}  {results[0.0]:>8.3f}  {results[0.05]:>8.3f}  "
              f"{results[0.10]:>8.3f}  {results[0.20]:>8.3f}  "
              f"{baseline - min_K:>15.3f}  {improvement:>7.2f}%", flush=True)


def m2_near_k(N: int = 437, d: int = 4, trials: int = 200):
    """M2: 노이즈 없이 k-neighbor enumeration 으로 SR 흉내."""
    print(f"\n# M2: k-neighbor enumeration (N={N}, d={d}, {trials} trials, noise-free)")
    print(f"  {'k_neighbors':>12}  {'hybrid K':>10}  {'vs baseline':>12}")
    baseline = None
    for k_n in [0, 1, 2, 3, 5, 10]:
        K = measure_hybrid_extended(
            N, d, {}, trials=trials, k_neighbors=k_n,
        )
        if baseline is None:
            baseline = K
        diff = K - baseline
        print(f"  {k_n:>12}  {K:>10.3f}  {diff:>+12.3f}", flush=True)


def main(argv):
    if "--m1" in argv:
        m1_larger_d()
    elif "--m2" in argv:
        m2_near_k()
    else:
        m1_larger_d()
        m2_near_k()


if __name__ == "__main__":
    main(sys.argv)
