"""
Stochastic resonance (SR) search — noise-as-feature 가설 검증.

가설: 적정 노이즈 수준이 인수분해를 *돕는다*.
- 노이즈 0: standard 분포, multi-base 누적 ~ K_λ^ideal
- 노이즈 ↑: 측정 분포 변형 → 일부 trial 이 *우연히 더 좋은* convergent 줄 수 있음
- 너무 ↑: 정보 손실 → 회수율 하락

U-shape (= minimum at σ > 0) 이면 quantum SR 인수분해 첫 evidence.

실험 두 가지:
- F1: ideal K_λ (multi-base lcm to λ(N)) vs noise level
- F2: hybrid K (runs to factor) vs noise level

실행:
    python -m experiments.sr_search
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
from noise import simulate_period_finding_noisy
from shor import simulate_period_finding
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


def lambda_of(N):
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            return math.lcm(p - 1, N // p - 1)


def measure_klambda_alg_noisy(
    N: int, noise_kwargs: dict, trials: int = 100, max_bases: int = 50, seed: int = 0,
) -> float:
    """노이즈 하의 (C) 알고리즘이 L=λ(N) 도달까지 base 수. (F1 검증)"""
    lam = lambda_of(N)
    Ks = []
    for t in range(trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)
        state = MultiBaseState()
        K = max_bases
        for k in range(1, max_bases + 1):
            for _retry in range(50):
                a = rng_py.randrange(2, N)
                if math.gcd(a, N) == 1:
                    break
            else:
                continue
            true_r = classical_order(a, N)

            # Fast path
            if state.L > 1 and pow(a, state.L, N) == 1:
                state.update(a, true_r)
                if state.L == lam:
                    K = k
                    break
                continue

            # Slow path with noise
            m = simulate_period_finding_noisy(a, N, rng=rng_np, **noise_kwargs)
            cands = set(convergent_denominators(m.k, m.Q, N - 1))
            if state.L > 1:
                cands.update(divisors(state.L))
            valid = [d for d in cands if d > 0 and pow(a, d, N) == 1]
            if valid:
                r = minimize_order(a, N, min(valid))
                if r > 0 and r == true_r:
                    state.update(a, r)
                    if state.L == lam:
                        K = k
                        break
        Ks.append(K)
    return sum(Ks) / len(Ks)


def measure_hybrid_K_noisy(
    N: int, d: int, noise_kwargs: dict, trials: int = 100, max_runs: int = 20,
    seed: int = 0,
) -> float:
    """노이즈 하의 hybrid 알고리즘의 인수까지 runs. (F2 검증)"""
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
                cands = set(convergent_denominators(ki, Q, N - 1))
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
                if K_found < max_runs:
                    break
            if K_found < max_runs:
                break
            # factor_from_exponent (L 짝수일 때만 의미 — Regev setup 에서는 보통 X)
            if state.L > 1:
                rng_f = random.Random(t)
                res = factor_from_exponent(N, state.L, rng_f, max_attempts=5)
                if res and 1 < res.factor < N:
                    K_found = K
                    break
        Ks.append(K_found)
    return sum(Ks) / len(Ks)


def f1_search(N: int = 437, trials: int = 50):
    """F1: 노이즈 수준에 따른 K_λ 변화."""
    print(f"\n# F1: K_λ vs noise (N={N}, {trials} trials)")
    print(f"  baseline (noise-free): ", end="", flush=True)
    base = measure_klambda_alg_noisy(N, {}, trials=trials)
    print(f"{base:.2f}")

    print(f"\n  phase_sigma σ:")
    for sigma in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
        K = measure_klambda_alg_noisy(N, {"phase_sigma": sigma}, trials=trials)
        marker = " ← MIN!" if K < base * 0.95 else ""
        print(f"    σ={sigma:>4.2f}: K_λ = {K:>6.2f}{marker}", flush=True)

    print(f"\n  depolarizing p:")
    for p in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]:
        K = measure_klambda_alg_noisy(N, {"depolarizing": p}, trials=trials)
        marker = " ← MIN!" if K < base * 0.95 else ""
        print(f"    p={p:>4.2f}: K_λ = {K:>6.2f}{marker}", flush=True)


def f2_search(N: int = 437, d: int = 4, trials: int = 50):
    """F2: 노이즈 수준에 따른 hybrid K 변화 (U-shape 검색)."""
    print(f"\n# F2: hybrid K vs noise (N={N}, d={d}, {trials} trials)")
    print(f"  baseline (noise-free): ", end="", flush=True)
    base = measure_hybrid_K_noisy(N, d, {}, trials=trials)
    print(f"{base:.2f}")

    print(f"\n  phase_sigma σ:")
    for sigma in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
        K = measure_hybrid_K_noisy(N, d, {"phase_sigma": sigma}, trials=trials)
        marker = " ← MIN!" if K < base * 0.95 else ""
        print(f"    σ={sigma:>4.2f}: hybrid K = {K:>6.2f}{marker}", flush=True)

    print(f"\n  depolarizing p:")
    for p in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]:
        K = measure_hybrid_K_noisy(N, d, {"depolarizing": p}, trials=trials)
        marker = " ← MIN!" if K < base * 0.95 else ""
        print(f"    p={p:>4.2f}: hybrid K = {K:>6.2f}{marker}", flush=True)

    print(f"\n  amplitude_damp γ:")
    for gamma in [0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01]:
        K = measure_hybrid_K_noisy(N, d, {"amplitude_damp": gamma}, trials=trials)
        marker = " ← MIN!" if K < base * 0.95 else ""
        print(f"    γ={gamma:>6.4f}: hybrid K = {K:>6.2f}{marker}", flush=True)


def main(argv):
    if "--f1" in argv:
        f1_search()
        return
    if "--f2" in argv:
        f2_search()
        return
    # 둘 다
    f1_search()
    f2_search()


if __name__ == "__main__":
    main(sys.argv)
