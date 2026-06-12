"""
노이즈 하의 알고리즘 K_λ (실제 다중-base 알고리즘) 측정 (Phase 2 정리 3 검증).

K_λ^alg(η) := 노이즈 η 하에서 (C) 알고리즘이 L=λ(N) 에 도달할 때까지 뽑은 base 수.

정리 3 예측: E[K_λ^alg(η)] ≤ E[K_λ^ideal] / g_M(η)
where g_M(η) 는 experiments/g_eta.py 가 측정한 worst-case 회수율.

실행:
    python -m experiments.k_lambda_alg
"""

from __future__ import annotations
import math
import random
import sys
from dataclasses import dataclass

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
)
from noise import simulate_period_finding_noisy


def run_algorithm_to_lambda(
    N: int, lam: int, noise_kwargs: dict,
    max_bases: int = 200, seed: int = 0,
) -> int:
    """(C) 알고리즘을 L = λ(N) 도달까지 실행, 사용한 base 수 반환."""
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    state = MultiBaseState()

    for K in range(1, max_bases + 1):
        for _retry in range(50):
            a = rng_py.randrange(2, N)
            if math.gcd(a, N) == 1:
                break
        else:
            continue

        true_r = classical_order(a, N)

        # Fast path: r_a | L_before (covered) — 회수 결정적 (정리 1)
        if state.L > 1 and pow(a, state.L, N) == 1:
            state.update(a, true_r)
            if state.L == lam:
                return K
            continue

        # Slow path: 측정 + (C) 후처리
        m = simulate_period_finding_noisy(a, N, rng=rng_np, **noise_kwargs)
        cands = set(convergent_denominators(m.k, m.Q, N - 1))
        if state.L > 1:
            cands.update(divisors(state.L))
        valid_ds = [d for d in cands if d > 0 and pow(a, d, N) == 1]
        if not valid_ds:
            continue
        r_recovered = minimize_order(a, N, min(valid_ds))
        if r_recovered == true_r:
            state.update(a, true_r)
            if state.L == lam:
                return K

    return max_bases


def measure_k_lambda_alg(
    N: int, noise_kwargs: dict, trials: int = 100, seed: int = 0,
) -> tuple[float, int]:
    """K_λ^alg 의 평균과 max."""
    lam = math.lcm(*[p - 1 for p in _factor_semiprime(N)])
    samples = [
        run_algorithm_to_lambda(N, lam, noise_kwargs, seed=seed + t * 1000)
        for t in range(trials)
    ]
    return sum(samples) / len(samples), max(samples)


def _factor_semiprime(N: int) -> tuple[int, int]:
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            return p, N // p
    raise ValueError(f"{N} 가 반소수가 아님")


def main(argv):
    Ns = [437] if len(argv) <= 1 else [int(x) for x in argv[1:]]

    noise_setups = [
        ("noise-free",       {}),
        ("depol p=0.1",      {"depolarizing": 0.1}),
        ("depol p=0.3",      {"depolarizing": 0.3}),
        ("depol p=0.5",      {"depolarizing": 0.5}),
        ("depol p=0.7",      {"depolarizing": 0.7}),
        ("bias_zero p=0.5",  {"bias_zero": 0.5}),
        ("modexp q=0.3",     {"modexp_error": 0.3}),
        ("modexp q=0.5",     {"modexp_error": 0.5}),
        ("phase σ=1.0",      {"phase_sigma": 1.0}),
        ("phase σ=2.0",      {"phase_sigma": 2.0}),
    ]
    trials = 100

    for N in Ns:
        print(f"\n# K_λ^alg 측정: N={N}, {trials} trials")
        print(f"  {'노이즈':<20} {'mean':>6} {'max':>4}  {'mean/baseline':>14}")
        baseline = None
        for label, kwargs in noise_setups:
            mean, mx = measure_k_lambda_alg(N, kwargs, trials=trials)
            if baseline is None:
                baseline = mean
            ratio = mean / baseline if baseline > 0 else 0
            print(f"  {label:<20} {mean:>6.2f} {mx:>4}  {ratio:>14.3f}")


if __name__ == "__main__":
    main(sys.argv)
