"""
Regev 2023 의 다중-base 회로 (단순화) + (C) 좌표별 후처리 검증.

가정 (Phase 5 정리 4): Regev 의 측정 (k_1, …, k_d) 의 *각 좌표 marginal* 이 Shor 의
단일 base 분포 `k_i ≈ j_i · Q / r_{a_i}` 와 동일. 좌표 간 결합 (joint correlation) 은
무시.

이 가정 하의 시뮬: 매 run 이 d 개 독립 Shor 측정. (C) 를 각 좌표에 적용 → 누적 L.

실행:
    python -m experiments.regev_c
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
    order_from_measurement,
)
from shor import simulate_period_finding


def lambda_semiprime(N: int) -> int:
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            return math.lcm(p - 1, N // p - 1)
    raise ValueError


def regev_runs_to_lambda(
    N: int, d: int, max_runs: int = 50, seed: int = 0,
) -> int:
    """(단순화) Regev: 매 run = d 개 독립 Shor 측정 + (C) 좌표별 후처리.

    L = λ(N) 에 도달할 때까지의 *run 수* 반환. K_λ^Regev-(C).
    """
    lam = lambda_semiprime(N)
    state = MultiBaseState()
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    for K in range(1, max_runs + 1):
        for _i in range(d):
            for _retry in range(50):
                a = rng_py.randrange(2, N)
                if math.gcd(a, N) == 1:
                    break
            else:
                continue
            m = simulate_period_finding(a, N, rng=rng_np)
            d_recovered = order_from_measurement(a, N, m.k, m.Q, state)
            if d_recovered > 0:
                state.update(a, d_recovered)
                if state.L == lam:
                    return K
    return max_runs


def measure(N: int, d: int, trials: int = 200, seed: int = 0) -> dict:
    samples = [regev_runs_to_lambda(N, d, seed=seed + t * 1000)
               for t in range(trials)]
    samples.sort()
    return {
        "N": N, "d": d, "trials": trials,
        "mean_runs": sum(samples) / len(samples),
        "mean_bases": sum(samples) / len(samples) * d,
        "p99_runs": samples[int(len(samples) * 0.99)],
        "max_runs": samples[-1],
    }


def main(argv):
    Ns = [int(x) for x in argv[1:]] if len(argv) > 1 else [77, 143, 437, 1147]

    print(f"# Regev-(C) 좌표별 후처리: K_λ^Regev (in runs) vs (in bases)")
    print(f"# d 선택: Regev 의 √(n+4) ≈ √(log₂ N) + 0.6 (소수 N 의 경우 ceil)")
    print(f"# trials = 200, noise-free\n")
    print(f"  {'N':>5}  {'d':>2}  {'runs (mean)':>12}  {'bases = runs·d':>15}  "
          f"{'p99(runs)':>10}  {'max(runs)':>10}")

    for N in Ns:
        n_bits = math.ceil(math.log2(N))
        d = max(1, int(math.sqrt(n_bits + 4) + 0.5))
        r = measure(N, d)
        print(f"  {N:>5}  {d:>2}  {r['mean_runs']:>12.2f}  "
              f"{r['mean_bases']:>15.2f}  {r['p99_runs']:>10}  {r['max_runs']:>10}")


if __name__ == "__main__":
    main(sys.argv)
