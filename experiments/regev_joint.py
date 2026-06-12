"""
Joint-constrained Regev-like 측정 시뮬레이션 — Phase 5 본실행.

Regev 2023 의 측정은 좌표별 독립이 아니라 *joint linear constraint*
`Σ_i b_i · k_i ≈ 0 (mod r)` 를 만족 (b_i 는 quadratic character 정보).

본 모듈은 이 제약을 numpy 로 단순화 모델링:
1. d 개 독립 Shor 측정 k_1, ..., k_d 샘플링.
2. 제약 `Σ b_i k_i ≡ 0 (mod target_r)` 에 가장 가까운 점으로 *직교 projection*.

이 분포로부터 (C) 좌표별 후처리 적용. 독립 Shor (regev_c.py) 대비 회수율 차이 측정.

실행:
    python -m experiments.regev_joint
"""

from __future__ import annotations
import math
import random
import sys

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    order_from_measurement,
)
from shor import simulate_period_finding
from noise import simulate_period_finding_noisy


def lambda_semiprime(N: int) -> int:
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            return math.lcm(p - 1, N // p - 1)
    raise ValueError


def regev_joint_sample(
    bases: list[int], N: int, Q: int, b_vec: list[int], target_r: int,
    rng_np: np.random.Generator, noise_kwargs: dict | None = None,
) -> list[int]:
    """d 개 base 의 joint-constrained 측정 샘플.

    1. 각 base 에서 독립 Shor 측정 k_i.
    2. delta = (Σ b_i k_i) mod target_r 계산 (signed → [-target_r/2, target_r/2)).
    3. 제약 만족하도록 최소 ℓ² 거리 projection: k_i ← k_i − delta · b_i / Σ b_j².

    반환: 제약 만족하는 (k_1, ..., k_d) ∈ [0, Q)^d.
    """
    raw = []
    for a in bases:
        if noise_kwargs:
            m = simulate_period_finding_noisy(a, N, rng=rng_np, **noise_kwargs)
        else:
            m = simulate_period_finding(a, N, rng=rng_np)
        raw.append(m.k)

    deviation = sum(b * k for b, k in zip(b_vec, raw)) % target_r
    if deviation > target_r // 2:
        deviation -= target_r

    norm_b_sq = sum(b * b for b in b_vec)
    if norm_b_sq == 0:
        return raw

    corrected = []
    for k_i, b_i in zip(raw, b_vec):
        delta_k = -deviation * b_i / norm_b_sq
        corrected.append(int(round(k_i + delta_k)) % Q)
    return corrected


def regev_joint_runs_to_lambda(
    N: int, d: int, max_runs: int = 50, seed: int = 0,
    noise_kwargs: dict | None = None,
) -> int:
    """Joint-constrained Regev + (C) 좌표별 후처리. K_λ 측정 (in runs)."""
    lam = lambda_semiprime(N)
    state = MultiBaseState()
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    Q = 1 << (2 * max(1, (N - 1).bit_length()))

    for K in range(1, max_runs + 1):
        # d 개 base 와 b 계수 무작위 선택
        bases = []
        while len(bases) < d:
            a = rng_py.randrange(2, N)
            if math.gcd(a, N) == 1:
                bases.append(a)
        b_vec = [rng_py.randrange(1, N) for _ in range(d)]

        ks = regev_joint_sample(
            bases, N, Q, b_vec, lam, rng_np, noise_kwargs=noise_kwargs,
        )

        # 좌표별 (C) 후처리
        for a, k in zip(bases, ks):
            d_recovered = order_from_measurement(a, N, k, Q, state)
            if d_recovered > 0:
                state.update(a, d_recovered)
                if state.L == lam:
                    return K
    return max_runs


def regev_independent_runs_to_lambda(
    N: int, d: int, max_runs: int = 50, seed: int = 0,
) -> int:
    """(비교 baseline) 독립 Shor 측정 (regev_c.py 와 동일 방식)."""
    lam = lambda_semiprime(N)
    state = MultiBaseState()
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    for K in range(1, max_runs + 1):
        for _ in range(d):
            for _retry in range(50):
                a = rng_py.randrange(2, N)
                if math.gcd(a, N) == 1:
                    break
            else:
                continue
            m = simulate_period_finding(a, N, rng=rng_np)
            d_rec = order_from_measurement(a, N, m.k, m.Q, state)
            if d_rec > 0:
                state.update(a, d_rec)
                if state.L == lam:
                    return K
    return max_runs


def compare(N: int, trials: int = 200):
    n_bits = math.ceil(math.log2(N))
    d = max(1, int(math.sqrt(n_bits + 4) + 0.5))

    indep = [regev_independent_runs_to_lambda(N, d, seed=t * 1000)
             for t in range(trials)]
    joint = [regev_joint_runs_to_lambda(N, d, seed=t * 1000)
             for t in range(trials)]

    return {
        "N": N, "d": d, "trials": trials,
        "indep_mean": sum(indep) / trials,
        "joint_mean": sum(joint) / trials,
        "indep_max": max(indep),
        "joint_max": max(joint),
        "ratio": (sum(joint) / sum(indep)) if sum(indep) > 0 else float("nan"),
    }


def measure_noisy(N: int, d: int, noise_kwargs: dict, trials: int = 100):
    """Joint-constrained + 노이즈 하의 (C) 좌표별 K_λ."""
    samples = [
        regev_joint_runs_to_lambda(
            N, d, seed=t * 1000, noise_kwargs=noise_kwargs,
        )
        for t in range(trials)
    ]
    return sum(samples) / len(samples), max(samples)


def main(argv):
    if "--noise" in argv:
        # Joint + 노이즈 견고함 측정
        N = int(argv[argv.index("--noise") + 1]) if len(argv) > argv.index("--noise") + 1 else 437
        n_bits = math.ceil(math.log2(N))
        d = max(1, int(math.sqrt(n_bits + 4) + 0.5))

        noise_setups = [
            ("noise-free",       {}),
            ("depol p=0.3",      {"depolarizing": 0.3}),
            ("depol p=0.5",      {"depolarizing": 0.5}),
            ("modexp q=0.3",     {"modexp_error": 0.3}),
            ("phase σ=1.0",      {"phase_sigma": 1.0}),
        ]
        print(f"# Joint Regev + (C) 좌표별 — 노이즈 견고함, N={N}, d={d}, 100 trials\n")
        print(f"  {'noise':<18}  {'mean runs':>10}  {'max':>4}")
        for label, kwargs in noise_setups:
            mean, mx = measure_noisy(N, d, kwargs)
            print(f"  {label:<18}  {mean:>10.2f}  {mx:>4}")
        return

    Ns = [int(x) for x in argv[1:]] if len(argv) > 1 else [77, 143, 437, 1147]
    trials = 200

    print(f"# Joint-constrained Regev + (C) 좌표별 vs 독립 Shor + (C) 좌표별")
    print(f"# trials = {trials}, noise-free, 제약 = Σ b_i k_i ≡ 0 mod λ(N)\n")
    print(f"  {'N':>5}  {'d':>2}  {'indep mean':>10}  {'joint mean':>10}  "
          f"{'ratio':>6}  {'indep/joint max':>16}")

    for N in Ns:
        r = compare(N, trials=trials)
        print(f"  {N:>5}  {r['d']:>2}  {r['indep_mean']:>10.2f}  "
              f"{r['joint_mean']:>10.2f}  {r['ratio']:>6.3f}  "
              f"{r['indep_max']:>7d}/{r['joint_max']:<8d}")


if __name__ == "__main__":
    main(sys.argv)
