"""
노이즈 모델별 효과적 추출 확률 g(η) 측정 (Phase 2).

정의: g_M(η) := P[(C) 가 새 base a 의 위수 r_a 를 회복 | r_a ∤ L_before, 노이즈 η, 모델 M].

방법:
  · 다양한 L_before (= λ(N) 의 약수) 를 강제로 세팅
  · r_a ∤ L_before 인 base a 만 선택
  · 노이즈 η 하에서 (C) 적용, 회복 성공률 측정
  · η = 0 의 g_0 (baseline) 와 비교

g_M(η) / g_0 의 (1-p) 식 가설 검증.

실행:
    python -m experiments.g_eta
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


def coprime_bases_with_r_not_dividing(N: int, L: int, rng: random.Random,
                                       max_tries: int = 200) -> int | None:
    """gcd(a, N)=1 이고 r_a ∤ L 인 base a 무작위 선택."""
    for _ in range(max_tries):
        a = rng.randrange(2, N)
        if math.gcd(a, N) != 1:
            continue
        if L > 1 and L % classical_order(a, N) == 0:
            continue  # r_a | L: 추출과 무관, skip
        return a
    return None


def measure_g(
    N: int, L_set: int, kwargs: dict, trials: int = 500, seed: int = 0,
) -> float:
    """주어진 L_before = L_set 에서 노이즈 kwargs 하의 g 측정."""
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    successes = 0
    valid = 0

    for _ in range(trials):
        a = coprime_bases_with_r_not_dividing(N, L_set, rng_py)
        if a is None:
            continue
        true_r = classical_order(a, N)
        m = simulate_period_finding_noisy(a, N, rng=rng_np, **kwargs)

        cands = set(convergent_denominators(m.k, m.Q, N - 1))
        if L_set > 1:
            cands.update(divisors(L_set))
        valid_ds = [d for d in cands if d > 0 and pow(a, d, N) == 1]
        success = (bool(valid_ds)
                   and minimize_order(a, N, min(valid_ds)) == true_r)

        valid += 1
        if success:
            successes += 1

    return successes / valid if valid > 0 else 0.0


def main(argv):
    N = int(argv[1]) if len(argv) > 1 else 437
    L_set = int(argv[2]) if len(argv) > 2 else 1  # 초기 L=1 (가장 어려운 조건)

    print(f"# g(η) 측정: N={N}, L_before={L_set}")
    print(f"# trials=500 per condition\n")

    # baseline g_0
    g0 = measure_g(N, L_set, {}, trials=500)
    print(f"  g_0 (noise-free): {g0:.4f}")

    # depolarizing curves
    print("\n  depolarizing p:")
    for p in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]:
        g = measure_g(N, L_set, {"depolarizing": p}, trials=500)
        ratio = g / g0 if g0 > 0 else 0
        predicted = 1 - p  # hypothesis
        print(f"    p={p:.2f}: g={g:.4f}, g/g_0={ratio:.4f}, (1-p)={predicted:.4f}")

    # bias_zero
    print("\n  bias_zero p:")
    for p in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]:
        g = measure_g(N, L_set, {"bias_zero": p}, trials=500)
        ratio = g / g0 if g0 > 0 else 0
        predicted = 1 - p
        print(f"    p={p:.2f}: g={g:.4f}, g/g_0={ratio:.4f}, (1-p)={predicted:.4f}")

    # modexp_error
    print("\n  modexp_error q:")
    for q in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]:
        g = measure_g(N, L_set, {"modexp_error": q}, trials=500)
        ratio = g / g0 if g0 > 0 else 0
        predicted = 1 - q
        print(f"    q={q:.2f}: g={g:.4f}, g/g_0={ratio:.4f}, (1-q)={predicted:.4f}")

    # phase_sigma
    print("\n  phase_sigma σ:")
    for sigma in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        g = measure_g(N, L_set, {"phase_sigma": sigma}, trials=500)
        ratio = g / g0 if g0 > 0 else 0
        print(f"    σ={sigma:.2f}: g={g:.4f}, g/g_0={ratio:.4f}")

    # amplitude_damp
    print("\n  amplitude_damp γ:")
    for gamma in [0.0, 0.001, 0.005, 0.01, 0.05, 0.1]:
        g = measure_g(N, L_set, {"amplitude_damp": gamma}, trials=500)
        ratio = g / g0 if g0 > 0 else 0
        print(f"    γ={gamma:.4f}: g={g:.4f}, g/g_0={ratio:.4f}")


if __name__ == "__main__":
    main(sys.argv)
