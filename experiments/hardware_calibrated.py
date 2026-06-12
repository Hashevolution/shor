"""
IBM Quantum Eagle (127q) 하드웨어 사양에 맞춘 (C) 시뮬레이션 — Phase 4 (hardware proxy).

본 실험은 실제 hardware run 이 아니라, 공개된 IBM Quantum 디바이스 평균 사양을
우리 noise.py 5종 노이즈 모델로 매핑한 numpy 시뮬레이션. 실제 hardware 실행은 qiskit
+ IBM Q 계정이 필요 (본 저장소 정책 외).

IBM Eagle 127q 평균 사양 (2024-2025):
- T1 ≈ 150 μs (relaxation)
- T2 ≈ 100 μs (dephasing)
- 1-qubit gate error ≈ 0.03 %
- 2-qubit (ECR) gate error ≈ 1 %
- Readout error ≈ 2 % per bit

N=15 Shor 회로 추정:
- 4 work + 8 counting = 12 qubits
- 회로 깊이 ≈ 50 μs, ~30 gates 중 ~10 two-qubit
- t = 8, Q = 256

mapping 으로 noise.py 파라미터:
- readout_flip p = 0.02   (per-bit readout error, 그대로)
- modexp_error q = 0.05   (10 two-qubit × 1% ≈ 10% cumulative; 부분 corruption q=0.05 사용)
- amplitude_damp γ = 0.002 (50μs / 150μs T1 → 약 30% decay over Q=256; γ·Q ≈ 0.5)
- phase_sigma σ = 0.3      (50μs / 100μs T2 → 약 σ ≈ 0.3 rad equivalent dephasing)
- depolarizing p = 0.01    (residual gate errors as effective uniform noise)

실행:
    python -m experiments.hardware_calibrated
"""

from __future__ import annotations
import math
import random
import sys

import numpy as np

from demo import verify_c_determinism
from experiments.k_lambda_alg import measure_k_lambda_alg


# Eagle 127q 평균 사양 → 우리 노이즈 모델로 매핑
HARDWARE_CALIBRATED = {
    "readout_flip": 0.02,
    "modexp_error": 0.05,
    "amplitude_damp": 0.002,
    "phase_sigma": 0.3,
    "depolarizing": 0.01,
}


def main(argv):
    print("# IBM Eagle 127q 하드웨어 사양 매핑 — N=15 Shor 시뮬")
    print(f"# 노이즈 (5종 동시): {HARDWARE_CALIBRATED}")
    print()

    # (1) Theorem 1 검증 (covered, violations)
    print("## Theorem 1: covered / violations 측정 (N=15, 500 trials)")
    verify_c_determinism(
        N=15, trials=500,
        noise_subset=None,  # 우리는 mixed setup 을 별도로 추가
    )
    print()

    # (2) Mixed setup 으로 직접 verify_c_determinism
    print("## Theorem 1: hardware-calibrated mixed noise (N=15, 500 trials)")
    print("   (verify_c_determinism 의 mixed setup 직접 호출)")
    from fractions import Fraction
    from classical import classical_order
    from multi_base import (
        MultiBaseState, convergent_denominators, divisors, minimize_order,
    )
    from noise import simulate_period_finding_noisy

    state = MultiBaseState()
    rng_py = random.Random(0)
    rng_np = np.random.default_rng(0)
    covered = violations = lucky = missed = success_total = 0
    N = 15
    trials = 500

    for _ in range(trials):
        for _retry in range(50):
            a = rng_py.randrange(2, N)
            if math.gcd(a, N) == 1:
                break
        else:
            continue

        true_r = classical_order(a, N)
        L_before = state.L
        cond = (L_before > 1) and (L_before % true_r == 0)

        m = simulate_period_finding_noisy(a, N, rng=rng_np, **HARDWARE_CALIBRATED)

        cands = set(convergent_denominators(m.k, m.Q, N - 1))
        if L_before > 1:
            cands.update(divisors(L_before))
        valid = [d for d in cands if d > 0 and pow(a, d, N) == 1]
        success = bool(valid) and minimize_order(a, N, min(valid)) == true_r

        if cond:
            covered += 1
            if not success:
                violations += 1
        else:
            if success:
                lucky += 1
            else:
                missed += 1
        if success:
            state.update(a, true_r)
            success_total += 1

    print(f"  covered={covered}, violations={violations}, lucky={lucky}, missed={missed}")
    print(f"  success rate: {success_total/trials:.1%}")
    print()

    # (3) Theorem 3 검증: K_λ^alg under hardware-calibrated noise
    print("## Theorem 3: K_λ^alg (N=15, 100 trials)")
    mean, mx = measure_k_lambda_alg(15, HARDWARE_CALIBRATED, trials=100)
    print(f"  mean K_λ^alg = {mean:.2f}, max = {mx}")

    # baseline
    mean_clean, mx_clean = measure_k_lambda_alg(15, {}, trials=100)
    print(f"  baseline (noise-free): mean = {mean_clean:.2f}, max = {mx_clean}")
    print(f"  ratio = {mean / mean_clean:.2f}x")


if __name__ == "__main__":
    main(sys.argv)
