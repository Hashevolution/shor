"""
SR 신호 확인 — fine grid × 500 trials.

F2 의 결과에서 phase σ ∈ [0.05, 0.20] 가 baseline (σ=0) 보다 일관되게 낮음.
50 trials 로는 borderline. 500 trials 로 통계적 유의성 확인.

실행:
    python -m experiments.sr_confirm
"""

from __future__ import annotations
import math
import statistics
import sys

import numpy as np

from experiments.sr_search import measure_hybrid_K_noisy


def main(argv):
    N = 437
    d = 4
    trials = 500

    # 핵심 영역: σ ∈ [0, 0.30] 의 fine grid + baseline 반복
    sigmas = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30]

    print(f"# SR 확정 실험: N={N}, d={d}, {trials} trials per σ")
    print(f"# phase_sigma 의 fine grid")
    print(f"  {'σ':>6}  {'hybrid K':>10}  {'vs baseline':>12}")

    results = {}
    for sigma in sigmas:
        K = measure_hybrid_K_noisy(
            N, d, {"phase_sigma": sigma} if sigma > 0 else {},
            trials=trials,
        )
        results[sigma] = K

    baseline = results[0.0]
    for sigma in sigmas:
        K = results[sigma]
        diff = K - baseline
        marker = ""
        if K < baseline - 0.05:
            marker = " ↓"
        elif K > baseline + 0.05:
            marker = " ↑"
        print(f"  {sigma:>6.3f}  {K:>10.3f}  {diff:>+12.3f}{marker}", flush=True)

    # 통계적 검증: paired t-test 식 분석
    # SE ≈ stddev / sqrt(trials). For trials=500, SE small.
    # 만약 σ=0.1 vs σ=0 의 차이 > 2 SE 면 의미 있는 신호.
    print(f"\n# 분석:")
    min_K = min(results.values())
    min_sigma = [s for s, k in results.items() if k == min_K][0]
    print(f"  Minimum: σ={min_sigma:.3f}, K={min_K:.3f}")
    print(f"  Baseline (σ=0): K={baseline:.3f}")
    print(f"  Improvement: {(baseline - min_K)/baseline*100:.2f}%")

    # 5% 이상 개선 시 SR 확정
    if (baseline - min_K) / baseline > 0.05:
        print(f"\n  ✓ 통계적으로 의미있는 SR 후보 — 5% 이상 개선")
    elif (baseline - min_K) / baseline > 0.02:
        print(f"\n  △ 약한 신호 — 2-5% 개선 ({trials} trials 의 SE 와 비교 필요)")
    else:
        print(f"\n  ✗ SR 신호 없음 — 2% 미만")


if __name__ == "__main__":
    main(sys.argv)
