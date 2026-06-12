"""
AOP (Anti-Optimization Principle) 검증.

가설: SR 발현 = algorithmic slack (K_baseline > floor) 만 있으면 됨.
N 자체 가 아니라 *slack* 이 결정.

실험: N=2491 (이전 SR=0%) 에서 d=1 (slack 재생) 시 SR 발현?
  Yes → AOP 확정 — slack 만 있으면 N 무관 SR
  No  → 다른 mechanism (N-specific) 필요

추가: N ∈ {437, 1147, 2491} × d ∈ {1, 2, 3} grid 로 SR vs K_baseline 정리.

실행:
    python -m experiments.sr_aop
"""

from __future__ import annotations
import math
import sys

from experiments.sr_amplify import measure_hybrid_extended


def slack_test(N: int, d: int, trials: int = 300, seed: int = 0):
    """주어진 (N, d) 에서 K_baseline + SR 측정."""
    K0 = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
    K1 = measure_hybrid_extended(N, d, {"phase_sigma": 0.05}, trials=trials, seed=seed)
    sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
    slack = K0 - 1.0  # "room above floor"
    return K0, K1, sr, slack


def main(argv):
    print("# AOP (Anti-Optimization Principle) 검증")
    print("# 가설: SR 발현 = K_baseline > 1 (slack) 만 있으면, N 무관\n")

    # Grid: N × d
    Ns = [437, 1147, 2491]
    ds = [1, 2, 3]
    trials = 300

    print(f"  {'N':>5} {'d':>2}  {'K_base':>8}  {'K_noise':>8}  {'SR %':>8}  "
          f"{'slack':>7}")

    results = []
    for N in Ns:
        for d in ds:
            try:
                K0, K1, sr, slack = slack_test(N, d, trials=trials)
                results.append((N, d, K0, K1, sr, slack))
                marker = ""
                if sr > 3:
                    marker = " ★ strong SR"
                elif sr > 1:
                    marker = " ↓ SR"
                elif sr < -1:
                    marker = " ↑ anti-SR"
                print(f"  {N:>5} {d:>2}  {K0:>8.3f}  {K1:>8.3f}  "
                      f"{sr:>+7.2f}%  {slack:>+6.2f}{marker}", flush=True)
            except Exception as e:
                print(f"  {N:>5} {d:>2}  ERROR: {e}", flush=True)

    # AOP 평가
    print(f"\n=== AOP 평가 ===")
    print(f"  가설: SR 가 slack (K_baseline-1) 단조 함수 인가?")
    print()

    # 모든 (N, d) 의 (slack, SR) 정렬
    sorted_results = sorted(results, key=lambda r: r[5])
    print(f"  slack 순:")
    for N, d, K0, K1, sr, slack in sorted_results:
        marker = " ✓ SR" if sr > 1 else (" ✗ no SR" if abs(sr) < 1 else " anti")
        print(f"    slack={slack:>5.2f}  N={N:>5} d={d:>2}  SR={sr:>+6.2f}%{marker}",
              flush=True)

    # 핵심 검정: N=2491 d=1 에서 SR > 0?
    n2491_d1 = [r for r in results if r[0] == 2491 and r[1] == 1]
    if n2491_d1:
        _, _, _, _, sr_2491_d1, slack_2491_d1 = n2491_d1[0]
        print()
        print(f"  *결정적*: N=2491 d=1 의 slack={slack_2491_d1:.2f}, SR={sr_2491_d1:+.2f}%")
        if sr_2491_d1 > 2:
            print(f"  → **AOP 확정** — slack 만 있으면 큰 N 에서도 SR 발현")
            print(f"     'Anti-optimization' framing = clean conceptual finding")
        elif sr_2491_d1 > 0.5:
            print(f"  → AOP 부분 지지 — small SR 있으나 작음")
        else:
            print(f"  → AOP 기각 — slack 만으로 SR 발현 안 함. N-specific mechanism 필요")


if __name__ == "__main__":
    main(sys.argv)
