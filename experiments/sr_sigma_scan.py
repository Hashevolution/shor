"""
H12 검증 — N 별 σ_opt 추적.

가설 H12c (사용자 직관):
  "작은 자물쇠는 열쇠를 작게 흔들어가며 맞춰주고,
   큰 자물쇠는 조금더 크게 흔들어가며 맞춰주고,
   각 자물쇠의 크기에 비례한 흔들림이 열쇠를 잘맞게 해준다."

= σ_opt 가 N 에 따라 증가하는가?

실험: N=437, 1147 각각의 fine σ grid scan.

실행:
    python -m experiments.sr_sigma_scan
"""

from __future__ import annotations
import math
import sys

from experiments.sr_amplify import measure_hybrid_extended


def sigma_scan(N: int, d: int, sigmas: list[float], trials: int = 200,
               seed: int = 0):
    """주어진 (N, d) 에서 σ scan, σ_opt 식별."""
    print(f"\n# σ scan: N={N}, d={d}, {trials} trials/σ")
    print(f"  {'σ':>7}  {'hybrid K':>10}  {'SR vs σ=0':>10}")

    baseline = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
    print(f"  {0.0:>7.3f}  {baseline:>10.3f}  {'baseline':>10}", flush=True)

    results = {0.0: baseline}
    for sigma in sigmas:
        K = measure_hybrid_extended(
            N, d, {"phase_sigma": sigma}, trials=trials, seed=seed,
        )
        sr = (baseline - K) / baseline * 100 if baseline > 0 else 0
        marker = ""
        if sr > 5:
            marker = " ★ strong SR"
        elif sr > 1:
            marker = " ↓ SR"
        elif sr < -1:
            marker = " ↑ anti-SR"
        print(f"  {sigma:>7.3f}  {K:>10.3f}  {sr:>+9.2f}%{marker}", flush=True)
        results[sigma] = K

    # Find σ_opt
    min_K = min(results.values())
    sigma_opt = [s for s, k in results.items() if k == min_K][0]
    SR_max = (baseline - min_K) / baseline * 100
    print(f"\n  σ_opt({N}, d={d}) = {sigma_opt:.3f}")
    print(f"  SR_max          = {SR_max:.2f}%")
    return sigma_opt, SR_max


def main(argv):
    # N=437 fine scan
    sigmas_small = [0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50]
    s_437_4, sr_437_4 = sigma_scan(N=437, d=4, sigmas=sigmas_small, trials=300)

    # N=1147 fine scan (H12c 의 핵심 검증)
    sigmas_med = [0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 1.0]
    s_1147_2, sr_1147_2 = sigma_scan(N=1147, d=2, sigmas=sigmas_med, trials=150)

    # 비교
    print(f"\n\n=== H12c 검증 ===")
    print(f"  σ_opt(437, d=4)   = {s_437_4:.3f}")
    print(f"  σ_opt(1147, d=2)  = {s_1147_2:.3f}")
    print(f"  σ ratio           = {s_1147_2 / s_437_4 if s_437_4 > 0 else 0:.2f}")
    print(f"  N ratio           = {1147/437:.2f}")

    if s_1147_2 > s_437_4 * 1.5:
        print(f"  → **H12c 지지** — σ_opt 가 N 따라 *명확히 증가*")
        print(f"    사용자 직관 (자물쇠 크기 ∝ 흔들림 크기) 확정")
    elif s_1147_2 > s_437_4 * 1.1:
        print(f"  → H12 약 지지 — 증가는 있으나 polynomial 인지는 불확실")
    elif abs(s_1147_2 - s_437_4) < 0.02:
        print(f"  → H12 기각, σ_opt 거의 일정 — H1-H4 가 충분 설명")
    elif s_1147_2 < s_437_4 * 0.7:
        print(f"  → H12b 지지 — σ_opt 가 N 따라 *감소*")
    else:
        print(f"  → 결과 모호, 더 많은 N 또는 더 많은 trials 필요")


if __name__ == "__main__":
    main(sys.argv)
