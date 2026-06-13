"""
N=1147 d=1 의 *빠른* preview confirm (trials=100).

각 seed 약 5-10분. 5 seeds × 10분 = ~50분.
trials 300 의 더 작은 변동성 가지나 빠른 신호 확인 가능.

실행:
    python -u -m experiments.confirm_1147_d1_quick
"""

from __future__ import annotations
import math
import statistics
import sys

from experiments.sr_amplify import measure_hybrid_extended


def main(argv):
    N = 1147
    d = 1
    trials = 100  # 빠르게
    n_seeds = 5

    print(f"# N=1147 d=1 quick confirm ({trials} trials × {n_seeds} seeds)")
    print(f"# 예상 시간: 각 seed ~5-10분, 총 30-50분\n", flush=True)
    print(f"  {'seed':>5}  {'K(σ=0)':>10}  {'K(σ=.05)':>10}  {'SR %':>8}", flush=True)

    srs = []
    K0s = []
    for seed in range(1, n_seeds + 1):
        K0 = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
        K1 = measure_hybrid_extended(N, d, {"phase_sigma": 0.05}, trials=trials, seed=seed)
        sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
        srs.append(sr)
        K0s.append(K0)
        marker = " ★" if sr > 2 else (" ↓" if sr > 0.5 else "")
        print(f"  {seed:>5}  {K0:>10.3f}  {K1:>10.3f}  {sr:>+7.2f}%{marker}",
              flush=True)

    print()
    mean_sr = statistics.mean(srs)
    sd_sr = statistics.stdev(srs) if len(srs) > 1 else 0
    print(f"  Mean SR = {mean_sr:+.2f}%  (sd={sd_sr:.2f})")
    print(f"  Mean K_base = {statistics.mean(K0s):.3f}")
    n_positive = sum(1 for sr in srs if sr > 0)
    print(f"  Positive: {n_positive}/{n_seeds}")

    if mean_sr > 1.5 and n_positive >= 4:
        print(f"\n  ★ Quick preview 지지 — 정식 300 trials 진행 권고")
    elif mean_sr > 0.5:
        print(f"\n  weak — 정식 trials 후 결정")
    else:
        print(f"\n  ✗ peak 2.60% 의심 — fluke 가능")


if __name__ == "__main__":
    main(sys.argv)
