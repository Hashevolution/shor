"""
N=1147 d=1 의 2.60% SR peak confirm.

5 seeds × 300 trials × σ ∈ {0, 0.05}
예상 시간: 30-45분

실행:
    python -m experiments.confirm_1147_d1
"""

from __future__ import annotations
import math
import statistics
import sys

from experiments.sr_amplify import measure_hybrid_extended


def main(argv):
    N = 1147
    d = 1
    trials = 300
    n_seeds = 5

    print(f"# N=1147 d=1 의 2.60% peak confirm")
    print(f"# {n_seeds} seeds × {trials} trials × σ ∈ {{0, 0.05}}")
    print(f"# 예상 시간: 30-45분\n")
    print(f"  {'seed':>5}  {'K(σ=0)':>10}  {'K(σ=.05)':>10}  {'SR %':>8}  {'Δ':>8}")

    srs = []
    deltas = []
    K0s = []
    K1s = []
    for seed in range(1, n_seeds + 1):
        K0 = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
        K1 = measure_hybrid_extended(N, d, {"phase_sigma": 0.05}, trials=trials, seed=seed)
        sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
        delta = K0 - K1
        srs.append(sr)
        deltas.append(delta)
        K0s.append(K0)
        K1s.append(K1)
        marker = ""
        if sr > 3:
            marker = " ★ strong"
        elif sr > 1:
            marker = " ↓ SR"
        elif sr < -1:
            marker = " ↑ anti"
        print(f"  {seed:>5}  {K0:>10.3f}  {K1:>10.3f}  {sr:>+7.2f}%  {delta:>+8.3f}{marker}",
              flush=True)

    print()
    mean_sr = statistics.mean(srs)
    sd_sr = statistics.stdev(srs) if len(srs) > 1 else 0
    se_sr = sd_sr / math.sqrt(len(srs)) if len(srs) > 1 else 0
    t = mean_sr / se_sr if se_sr > 0 else float('nan')
    mean_K0 = statistics.mean(K0s)
    mean_K1 = statistics.mean(K1s)

    print(f"  Mean K(σ=0)  = {mean_K0:.3f}")
    print(f"  Mean K(σ=.05) = {mean_K1:.3f}")
    print(f"  Mean SR  = {mean_sr:+.2f}%  (sd={sd_sr:.2f}, SE={se_sr:.2f}, t={t:.2f})")
    print(f"  Mean Δ   = {statistics.mean(deltas):+.4f}")

    n_positive = sum(1 for sr in srs if sr > 0)
    print(f"  Positive direction: {n_positive}/{n_seeds}")

    # Confirm 판정
    print(f"\n=== Confirm 판정 ===")
    if mean_sr > 2 and t > 2:
        print(f"  ★★★ STRONGLY CONFIRMED — N=1147 d=1 SR ~{mean_sr:.1f}% real, t={t:.1f}")
        print(f"  → AOP 의 peak 진짜")
    elif mean_sr > 1.5:
        print(f"  ★★ confirmed — substantial SR ({mean_sr:.1f}%)")
        print(f"  → N=1147 d=1 가 sub-optimal peak")
    elif mean_sr > 0.5:
        print(f"  ★ weak SR — direction 일치하나 작음")
    else:
        print(f"  ✗ 2.60% 가 fluke 가능")


if __name__ == "__main__":
    main(sys.argv)
