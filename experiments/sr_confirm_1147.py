"""
17.86% confirm — N=1147 d=2 의 σ=0.01 결과 1000 trials 로 확정.

σ scan 에서 150 trials/σ 로 17.86% 측정. 너무 적어서 sample fluke 가능성.
1000 trials × multiple seeds 로 통계적 확정.

실행:
    python -m experiments.sr_confirm_1147
"""

from __future__ import annotations
import math
import statistics
import sys

from experiments.sr_amplify import measure_hybrid_extended


def main(argv):
    N = 1147
    d = 2

    print(f"# 17.86% SR_max confirm: N={N}, d={d}")
    print(f"# σ=0.01 (σ_opt) 측정, 1000 trials × multiple seeds\n")
    print(f"  {'seed':>5}  {'σ=0':>10}  {'σ=0.01':>10}  {'SR %':>8}  {'Δ':>8}")

    n_seeds = 4
    srs = []
    deltas = []
    for seed in range(1, n_seeds + 1):
        K0 = measure_hybrid_extended(N, d, {}, trials=1000, seed=seed)
        K1 = measure_hybrid_extended(N, d, {"phase_sigma": 0.01}, trials=1000, seed=seed)
        sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
        delta = K0 - K1
        srs.append(sr)
        deltas.append(delta)
        print(f"  {seed:>5}  {K0:>10.3f}  {K1:>10.3f}  {sr:>+7.2f}%  {delta:>+8.3f}",
              flush=True)

    print()
    mean_sr = statistics.mean(srs)
    sd_sr = statistics.stdev(srs) if len(srs) > 1 else 0
    se_sr = sd_sr / math.sqrt(len(srs)) if len(srs) > 1 else 0
    t = mean_sr / se_sr if se_sr > 0 else float('nan')

    print(f"  Mean SR = {mean_sr:+.2f}%, sd = {sd_sr:.2f}, SE = {se_sr:.2f}, t = {t:.2f}")
    print(f"  Mean Δ  = {statistics.mean(deltas):+.4f}")

    n_positive = sum(1 for sr in srs if sr > 0)
    print(f"  Positive direction: {n_positive}/{n_seeds}")

    # 분기 평가
    print(f"\n=== Confirm 결과 ===")
    if mean_sr > 12 and t > 2:
        print(f"  ★★★ STRONGLY CONFIRMED — SR ~{mean_sr:.0f}% real, t={t:.1f}")
        print(f"  → arXiv 'polynomial SR' finding 의 핵심")
    elif mean_sr > 8:
        print(f"  ★★ CONFIRMED — substantial SR ({mean_sr:.1f}%)")
        print(f"  → polynomial scaling 강한 evidence")
    elif mean_sr > 3:
        print(f"  ★ partial — SR real but smaller than scan 추정")
    else:
        print(f"  ✗ NOT CONFIRMED — 17.86% scan 결과는 fluke 가능")


if __name__ == "__main__":
    main(sys.argv)
