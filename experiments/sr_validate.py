"""
V1: 다른 seed 로 SR 재현 검증.
V3: 2000 trials 로 sweet spot 의 statistical 강화.

목적: SR effect 가 seed=0 의 fluke 가 아닌지 확정.

실행:
    python -m experiments.sr_validate
"""

from __future__ import annotations
import math
import statistics
import sys

from experiments.sr_amplify import measure_hybrid_extended


def v1_seed_replication(N: int, d: int, n_seeds: int = 5, trials: int = 200):
    """V1: 5 개 독립 seed 에서 SR 측정. 일관되면 진짜."""
    print(f"\n# V1: N={N}, d={d}, {n_seeds} seeds × {trials} trials")
    print(f"  {'seed':>5}  {'σ=0':>10}  {'σ=0.05':>10}  {'SR %':>8}")

    srs = []
    for seed in range(1, n_seeds + 1):
        K0 = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
        K1 = measure_hybrid_extended(
            N, d, {"phase_sigma": 0.05}, trials=trials, seed=seed,
        )
        sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
        srs.append(sr)
        print(f"  {seed:>5}  {K0:>10.3f}  {K1:>10.3f}  {sr:>+7.2f}%", flush=True)

    mean_sr = statistics.mean(srs)
    if len(srs) > 1:
        sd_sr = statistics.stdev(srs)
        # 1-sample t-test against 0
        t = mean_sr / (sd_sr / math.sqrt(len(srs)))
    else:
        sd_sr = 0
        t = float('nan')

    print(f"\n  mean SR = {mean_sr:+.2f}%, sd = {sd_sr:.2f}, t = {t:.2f}")
    if mean_sr > 1.0 and t > 2:
        print(f"  ✓ V1 PASS — SR replicates across seeds (significant)")
    elif mean_sr > 0:
        print(f"  △ V1 weak — direction consistent but borderline")
    else:
        print(f"  ✗ V1 FAIL — no consistent SR")
    return srs


def v3_high_trials(N: int, d: int, trials: int = 2000, seed: int = 0):
    """V3: 2000 trials 로 sweet spot 의 SE 절반."""
    print(f"\n# V3: N={N}, d={d}, {trials} trials (high statistics)")
    print(f"  {'σ':>6}  {'hybrid K':>10}  {'SE est':>8}")

    K0 = measure_hybrid_extended(N, d, {}, trials=trials, seed=seed)
    sigmas = [0.025, 0.05, 0.075, 0.10, 0.20]
    print(f"  {0.0:>6.3f}  {K0:>10.4f}  {'-':>8}  (baseline)", flush=True)

    diffs = []
    for sigma in sigmas:
        K = measure_hybrid_extended(
            N, d, {"phase_sigma": sigma}, trials=trials, seed=seed,
        )
        diff = K - K0
        diffs.append(diff)
        # SE estimate (rough): assume std ~ 1.5, then SE ~ 1.5/sqrt(2000) ~ 0.034
        se_est = 1.5 / math.sqrt(trials)
        sig = "✓" if abs(diff) > 2 * se_est else " "
        print(f"  {sigma:>6.3f}  {K:>10.4f}  {se_est:>8.4f}  {sig} Δ = {diff:+.4f}",
              flush=True)

    mean_drop = -statistics.mean(diffs)
    print(f"\n  Mean reduction across σ ∈ {sigmas[0]}..{sigmas[-1]}: "
          f"{mean_drop:.4f} ({mean_drop/K0*100:.2f}%)")


def main(argv):
    print("=" * 60)
    print("SR 진위 확정 (V1 + V3)")
    print("=" * 60)

    # V1: seed 변경 — N=437 d=4 와 N=1147 d=2 둘 다
    v1_seed_replication(N=437, d=4, n_seeds=5, trials=300)
    v1_seed_replication(N=1147, d=2, n_seeds=5, trials=100)

    # V3: 2000 trials 강한 통계
    v3_high_trials(N=437, d=4, trials=2000)


if __name__ == "__main__":
    main(sys.argv)
