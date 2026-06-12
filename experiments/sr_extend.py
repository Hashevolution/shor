"""
SR 확장: M3 (multi-noise combo) + N 확장.

M3: phase σ=0.05 + 작은 depol/amp 조합. 독립이면 효과 합산.
N 확장: N=1147, 2491 에서 SR 존재 여부. d 의 sweet spot 변화 확인.

실행:
    python -m experiments.sr_extend
"""

from __future__ import annotations
import math
import sys

from experiments.sr_amplify import measure_hybrid_extended


def m3_combo(N: int = 437, d: int = 4, trials: int = 500):
    """M3: phase σ=0.05 + tiny depol/amp 조합 SR."""
    print(f"\n# M3: 노이즈 조합 (N={N}, d={d}, {trials} trials)")
    print(f"  {'노이즈':<32}  {'hybrid K':>10}  {'vs baseline':>12}")

    baseline = measure_hybrid_extended(N, d, {}, trials=trials)
    print(f"  {'noise-free (baseline)':<32}  {baseline:>10.3f}  {'+0.000':>12}", flush=True)

    setups = [
        ("phase σ=0.05", {"phase_sigma": 0.05}),
        ("phase σ=0.05 + depol p=0.02", {"phase_sigma": 0.05, "depolarizing": 0.02}),
        ("phase σ=0.05 + depol p=0.05", {"phase_sigma": 0.05, "depolarizing": 0.05}),
        ("phase σ=0.05 + amp γ=0.00005", {"phase_sigma": 0.05, "amplitude_damp": 0.00005}),
        ("phase σ=0.05 + amp γ=0.0001",  {"phase_sigma": 0.05, "amplitude_damp": 0.0001}),
        ("phase σ=0.10 + depol p=0.02", {"phase_sigma": 0.10, "depolarizing": 0.02}),
        ("phase σ=0.10 + amp γ=0.0001", {"phase_sigma": 0.10, "amplitude_damp": 0.0001}),
    ]
    for label, kw in setups:
        K = measure_hybrid_extended(N, d, kw, trials=trials)
        diff = K - baseline
        marker = " ↓" if diff < -0.01 else (" ↑" if diff > 0.01 else "")
        print(f"  {label:<32}  {K:>10.3f}  {diff:>+12.3f}{marker}", flush=True)


def n_extend(trials: int = 200):
    """N 확장: SR 존재성 + sweet spot d 확인."""
    print(f"\n# N 확장: 큰 N 에서 SR 의 sweet spot (각 200 trials)")
    print(f"  {'N':>5} {'d':>2}  {'σ=0':>8}  {'σ=0.05':>8}  {'σ=0.10':>8}  "
          f"{'σ=0.20':>8}  {'SR %':>7}")
    for N in [1147, 2491]:
        for d in [2, 4, 8]:
            results = {}
            for sigma in [0.0, 0.05, 0.10, 0.20]:
                K = measure_hybrid_extended(
                    N, d,
                    {"phase_sigma": sigma} if sigma > 0 else {},
                    trials=trials,
                )
                results[sigma] = K
            baseline = results[0.0]
            min_K = min(results.values())
            sr = (baseline - min_K) / baseline * 100 if baseline > 0 else 0
            print(f"  {N:>5} {d:>2}  {results[0.0]:>8.3f}  "
                  f"{results[0.05]:>8.3f}  {results[0.10]:>8.3f}  "
                  f"{results[0.20]:>8.3f}  {sr:>6.2f}%", flush=True)


def main(argv):
    if "--m3" in argv:
        m3_combo()
    elif "--n" in argv:
        n_extend()
    else:
        m3_combo()
        n_extend()


if __name__ == "__main__":
    main(sys.argv)
