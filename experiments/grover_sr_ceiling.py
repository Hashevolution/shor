"""
grover_sr_ceiling.py — Grover ceiling cell (K_base ≈ 1, no boundary).

목적:
  Grover SR regime map 의 ceiling cell (K_base ≈ 1) 측정.
  Paper §3.6 의 ceiling cell ((437, 8) K_base=1.19 → SR=+0.00%, (1147, 8) K_base=1.06
  → SR=+0.47%) 와 비교 가능한 Grover analog.

설계:
  - k_iter = 6 (theoretical p = 0.997, K_theory = 1.003)
  - 작은 ϕ_oracle range (U(-π/64, π/64))
  - σ ∈ {0.025, 0.050, 0.100} (focused scan 의 plateau 영역만)
  - 5 seeds × 200 trials

예측 (§3.6 regime map):
  - K_baseline ≈ 1.0-1.2
  - mean SR ≈ 0%, |per-seed SR| < 1-2% (no boundary trials to flip)

Reproduction:
  python -u -m experiments.grover_sr_ceiling
"""
from __future__ import annotations
import math
import statistics
import time
from pathlib import Path

import numpy as np

from experiments.grover_sr_focused import (
    N_SEARCH, M_MARKED, MAX_RUNS, PHI_RANGE,
    seed_oracle_phase, measure_cell,
)


K_ITER_CEIL = 6
SIGMAS_CEIL = [0.000, 0.025, 0.050, 0.100]
N_SEEDS = 5
TRIALS = 200

RESULTS_FILE = Path("experiments/grover_sr_ceiling_results.txt")


def main():
    t0 = time.time()
    lines = []
    header = (
        f"# Grover SR ceiling cell (K_base ≈ 1, §3.6 ceiling regime)\n"
        f"# N_search={N_SEARCH} M={M_MARKED} k_iter={K_ITER_CEIL}\n"
        f"# σ ∈ {SIGMAS_CEIL} (focused on plateau region only)\n"
        f"# {N_SEEDS} seeds × {TRIALS} trials/cell × max_runs={MAX_RUNS}\n"
        f"# ϕ_oracle range: U(-π/64, π/64) ≈ ±{PHI_RANGE:.4f}\n"
    )
    print(header)
    lines.append(header)

    theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))
    p_theory = math.sin((2 * K_ITER_CEIL + 1) * theta) ** 2
    th_line = f"# Theoretical p (ϕ=0): {p_theory:.5f}, K_theory ≈ {1/p_theory:.4f}\n\n"
    print(th_line)
    lines.append(th_line)

    K_base = {}
    SR_table = {s: {} for s in SIGMAS_CEIL if s > 0}

    print("## Baseline σ=0.000\n")
    lines.append("## Baseline σ=0.000\n")
    for seed in range(1, N_SEEDS + 1):
        mean_K, sd_K, _ = measure_cell(K_ITER_CEIL, 0.0, seed)
        K_base[seed] = mean_K
        line = f"seed {seed} (ϕ={seed_oracle_phase(seed):+.4f}): K_base={mean_K:.3f}±{sd_K:.2f}\n"
        print(line, end="")
        lines.append(line)
    lines.append("\n")

    for sigma in SIGMAS_CEIL:
        if sigma == 0.0:
            continue
        section = f"\n## σ = {sigma}\n"
        print(section, end="")
        lines.append(section)
        for seed in range(1, N_SEEDS + 1):
            mean_K, sd_K, _ = measure_cell(K_ITER_CEIL, sigma, seed)
            sr = (K_base[seed] - mean_K) / K_base[seed] * 100 if K_base[seed] > 0 else 0.0
            SR_table[sigma][seed] = sr
            line = f"seed {seed}: K={mean_K:.3f} SR={sr:+.3f}%\n"
            print(line, end="")
            lines.append(line)

        srs = list(SR_table[sigma].values())
        mean_sr = statistics.mean(srs)
        sd_sr = statistics.stdev(srs) if len(srs) >= 2 else 0.0
        se = sd_sr / math.sqrt(len(srs)) if srs else 0.0
        t = mean_sr / se if se > 0 else 0.0
        n_pos = sum(1 for s in srs if s > 0)
        n_neg = sum(1 for s in srs if s < 0)
        agg = (
            f"  → mean SR={mean_sr:+.3f}% sd={sd_sr:.3f} SE={se:.3f} "
            f"t={t:+.2f}  direction={n_pos}+/{n_neg}-\n"
        )
        print(agg, end="")
        lines.append(agg)

    K_base_mean = statistics.mean(K_base.values())
    K_base_sd = statistics.stdev(K_base.values())
    kbline = f"\n# K_baseline: {K_base_mean:.3f} ± {K_base_sd:.3f}\n"
    print(kbline, end="")
    lines.append(kbline)

    elapsed = time.time() - t0
    footer = f"# Elapsed: {elapsed:.1f}s\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
