"""
grover_sr_focused.py — Grover SR focused scan (Benzi plateau 검증).

목적:
  grover_sr.py 의 1차 패스 결과 (mean SR null at σ=0.05, strongly positive at σ=0.15)
  를 정밀화하기 위한 focused scan.

  - ϕ_oracle 범위 축소 (U(-π/64, π/64)) → K_base 분포 tight.
  - k_iter = 2 (active boundary, K_base ~ 2-3 regime, §3.6 의 Goldilocks 와 매칭).
  - σ 범위 7-점 fine grid: {0.005, 0.025, 0.050, 0.075, 0.100, 0.150, 0.200}.
    Shor §3.6 의 plateau scan 범위와 동일.
  - 5 seeds × 200 trials = 1000 K 측정 per σ → 총 7000.

판정:
  - σ ∈ [0.005, 0.100] mean SR ≈ const (plateau) → §3.6 Benzi shape ✓
  - σ=0.150 mean SR shift (decline or rise) → 어떤 mechanism 인지 판정.
  - per-seed direction stochasticity 가 σ=0.005-0.100 plateau 에서 보존되는지.

Reproduction:
  python -u -m experiments.grover_sr_focused
"""
from __future__ import annotations
import math
import statistics
import time
from pathlib import Path

import numpy as np


N_SEARCH = 64
M_MARKED = 1
K_ITER = 2  # active boundary cell (K_base ~ 2.9 theoretical)
SIGMAS = [0.000, 0.005, 0.025, 0.050, 0.075, 0.100, 0.150, 0.200]
N_SEEDS = 5
TRIALS = 200
MAX_RUNS = 20

# Reduced oracle imperfection range: |k·ϕ| < 0.1 effective deviation per iteration
PHI_RANGE = math.pi / 64

RESULTS_FILE = Path("experiments/grover_sr_focused_results.txt")


def grover_p_marked(k_iter: int, sigma: float, phi_oracle: float, rng: np.random.Generator) -> float:
    theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))
    a_m = math.sin(theta)
    a_u = math.cos(theta)
    base_angle = 2 * theta + phi_oracle
    for _ in range(k_iter):
        angle = base_angle
        if sigma > 0:
            angle += rng.normal(0.0, sigma)
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        new_m = cos_a * a_m + sin_a * a_u
        new_u = -sin_a * a_m + cos_a * a_u
        a_m, a_u = new_m, new_u
    p = a_m * a_m
    return min(max(p, 0.0), 1.0)


def grover_one_trial(k_iter: int, sigma: float, phi_oracle: float, rng: np.random.Generator) -> int:
    for K in range(1, MAX_RUNS + 1):
        p = grover_p_marked(k_iter, sigma, phi_oracle, rng)
        if rng.random() < p:
            return K
    return MAX_RUNS


def seed_oracle_phase(seed: int) -> float:
    rng = np.random.default_rng(seed * 7919 + 12347)
    return float(rng.uniform(-PHI_RANGE, PHI_RANGE))


def measure_cell(k_iter: int, sigma: float, seed: int) -> tuple[float, float, list[int]]:
    phi_oracle = seed_oracle_phase(seed)
    rng = np.random.default_rng(seed * 65537 + int(sigma * 1e7) * 31 + k_iter)
    Ks = [grover_one_trial(k_iter, sigma, phi_oracle, rng) for _ in range(TRIALS)]
    mean_K = statistics.mean(Ks)
    sd_K = statistics.stdev(Ks) if len(Ks) >= 2 else 0.0
    return mean_K, sd_K, Ks


def main():
    t_start = time.time()
    lines = []
    header = (
        f"# Grover SR focused scan (Benzi plateau 검증)\n"
        f"# N_search={N_SEARCH} M={M_MARKED} k_iter={K_ITER}\n"
        f"# σ ∈ {SIGMAS}\n"
        f"# {N_SEEDS} seeds × {TRIALS} trials/cell × max_runs={MAX_RUNS}\n"
        f"# ϕ_oracle range: U(-π/64, π/64) ≈ U(-{PHI_RANGE:.4f}, +{PHI_RANGE:.4f})\n"
    )
    print(header)
    lines.append(header)

    theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))
    p_theory = math.sin((2 * K_ITER + 1) * theta) ** 2
    theory_line = f"# Theoretical p (ϕ=0): {p_theory:.4f}, K_theory ≈ {1/p_theory:.3f}\n\n"
    print(theory_line)
    lines.append(theory_line)

    # K_base[seed] from σ=0, then SR per (seed, σ)
    K_base: dict[int, float] = {}
    SR_table: dict[float, dict[int, float]] = {sigma: {} for sigma in SIGMAS if sigma > 0}
    K_distributions: dict[int, dict[float, list[int]]] = {seed: {} for seed in range(1, N_SEEDS + 1)}

    # σ=0 baseline first
    print("## Baseline σ=0.000\n")
    lines.append("## Baseline σ=0.000\n")
    for seed in range(1, N_SEEDS + 1):
        mean_K, sd_K, Ks = measure_cell(K_ITER, 0.0, seed)
        K_base[seed] = mean_K
        K_distributions[seed][0.0] = Ks
        line = f"seed {seed} (ϕ={seed_oracle_phase(seed):+.4f}): K_base={mean_K:.3f}±{sd_K:.2f}\n"
        print(line, end="")
        lines.append(line)
    lines.append("\n")
    print()

    # σ > 0 scan
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        section = f"## σ = {sigma}\n"
        print(section, end="")
        lines.append(section)
        for seed in range(1, N_SEEDS + 1):
            mean_K, sd_K, Ks = measure_cell(K_ITER, sigma, seed)
            sr = (K_base[seed] - mean_K) / K_base[seed] * 100 if K_base[seed] > 0 else 0.0
            SR_table[sigma][seed] = sr
            K_distributions[seed][sigma] = Ks
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
            f"t={t:+.2f}  direction={n_pos}+/{n_neg}-\n\n"
        )
        print(agg, end="")
        lines.append(agg)

    # Summary table (Benzi shape detection)
    print("## Summary (mean SR vs σ)\n")
    lines.append("## Summary (mean SR vs σ)\n")
    summary_hdr = "| σ | mean SR % | sd | SE | t | direction (n+/n-) |\n|---:|---:|---:|---:|---:|---:|\n"
    print(summary_hdr, end="")
    lines.append(summary_hdr)
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        srs = list(SR_table[sigma].values())
        mean_sr = statistics.mean(srs)
        sd_sr = statistics.stdev(srs) if len(srs) >= 2 else 0.0
        se = sd_sr / math.sqrt(len(srs)) if srs else 0.0
        t = mean_sr / se if se > 0 else 0.0
        n_pos = sum(1 for s in srs if s > 0)
        n_neg = sum(1 for s in srs if s < 0)
        row = f"| {sigma:.3f} | {mean_sr:+.3f} | {sd_sr:.3f} | {se:.3f} | {t:+.2f} | {n_pos}+/{n_neg}- |\n"
        print(row, end="")
        lines.append(row)

    # K_baseline mean
    K_base_mean = statistics.mean(K_base.values())
    K_base_sd = statistics.stdev(K_base.values())
    kbline = f"\n# K_baseline (mean ± sd over seeds): {K_base_mean:.3f} ± {K_base_sd:.3f}\n"
    print(kbline, end="")
    lines.append(kbline)

    elapsed = time.time() - t_start
    footer = f"# Elapsed: {elapsed:.1f}s\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
