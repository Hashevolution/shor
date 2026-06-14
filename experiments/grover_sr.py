"""
grover_sr.py — Grover algorithm SR probe (paper §3.6 regime map 일반화 시험).

목적:
  paper §3.6 의 trial-level boundary-flip SR mechanism 이 Shor / Regev / Hybrid
  framework 의 *밖* (= 다른 양자 알고리즘) 으로 확장되는지 검증.

설계:
  - Grover circuit on N_search = 64, M = 1 marked element.
  - k_iter ∈ {1, 2, 3, 6} 으로 regime map 4 cells 동시 측정:
    · k=6: ceiling (K_base ~ 1.0)
    · k=3: active boundary (K_base ~ 1.7)
    · k=2: active boundary (K_base ~ 2.9)
    · k=1: noise floor (K_base ~ 7.4)
  - per seed: random imperfect oracle phase ϕ_seed ∈ U(-π/16, π/16) → seed
    별 K-distribution 변이.
  - Per iteration phase noise N(0, σ²) on rotation angle. K-loop 매 run 마다
    독립 sample.

Protocol (paper §3.6 와 정합):
  - σ ∈ {0.000, 0.050, 0.150}
  - 5 seeds × 100 trials × 3 σ × 4 cells = 6000 K 측정

Reproduction:
  python -u -m experiments.grover_sr

예측 결과 (paper §3.6 regime map):
  - k=6 (ceiling): mean |SR| < 0.5%
  - k=3 (active): per-seed |SR| 1-3%, direction stochastic
  - k=2 (active): per-seed |SR| 1-5%, direction stochastic
  - k=1 (noise floor): variance > effect, |mean SR| < SE

만약 prediction 일치 → §3.6 regime map 이 Shor 범위 *밖* 으로 확장됨 (= 일반화).
"""
from __future__ import annotations
import math
import statistics
import time
from pathlib import Path

import numpy as np


N_SEARCH = 64
M_MARKED = 1
K_ITER_CELLS = [1, 2, 3, 6]  # 4 regime cells
SIGMAS = [0.000, 0.050, 0.150]
N_SEEDS = 5
TRIALS = 100
MAX_RUNS = 20

RESULTS_FILE = Path("experiments/grover_sr_results.txt")


def grover_p_marked(k_iter: int, sigma: float, phi_oracle: float, rng: np.random.Generator) -> float:
    """
    Single Grover circuit execution 의 marked-subspace 확률 p.

    모델 (2D subspace, |m⟩ vs |u⟩):
      - 초기상태: (sin θ, cos θ), θ = asin(√(M/N))
      - 각 iteration: rotation by 2θ + ϕ_oracle + N(0, σ²)
      - ϕ_oracle: per-seed fixed (imperfect oracle phase)
      - sigma: per-iteration independent Gaussian noise

    Returns: p_marked ∈ [0, 1].
    """
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
    if p < 0.0:
        p = 0.0
    elif p > 1.0:
        p = 1.0
    return p


def grover_one_trial(k_iter: int, sigma: float, phi_oracle: float, rng: np.random.Generator) -> int:
    """K-loop: 매 run 마다 새 noise realization, K = first marked success."""
    for K in range(1, MAX_RUNS + 1):
        p = grover_p_marked(k_iter, sigma, phi_oracle, rng)
        if rng.random() < p:
            return K
    return MAX_RUNS


def seed_oracle_phase(seed: int) -> float:
    """Per-seed random imperfect oracle phase ∈ U(-π/16, π/16). 고정 seed → 고정 ϕ."""
    rng = np.random.default_rng(seed * 7919 + 12347)
    return float(rng.uniform(-math.pi / 16, math.pi / 16))


def measure_cell(k_iter: int, sigma: float, seed: int) -> tuple[float, float]:
    """Return (mean K, sd K) over TRIALS trials for given (k_iter, sigma, seed)."""
    phi_oracle = seed_oracle_phase(seed)
    rng = np.random.default_rng(seed * 65537 + int(sigma * 1e6) * 31 + k_iter)
    Ks = [grover_one_trial(k_iter, sigma, phi_oracle, rng) for _ in range(TRIALS)]
    mean_K = statistics.mean(Ks)
    sd_K = statistics.stdev(Ks) if len(Ks) >= 2 else 0.0
    return mean_K, sd_K


def main():
    t_start = time.time()
    lines = []
    header = (
        f"# Grover SR probe (paper §3.6 regime map 일반화 시험)\n"
        f"# N_search={N_SEARCH} M={M_MARKED} σ ∈ {SIGMAS} k_iter ∈ {K_ITER_CELLS}\n"
        f"# {N_SEEDS} seeds × {TRIALS} trials/cell × max_runs={MAX_RUNS}\n"
        f"# Per-seed imperfect oracle phase ϕ_seed ∈ U(-π/16, π/16)\n"
    )
    print(header)
    lines.append(header)

    # Per-cell aggregate
    for k_iter in K_ITER_CELLS:
        section_hdr = f"\n## k_iter = {k_iter} (theoretical p = sin²((2k+1)·θ) = "
        theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))
        p_theory = math.sin((2 * k_iter + 1) * theta) ** 2
        section_hdr += f"{p_theory:.3f}, K_theory ≈ {1/p_theory:.2f})\n"
        print(section_hdr)
        lines.append(section_hdr)

        # collect K_base[seed] and K_sigma[seed][sigma]
        K_base = {}
        K_sigma = {sigma: {} for sigma in SIGMAS if sigma > 0}
        SR_seed = {sigma: {} for sigma in SIGMAS if sigma > 0}

        for seed in range(1, N_SEEDS + 1):
            row = f"seed {seed} (ϕ={seed_oracle_phase(seed):+.4f}): "
            for sigma in SIGMAS:
                mean_K, sd_K = measure_cell(k_iter, sigma, seed)
                if sigma == 0.0:
                    K_base[seed] = mean_K
                    row += f"K(σ=0)={mean_K:.3f}±{sd_K:.2f}  "
                else:
                    K_sigma[sigma][seed] = mean_K
                    sr = (K_base[seed] - mean_K) / K_base[seed] * 100 if K_base[seed] > 0 else 0.0
                    SR_seed[sigma][seed] = sr
                    row += f"K(σ={sigma})={mean_K:.3f} SR={sr:+.2f}%  "
            print(row)
            lines.append(row + "\n")

        # cross-seed aggregate per sigma
        for sigma in SIGMAS:
            if sigma == 0.0:
                continue
            srs = list(SR_seed[sigma].values())
            mean_sr = statistics.mean(srs)
            sd_sr = statistics.stdev(srs) if len(srs) >= 2 else 0.0
            se = sd_sr / math.sqrt(len(srs)) if srs else 0.0
            t = mean_sr / se if se > 0 else 0.0
            n_pos = sum(1 for s in srs if s > 0)
            n_neg = sum(1 for s in srs if s < 0)
            agg = (
                f"  σ={sigma}: mean SR={mean_sr:+.3f}% sd={sd_sr:.3f} SE={se:.3f} "
                f"t={t:+.2f}  direction={n_pos}+/{n_neg}-\n"
            )
            print(agg, end="")
            lines.append(agg)

        # K_base summary
        K_base_mean = statistics.mean(K_base.values())
        K_base_sd = statistics.stdev(K_base.values()) if len(K_base) >= 2 else 0.0
        kbline = f"  K_baseline (mean over seeds): {K_base_mean:.3f}±{K_base_sd:.2f}\n"
        print(kbline, end="")
        lines.append(kbline)

    elapsed = time.time() - t_start
    footer = f"\n# Elapsed: {elapsed:.1f}s\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
