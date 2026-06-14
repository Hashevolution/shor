"""
grover_sigma_curve_model.py — Grover σ-curve closed-form 검증.

목적:
  Grover circuit 의 σ-curve 가 closed-form 으로 정확히 예측되는지 verifies.

  Closed form derivation (per-iter iid Gaussian noise on rotation angle):
    각 iter: rotation angle = base + ε_i, ε_i ~ N(0, σ²) iid.
    초기상태: (sin θ, cos θ), θ = asin(√(M/N))
    base_angle = 2θ + ϕ_oracle (per-seed)
    After k iter: total angle accumulated T_k = k·base_angle + Σε_i, Σε_i ~ N(0, kσ²).
    Final marked amplitude = sin(θ + T_k) = sin(μ + Σε_i) where μ = (2k+1)θ + kϕ_oracle.
    p_marked = sin²(μ + Σε_i).
    E[p_marked] = (1 - E[cos(2(μ + Σε_i))])/2
                = (1 - cos(2μ)·E[cos(2·Σε_i)])/2
                = (1 - cos(2μ)·exp(-2kσ²))/2    [Σε ~ N(0, kσ²), char fn at 2]

  Predicted K_mean(σ) = 1 / p_mean(σ) (geometric draws).
  Predicted SR % = (K_0 - K_σ) / K_0 × 100 = 1 - p_0/p_σ at the mean level.

검증:
  - grover_sr_focused 의 5-seed × 7-σ 측정 결과 (K_ITER=2) 와 fit.
  - per-seed K_predicted vs K_measured 의 R².
  - cross-seed mean SR 의 predicted vs measured.

Reproduction:
  python -u -m experiments.grover_sigma_curve_model
"""
from __future__ import annotations
import math
import statistics
from pathlib import Path

import numpy as np

from experiments.grover_sr_focused import (
    N_SEARCH, M_MARKED, K_ITER, SIGMAS, N_SEEDS, TRIALS,
    seed_oracle_phase, measure_cell,
)


RESULTS_FILE = Path("experiments/grover_sigma_curve_model_results.txt")


def predicted_p(k_iter: int, sigma: float, phi: float) -> float:
    """Closed-form: p_mean(σ) = (1 - cos(2μ)·exp(-2kσ²))/2 where μ = (2k+1)θ + kϕ."""
    theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))
    mu = (2 * k_iter + 1) * theta + k_iter * phi
    decay = math.exp(-2 * k_iter * sigma * sigma)
    return (1.0 - math.cos(2 * mu) * decay) / 2.0


def predicted_K_mean(k_iter: int, sigma: float, phi: float, max_runs: int = 20) -> float:
    """
    Predicted mean K for geometric draws with success prob p, truncated at max_runs.
    E[K] for truncated geometric: (1 - (1-p)^max_runs) / p, but at max_runs cap.
    For p > 0.05 and max_runs=20, truncation negligible — use 1/p.
    """
    p = predicted_p(k_iter, sigma, phi)
    if p <= 0:
        return float(max_runs)
    # Exact truncated geometric expectation:
    # E[min(K, M)] = (1 - (1-p)^M) / p, where K is geometric.
    q = 1.0 - p
    qM = q ** max_runs
    return (1.0 - qM) / p


def main():
    lines = []
    header = (
        f"# Grover σ-curve closed-form 검증\n"
        f"# Model: E[p_marked](σ) = (1 - cos(2μ)·exp(-2kσ²))/2\n"
        f"#   μ = (2k+1)·θ + k·ϕ_oracle, θ = asin(√(M/N))\n"
        f"# N_search={N_SEARCH} M={M_MARKED} k_iter={K_ITER}\n"
        f"# σ ∈ {SIGMAS}\n"
        f"# Compare against measured (5 seeds × {TRIALS} trials × {len(SIGMAS)} σ)\n\n"
    )
    print(header)
    lines.append(header)

    theta = math.asin(math.sqrt(M_MARKED / N_SEARCH))

    # Per-seed analysis
    print("## Per-seed predictions\n")
    lines.append("## Per-seed predictions\n")
    K_measured_table = {}  # (seed, sigma) -> measured K
    K_predicted_table = {}  # (seed, sigma) -> predicted K

    for seed in range(1, N_SEEDS + 1):
        phi = seed_oracle_phase(seed)
        mu = (2 * K_ITER + 1) * theta + K_ITER * phi
        cos_2mu = math.cos(2 * mu)
        p_0 = predicted_p(K_ITER, 0.0, phi)

        seed_hdr = (
            f"### seed {seed} (ϕ={phi:+.4f}, μ={mu:.4f}, cos(2μ)={cos_2mu:+.4f}, "
            f"p_0={p_0:.4f}, K_0_pred={1/p_0:.3f})\n"
        )
        print(seed_hdr)
        lines.append(seed_hdr)

        # collect measured per σ (re-measure here to be self-contained)
        measured = {}
        for sigma in SIGMAS:
            mean_K, _, _ = measure_cell(K_ITER, sigma, seed)
            measured[sigma] = mean_K
            K_measured_table[(seed, sigma)] = mean_K

        K_0_meas = measured[0.0]
        for sigma in SIGMAS:
            K_pred = predicted_K_mean(K_ITER, sigma, phi)
            K_meas = measured[sigma]
            K_predicted_table[(seed, sigma)] = K_pred
            SR_pred = (K_0_meas - K_pred) / K_0_meas * 100 if K_0_meas > 0 else 0.0
            SR_meas = (K_0_meas - K_meas) / K_0_meas * 100 if K_0_meas > 0 else 0.0
            row = (
                f"  σ={sigma:.3f}: K_pred={K_pred:.3f} K_meas={K_meas:.3f} "
                f"diff={K_meas - K_pred:+.3f}  "
                f"SR_pred={SR_pred:+.3f}% SR_meas={SR_meas:+.3f}%\n"
            )
            print(row, end="")
            lines.append(row)
        lines.append("\n")

    # Cross-seed mean SR per σ
    print("\n## Cross-seed mean SR: predicted vs measured\n")
    lines.append("\n## Cross-seed mean SR: predicted vs measured\n")
    tbl_hdr = "| σ | mean SR predicted % | mean SR measured % | diff |\n|---:|---:|---:|---:|\n"
    print(tbl_hdr, end="")
    lines.append(tbl_hdr)
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        SR_preds = []
        SR_meass = []
        for seed in range(1, N_SEEDS + 1):
            K_0 = K_measured_table[(seed, 0.0)]
            K_pred = K_predicted_table[(seed, sigma)]
            K_meas = K_measured_table[(seed, sigma)]
            if K_0 > 0:
                SR_preds.append((K_0 - K_pred) / K_0 * 100)
                SR_meass.append((K_0 - K_meas) / K_0 * 100)
        SR_pred_mean = statistics.mean(SR_preds)
        SR_meas_mean = statistics.mean(SR_meass)
        row = f"| {sigma:.3f} | {SR_pred_mean:+.3f} | {SR_meas_mean:+.3f} | {SR_meas_mean - SR_pred_mean:+.3f} |\n"
        print(row, end="")
        lines.append(row)

    # Goodness-of-fit: per-seed R² of K_meas ~ K_pred
    print("\n## Per-seed goodness of fit (K_meas vs K_pred)\n")
    lines.append("\n## Per-seed goodness of fit (K_meas vs K_pred)\n")
    for seed in range(1, N_SEEDS + 1):
        Kms = [K_measured_table[(seed, s)] for s in SIGMAS]
        Kps = [K_predicted_table[(seed, s)] for s in SIGMAS]
        # R² coefficient
        K_mean = statistics.mean(Kms)
        ss_res = sum((m - p) ** 2 for m, p in zip(Kms, Kps))
        ss_tot = sum((m - K_mean) ** 2 for m in Kms)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(Kms))
        line = f"seed {seed}: R² = {r2:+.3f}, RMSE = {rmse:.3f} (over {len(SIGMAS)} σ values)\n"
        print(line, end="")
        lines.append(line)

    # Aggregate R² over all (seed, σ)
    all_Kms = []
    all_Kps = []
    for seed in range(1, N_SEEDS + 1):
        for sigma in SIGMAS:
            all_Kms.append(K_measured_table[(seed, sigma)])
            all_Kps.append(K_predicted_table[(seed, sigma)])
    agg_mean = statistics.mean(all_Kms)
    ss_res = sum((m - p) ** 2 for m, p in zip(all_Kms, all_Kps))
    ss_tot = sum((m - agg_mean) ** 2 for m in all_Kms)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse = math.sqrt(ss_res / len(all_Kms))
    agg = f"\n# Aggregate (all seeds × σ): R² = {r2:+.4f}, RMSE = {rmse:.3f}, n = {len(all_Kms)}\n"
    print(agg, end="")
    lines.append(agg)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
