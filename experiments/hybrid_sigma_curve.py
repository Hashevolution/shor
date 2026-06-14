"""
hybrid_sigma_curve.py — Hybrid (C)+b-trick closed-form fit (paper §3.6 setup).

목적:
  paper §3.6 의 active boundary cell (N=437, d=4) 의 Hybrid (C)+b-trick setup 에
  본 framework 의 closed form fit 적용.

  Pure single-base Shor 와 달리:
  - d = 4 multi-base (parallel measurements per run).
  - state.L 누적 (lcm of recovered orders): (C) augmentation 의 fast path.
  - factor_from_exponent (L 이 큰 약수 가지면 추가 회수 시도).

Closed form 가설:
  Hybrid 도 same form K(σ) ≈ K_∞ + (K_0 - K_∞) · exp(-σ²) 따를까?
  - K_0 = noise-free K_mean.
  - K_∞ = high-σ asymptote (≈ 1/p_uniform_d).
  - Decay rate = 1 (single QFT measurement per base).

  d-base joint per run: p_run(σ) = 1 - ∏_i (1 - p_i(σ)).
  Per-base p_i(σ) = ρ_i + (p_0_i - ρ_i)·exp(-σ²).
  → p_run 의 σ-curve 는 exp 의 polynomial → 여전히 monotone smooth.

설계:
  - Fixed setup (a, b)_i for i=1..4: same base_seed per cell.
  - 5 base_seeds × 8 σ × 200 trials.
  - 3-parameter fit K(σ²) = K_∞ + (K_0 - K_∞)·exp(-σ²).

Reproduction:
  python -u -m experiments.hybrid_sigma_curve
"""
from __future__ import annotations
import math
import random
import statistics
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    factor_from_exponent,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


N = 437
D = 4
SIGMAS = [0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500]
N_SEEDS = 5
TRIALS = 200
MAX_RUNS = 20

RESULTS_FILE = Path("experiments/hybrid_sigma_curve_results.txt")


def hybrid_one_trial_fixed(setup, sigma: float, rng_np: np.random.Generator,
                            trial_seed: int):
    """K-loop with FIXED setup (a, b) per trial. Only measurement noise varies."""
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    state = MultiBaseState()
    noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}

    for K in range(1, MAX_RUNS + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)
        for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
            cands = set(convergent_denominators(ki, Q, N - 1))
            if state.L > 1:
                cands.update(divisors(state.L))
            valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
            if valid:
                r = minimize_order(ai, N, min(valid))
                if r > 0 and r == classical_order(ai, N):
                    state.update(ai, r)
                    b_pow = pow(bi, r, N)
                    if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                        for delta in (-1, 1):
                            g = math.gcd((b_pow + delta) % N, N)
                            if 1 < g < N:
                                return K
        if state.L > 1:
            rng_f = random.Random(trial_seed)
            res = factor_from_exponent(N, state.L, rng_f, max_attempts=5)
            if res and 1 < res.factor < N:
                return K
    return MAX_RUNS


def measure_cell_fixed(setup, sigma: float, n_trials: int,
                        rng_seed_base: int) -> tuple[float, float]:
    Ks = []
    for t in range(n_trials):
        # Different RNG per trial (noise samples), but same setup.
        rng_np = np.random.default_rng(rng_seed_base * 991 + t * 17 + int(sigma * 1e6))
        K = hybrid_one_trial_fixed(setup, sigma, rng_np, trial_seed=rng_seed_base + t * 1000)
        Ks.append(K)
    return statistics.mean(Ks), statistics.stdev(Ks) if len(Ks) >= 2 else 0.0


def exp_fit(sigmas: list[float], Ks: list[float]) -> tuple[float, float, float, float]:
    """
    Fit K(σ) = K_inf + (K_0 - K_inf) · exp(-σ²).
    Solve via linear LS: let u = exp(-σ²), K = K_inf + (K_0 - K_inf)·u.
    Then K = K_inf·(1-u) + K_0·u.
    """
    us = [math.exp(-s * s) for s in sigmas]
    n = len(us)
    # Linear regression: K = a + b·u where a = K_inf, a+b = K_0.
    mean_u = sum(us) / n
    mean_K = sum(Ks) / n
    num = sum((u - mean_u) * (K_ - mean_K) for u, K_ in zip(us, Ks))
    den = sum((u - mean_u) ** 2 for u in us)
    if den <= 0:
        return mean_K, mean_K, 0.0, 0.0
    b = num / den
    a = mean_K - b * mean_u
    K_inf = a
    K_0 = a + b
    # R²
    ss_res = sum((K_ - (a + b * u)) ** 2 for u, K_ in zip(us, Ks))
    ss_tot = sum((K_ - mean_K) ** 2 for K_ in Ks)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse = math.sqrt(ss_res / n)
    return K_inf, K_0, r2, rmse


def main():
    t0 = time.time()
    lines = []
    header = (
        f"# Hybrid (C)+b-trick σ-curve closed-form fit (paper §3.6 setup)\n"
        f"# Setup: N={N}, d={D}, fixed (a, b) per cell.\n"
        f"# Closed form: K(σ) ≈ K_∞ + (K_0 - K_∞)·exp(-σ²)\n"
        f"# σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials × max_runs={MAX_RUNS}\n\n"
    )
    print(header)
    lines.append(header)

    all_pred = []
    all_meas = []
    per_seed_fits = []

    for seed in range(1, N_SEEDS + 1):
        rng_py = random.Random(seed)
        setup = regev_setup_bases(N, D, rng_py)
        sec = f"## seed {seed}: bases a={setup.a}, b={setup.b}\n"
        print(sec)
        lines.append(sec)

        K_measurements = []
        for sigma in SIGMAS:
            t_cell = time.time()
            K_mean, K_sd = measure_cell_fixed(setup, sigma, TRIALS, seed)
            ct = time.time() - t_cell
            row = f"  σ={sigma:.3f}: K_mean={K_mean:.3f}±{K_sd:.2f}  ({ct:.0f}s)\n"
            print(row, end="")
            lines.append(row)
            K_measurements.append(K_mean)

        # 3-parameter exp fit
        K_inf, K_0, r2, rmse = exp_fit(SIGMAS, K_measurements)
        per_seed_fits.append((K_0, K_inf, r2, rmse))
        fit_line = (
            f"  fit: K_0={K_0:.3f}, K_∞={K_inf:.3f}, "
            f"R²={r2:+.4f}, RMSE={rmse:.3f}\n\n"
        )
        print(fit_line, end="")
        lines.append(fit_line)

        # Add to aggregate
        for sigma, K_meas in zip(SIGMAS, K_measurements):
            K_pred = K_inf + (K_0 - K_inf) * math.exp(-sigma * sigma)
            all_pred.append(K_pred)
            all_meas.append(K_meas)

    # Aggregate R²
    if len(all_pred) >= 2:
        mm = statistics.mean(all_meas)
        ss_res = sum((m - p) ** 2 for m, p in zip(all_meas, all_pred))
        ss_tot = sum((m - mm) ** 2 for m in all_meas)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(all_meas))
        agg = (
            f"\n# Aggregate fit: R² = {r2:+.4f}, RMSE = {rmse:.3f}, "
            f"n = {len(all_meas)}\n"
        )
        print(agg, end="")
        lines.append(agg)

    # Summary
    print("\n## Per-seed exp fit summary\n")
    lines.append("\n## Per-seed exp fit summary\n")
    tbl = "| seed | K_0 | K_∞ | K_0 - K_∞ | R² | RMSE |\n|---:|---:|---:|---:|---:|---:|\n"
    print(tbl, end="")
    lines.append(tbl)
    for i, (K0, Kinf, r2, rmse) in enumerate(per_seed_fits, 1):
        row = f"| {i} | {K0:.3f} | {Kinf:.3f} | {K0 - Kinf:+.3f} | {r2:+.4f} | {rmse:.3f} |\n"
        print(row, end="")
        lines.append(row)

    elapsed = time.time() - t0
    footer = f"\n# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
