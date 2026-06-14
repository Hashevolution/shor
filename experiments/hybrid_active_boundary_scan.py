"""
hybrid_active_boundary_scan.py — paper §3.6 의 active boundary cells 탐색 + fit.

목적:
  v0.3.0 §3.6.bis 의 self-correction 의 가장 강한 직접 증거:
  paper §3.6 의 active boundary cell (K_base ~ 2, "boundary flip mechanism" 관찰
  된 영역) 의 *fixed setup* 들을 base_seed sweep 으로 찾고, 거기에 closed form
  `K(σ) = K_∞ + (K_0 - K_∞)·exp(-σ²)` 를 fit.

  이전 `hybrid_sigma_curve.py` 는 base_seed 1-5 (모두 K_base ≈ 1.05, ceiling cells)
  만 잡았음. 작은 dynamic range 때문에 per-seed R² 0.39-0.88 변동 큼.

  본 script:
  1. base_seed 1-100 sweep, fixed setup K_base (50 noise-free trials/seed)
  2. K_base ∈ [1.8, 2.5] 인 active boundary setups 식별
  3. 그 setups (최대 5개) 에 σ-curve 측정 (8 σ × 200 trials)
  4. closed form 3-parameter fit → per-seed R² + aggregate R²

기대:
  - Active boundary cells 에서 dynamic range 크므로 R² > 0.9 가능.
  - paper §3.6 의 boundary-flip 어휘 없이 동일 cell 의 σ-curve 가
    closed form 으로 정확히 설명됨을 입증.

Reproduction:
  python -u -m experiments.hybrid_active_boundary_scan
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
MAX_RUNS = 20

SCAN_SEEDS = 100
SCAN_TRIALS = 50
TARGET_KBASE_LOW = 1.8
TARGET_KBASE_HIGH = 2.5
MAX_SETUPS_TO_FIT = 5

FIT_SIGMAS = [0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500]
FIT_TRIALS = 200

RESULTS_FILE = Path("experiments/hybrid_active_boundary_results.txt")


def hybrid_one_trial(setup, sigma, rng_np, trial_seed):
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


def measure_K(setup, sigma, n_trials, rng_seed_base):
    Ks = []
    for t in range(n_trials):
        rng_np = np.random.default_rng(
            rng_seed_base * 991 + t * 17 + int(sigma * 1e6)
        )
        K = hybrid_one_trial(setup, sigma, rng_np,
                             trial_seed=rng_seed_base + t * 1000)
        Ks.append(K)
    return statistics.mean(Ks), statistics.stdev(Ks) if len(Ks) >= 2 else 0.0


def exp_fit(sigmas, Ks):
    """K(σ) = K_inf + (K_0 - K_inf) · exp(-σ²) linear LS fit."""
    us = [math.exp(-s * s) for s in sigmas]
    n = len(us)
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
    ss_res = sum((K_ - (a + b * u)) ** 2 for u, K_ in zip(us, Ks))
    ss_tot = sum((K_ - mean_K) ** 2 for K_ in Ks)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse = math.sqrt(ss_res / n)
    return K_inf, K_0, r2, rmse


def main():
    t0 = time.time()
    lines = []
    header = (
        f"# Hybrid active boundary cell scan + fit (paper §3.6 self-correction)\n"
        f"# N={N}, d={D}, σ scan: {FIT_SIGMAS}\n"
        f"# Sweep {SCAN_SEEDS} base_seeds × {SCAN_TRIALS} noise-free trials\n"
        f"# Target K_base ∈ [{TARGET_KBASE_LOW}, {TARGET_KBASE_HIGH}]\n"
        f"# Then σ-fit on up to {MAX_SETUPS_TO_FIT} active cells × {FIT_TRIALS} trials each\n\n"
    )
    print(header)
    lines.append(header)

    # ── Phase 1: K_base sweep ──
    print("## Phase 1: base_seed sweep for K_base estimation\n")
    lines.append("## Phase 1: base_seed sweep for K_base estimation\n")
    sweep_results = []
    t_sweep_start = time.time()
    for seed in range(1, SCAN_SEEDS + 1):
        rng_py = random.Random(seed)
        setup = regev_setup_bases(N, D, rng_py)
        Kmean, Ksd = measure_K(setup, 0.0, SCAN_TRIALS, seed)
        is_active = TARGET_KBASE_LOW <= Kmean <= TARGET_KBASE_HIGH
        marker = " ★" if is_active else ""
        line = (f"  seed {seed:>3}: K_base={Kmean:.3f}±{Ksd:.2f}"
                f" (a={setup.a}, b={setup.b}){marker}\n")
        print(line, end="")
        lines.append(line)
        sweep_results.append((seed, setup, Kmean, Ksd, is_active))
        # Progress
        if seed % 10 == 0:
            elapsed = time.time() - t_sweep_start
            eta = elapsed * (SCAN_SEEDS - seed) / seed
            print(f"  [progress] {seed}/{SCAN_SEEDS}, elapsed {elapsed:.0f}s, "
                  f"ETA {eta:.0f}s\n", end="")

    active_cells = [r for r in sweep_results if r[4]]
    summary = (
        f"\n# Sweep done in {time.time()-t_sweep_start:.0f}s. "
        f"Found {len(active_cells)} active boundary cells out of {SCAN_SEEDS} seeds.\n"
        f"# Active K_base distribution: "
        f"mean={statistics.mean([r[2] for r in active_cells]):.3f}, "
        f"min={min([r[2] for r in active_cells]):.3f}, "
        f"max={max([r[2] for r in active_cells]):.3f}\n\n"
    ) if active_cells else f"\n# No active cells found in {SCAN_SEEDS} seeds.\n\n"
    print(summary)
    lines.append(summary)

    if not active_cells:
        # Save sweep results and exit
        RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
        return

    # ── Phase 2: σ-curve fit on first MAX_SETUPS_TO_FIT active cells ──
    print(f"## Phase 2: σ-curve closed-form fit on top {MAX_SETUPS_TO_FIT} active cells\n")
    lines.append(f"## Phase 2: σ-curve closed-form fit on top {MAX_SETUPS_TO_FIT} active cells\n")

    # Pick cells closest to K_base=2.0
    active_cells_sorted = sorted(active_cells, key=lambda r: abs(r[2] - 2.0))
    fit_cells = active_cells_sorted[:MAX_SETUPS_TO_FIT]

    all_pred = []
    all_meas = []
    per_cell_fits = []

    for idx, (seed, setup, Kbase_est, _, _) in enumerate(fit_cells, 1):
        sec = f"\n### active cell {idx} (seed {seed}, K_base sweep est {Kbase_est:.3f})\n"
        print(sec, end="")
        lines.append(sec)
        sec2 = f"   a={setup.a}, b={setup.b}\n"
        print(sec2, end="")
        lines.append(sec2)

        K_meas_list = []
        for sigma in FIT_SIGMAS:
            t_c = time.time()
            Kmean, Ksd = measure_K(setup, sigma, FIT_TRIALS, seed)
            ct = time.time() - t_c
            row = f"  σ={sigma:.3f}: K_mean={Kmean:.3f}±{Ksd:.2f}  ({ct:.0f}s)\n"
            print(row, end="")
            lines.append(row)
            K_meas_list.append(Kmean)

        K_inf, K_0, r2, rmse = exp_fit(FIT_SIGMAS, K_meas_list)
        per_cell_fits.append((seed, Kbase_est, K_0, K_inf, r2, rmse))
        fit_line = (
            f"  fit: K_0={K_0:.3f}, K_∞={K_inf:.3f}, K_0-K_∞={K_0-K_inf:+.3f}, "
            f"R²={r2:+.4f}, RMSE={rmse:.3f}\n"
        )
        print(fit_line, end="")
        lines.append(fit_line)

        for sigma, K_meas in zip(FIT_SIGMAS, K_meas_list):
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
            f"\n# Aggregate fit (active boundary cells): R² = {r2:+.4f}, "
            f"RMSE = {rmse:.3f}, n = {len(all_meas)}\n"
        )
        print(agg, end="")
        lines.append(agg)

    # Per-cell summary
    print("\n## Per-cell summary\n")
    lines.append("\n## Per-cell summary\n")
    tbl = ("| cell | seed | K_base est | K_0 | K_∞ | K_0-K_∞ | R² | RMSE |\n"
           "|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    print(tbl, end="")
    lines.append(tbl)
    for i, (seed, kbase, K0, Kinf, r2, rmse) in enumerate(per_cell_fits, 1):
        row = (f"| {i} | {seed} | {kbase:.3f} | {K0:.3f} | {Kinf:.3f} | "
               f"{K0-Kinf:+.3f} | {r2:+.4f} | {rmse:.3f} |\n")
        print(row, end="")
        lines.append(row)

    elapsed = time.time() - t0
    footer = f"\n# Total elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
