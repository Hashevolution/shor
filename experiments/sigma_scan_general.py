"""
sigma_scan_general.py — Generalized σ scan for any (N, d).

목적:
  (437, 4) 외 다른 (N, d) cell 에서 boundary-flip mechanism universality 검증.
  K_mean + K-histogram 저장 + per-seed flip 분석.

설계:
  - N, d 가변
  - 파일명: sigma_scan_N{N}_d{d}_results.txt / histograms.txt
  - 매 cell 즉시 저장 (중단 보존)
  - 매 seed 후 결합 분석 + flip 표 (early stopping 판단)

실행:
  python -u -m experiments.sigma_scan_general 1147 2 3 100      # (N=1147, d=2, 3 seeds × 100 trials)
  python -u -m experiments.sigma_scan_general 1147 2 5 200      # full scan
  python -u -m experiments.sigma_scan_general 437 3 13 200      # K_base=2.71 boundary shift test
  python -u -m experiments.sigma_scan_general 4087 4 5 200      # 큰 N universality

시간 추정:
  per-trial time ∝ Q × K_base / d ≈ N² × K_base / d
  (437, 4): ~0.23s/trial
  (1147, 2): ~2.3s/trial
  (4087, 4): ~4s/trial
"""

from __future__ import annotations
import collections
import math
import random
import statistics
import sys
import time
from math import comb, erfc, sqrt
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    factor_from_exponent,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


FULL_SIGMA_GRID = [
    0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.035, 0.050, 0.075, 0.100,
    0.150, 0.200,
]

# Compact grid: 5 σ covering full mechanism structure
#   0.000: baseline
#   0.025: plateau early (V3 measurement point)
#   0.050: plateau center (mechanism anchor)
#   0.150: decline onset
#   0.200: overload (direction asymmetry check)
COMPACT_SIGMA_GRID = [0.000, 0.025, 0.050, 0.150, 0.200]

# Minimal grid: just plateau center + baseline
MINIMAL_SIGMA_GRID = [0.000, 0.050]

SIGMA_GRID = FULL_SIGMA_GRID  # default, overridden in main() if grid_mode arg given


def hybrid_one_trial(N, d, noise_kwargs, seed, max_runs=20):
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)
    state = MultiBaseState()
    for K in range(1, max_runs + 1):
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
            rng_f = random.Random(seed)
            res = factor_from_exponent(N, state.L, rng_f, max_attempts=5)
            if res and 1 < res.factor < N:
                return K
    return max_runs


def measure_cell(N, d, noise_kwargs, trials, base_seed):
    Ks = []
    for t in range(trials):
        K = hybrid_one_trial(N, d, noise_kwargs, seed=base_seed + t * 1000)
        Ks.append(K)
    return Ks


def Ks_to_histogram(Ks):
    return dict(collections.Counter(Ks))


def p_value_normal(t):
    if math.isnan(t):
        return float('nan')
    return 0.5 * erfc(t / sqrt(2))


def identify_dominant_flip(h0, h_sigma):
    K_range = sorted(set(h0) | set(h_sigma))
    diffs = {K: h_sigma.get(K, 0) - h0.get(K, 0) for K in K_range}
    gains = [(K, d) for K, d in diffs.items() if d > 0]
    losses = [(K, -d) for K, d in diffs.items() if d < 0]
    if not gains or not losses:
        return None
    K_to, gain_mag = max(gains, key=lambda x: x[1])
    K_from, loss_mag = max(losses, key=lambda x: x[1])
    return K_from, K_to, min(gain_mag, loss_mag)


def report_combined(K_means_all, label=""):
    all_seeds = sorted(K_means_all)
    n = len(all_seeds)
    if n == 0:
        return
    K_bases = [K_means_all[s].get(0.0, 0) for s in all_seeds if 0.0 in K_means_all[s]]
    print(f"\n  {'='*70}")
    print(f"  결합 분석 {label} (n_seeds = {n})")
    print(f"  {'='*70}")
    if K_bases:
        print(f"  K_base: mean={statistics.mean(K_bases):.4f}  "
              f"range=[{min(K_bases):.3f}, {max(K_bases):.3f}]")

    print(f"\n  Between-seed table:")
    print(f"  {'σ':>7}  {'mean SR%':>10}  {'sd':>7}  {'SE':>7}  {'t':>7}  {'p':>9}")
    for sigma in SIGMA_GRID:
        if sigma == 0.0:
            continue
        srs = []
        for seed in all_seeds:
            if 0.0 in K_means_all[seed] and sigma in K_means_all[seed]:
                K_base = K_means_all[seed][0.0]
                Kσ = K_means_all[seed][sigma]
                if K_base > 0:
                    srs.append((K_base - Kσ) / K_base * 100)
        if not srs:
            continue
        m = statistics.mean(srs)
        sd = statistics.stdev(srs) if len(srs) > 1 else 0
        se = sd / math.sqrt(len(srs)) if len(srs) > 1 else 0
        t = m / se if se > 0 else float('nan')
        p = p_value_normal(t)
        marker = ""
        if not math.isnan(t):
            if abs(t) > 2:
                marker = " ★★"
            elif abs(t) > 1.5:
                marker = " ★"
        print(f"  {sigma:>7.3f}  {m:>+9.3f}%  {sd:>7.3f}  {se:>7.3f}  "
              f"{t:>+7.2f}  {p:>9.4f}{marker}")

    # Per-seed plateau SR
    plateau = 0.050
    print(f"\n  Per-seed plateau SR (σ={plateau}):")
    for seed in all_seeds:
        if 0.0 in K_means_all[seed] and plateau in K_means_all[seed]:
            K_base = K_means_all[seed][0.0]
            Kσ = K_means_all[seed][plateau]
            sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
            mark = "+" if sr > 0 else ("-" if sr < 0 else " ")
            print(f"    seed {seed:>2}: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
                  f"SR={sr:+.3f}%  {mark}")


def report_flips(hists, K_means_all, plateau=0.050):
    print(f"\n  {'='*70}")
    print(f"  Per-seed flip analysis (σ={plateau})")
    print(f"  {'='*70}")
    seeds_with_both = sorted(
        s for s in hists
        if 0.0 in hists[s] and plateau in hists[s]
        and 0.0 in K_means_all.get(s, {}) and plateau in K_means_all.get(s, {})
    )
    pos_count = neg_count = 0
    boundaries_pos = collections.Counter()
    boundaries_neg = collections.Counter()
    for seed in seeds_with_both:
        h0 = hists[seed][0.0]
        hs = hists[seed][plateau]
        K_base = K_means_all[seed][0.0]
        Kσ = K_means_all[seed][plateau]
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
        flip = identify_dominant_flip(h0, hs)
        if flip is None:
            continue
        K_from, K_to, mag = flip
        direction = "+" if sr > 0 else "-"
        if sr > 0:
            pos_count += 1
            boundaries_pos[(K_from, K_to)] += 1
        else:
            neg_count += 1
            boundaries_neg[(K_from, K_to)] += 1
        print(f"    seed {seed:>2}: K_base={K_base:.3f}  SR={sr:+.2f}%  {direction}  "
              f"K={K_from}→K={K_to}  mag={mag}")
    print(f"\n  Direction: {pos_count} positive, {neg_count} negative")
    if boundaries_pos:
        print(f"  Positive flip 분포: {dict(boundaries_pos)}")
    if boundaries_neg:
        print(f"  Negative flip 분포: {dict(boundaries_neg)}")
    # K=1/K=2 universality
    k12 = sum(c for (kf, kt), c in (list(boundaries_pos.items()) + list(boundaries_neg.items()))
              if {kf, kt} == {1, 2})
    total_flips = pos_count + neg_count
    if total_flips:
        print(f"  K=1/K=2 boundary 비율: {k12}/{total_flips} = {k12/total_flips*100:.1f}%")


def main(argv):
    global SIGMA_GRID
    if len(argv) < 3:
        print(f"Usage: python -m experiments.sigma_scan_general N d [n_seeds] [trials] [grid_mode]")
        print(f"  grid_mode: full (12 σ, default), compact (5 σ), minimal (2 σ)")
        return
    N = int(argv[1])
    D = int(argv[2])
    n_seeds = int(argv[3]) if len(argv) > 3 else 3
    trials = int(argv[4]) if len(argv) > 4 else 100
    grid_mode = argv[5] if len(argv) > 5 else "full"

    if grid_mode == "compact":
        SIGMA_GRID = COMPACT_SIGMA_GRID
    elif grid_mode == "minimal":
        SIGMA_GRID = MINIMAL_SIGMA_GRID
    else:
        SIGMA_GRID = FULL_SIGMA_GRID

    # File suffix: compact mode 는 별도 파일로 (full 결과와 섞이지 않게)
    suffix = "_compact" if grid_mode == "compact" else ("_minimal" if grid_mode == "minimal" else "")
    K_file = Path(f"experiments/sigma_scan_N{N}_d{D}{suffix}_results.txt")
    H_file = Path(f"experiments/sigma_scan_N{N}_d{D}{suffix}_histograms.txt")

    # Initialize files (no overwrite)
    if not K_file.exists():
        with open(K_file, "w", encoding="utf-8") as f:
            f.write(f"# sigma scan at N={N} d={D}\n")
            f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"sigma   seed  K_mean\n")
    if not H_file.exists():
        with open(H_file, "w", encoding="utf-8") as f:
            f.write(f"# K histograms at N={N} d={D}\n")
            f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"seed\tsigma\tK\tcount\n")

    # Load existing (resume support)
    K_means_all = {}
    hists = {}
    with open(K_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("sigma"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                sigma = float(parts[0])
                seed = int(parts[1])
                K_mean = float(parts[2])
                K_means_all.setdefault(seed, {})[sigma] = K_mean
            except ValueError:
                continue
    with open(H_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("seed"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                seed = int(parts[0])
                sigma = float(parts[1])
                K = int(parts[2])
                count = int(parts[3])
                hists.setdefault(seed, {}).setdefault(sigma, {})[K] = count
            except ValueError:
                continue

    existing = sorted(K_means_all)
    requested = list(range(1, n_seeds + 1))
    to_run = [s for s in requested if s not in existing]

    Q_size = 1 << (2 * max(1, (N - 1).bit_length()))
    print(f"# sigma scan at N={N} d={D}")
    print(f"# Q = 2^{2 * (N-1).bit_length()} = {Q_size:,}")
    print(f"# 기존 seeds: {existing}")
    print(f"# 신규 seeds: {to_run}")
    print(f"# {trials} trials × {len(SIGMA_GRID)} σ values")
    print(f"# Files: {K_file.name}, {H_file.name}")

    # Time estimate (per-trial scaling): for N=437 d=4: 0.23s; ∝ Q
    base_time = 0.23 * (Q_size / (1 << 18)) * (4 / D) * 1.5  # rough estimate
    total_trials = len(to_run) * trials * len(SIGMA_GRID)
    est_sec = total_trials * base_time
    print(f"# 추정 시간: ~{est_sec / 3600:.1f} 시간 ({est_sec / 60:.0f} 분)")
    print(flush=True)

    if not to_run:
        print(f"# 모든 seeds 완료. 기존 결과만 분석.\n")
        report_combined(K_means_all, label="(완료)")
        report_flips(hists, K_means_all)
        return

    t_start = time.time()
    total_cells = len(to_run) * len(SIGMA_GRID)
    cell_idx = 0

    for seed in to_run:
        print(f"\n━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)
        for sigma in SIGMA_GRID:
            cell_idx += 1
            t_cell = time.time()
            noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}
            Ks = measure_cell(N, D, noise_kwargs, trials, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            hist = Ks_to_histogram(Ks)

            # Save immediately
            with open(K_file, "a", encoding="utf-8") as f:
                f.write(f"{sigma:.3f}   {seed}     {K_mean:.4f}\n")
            with open(H_file, "a", encoding="utf-8") as f:
                for K in sorted(hist):
                    f.write(f"{seed}\t{sigma:.3f}\t{K}\t{hist[K]}\n")

            K_means_all.setdefault(seed, {})[sigma] = K_mean
            hists.setdefault(seed, {})[sigma] = hist

            elapsed = time.time() - t_cell
            t_total = time.time() - t_start
            eta = t_total * (total_cells - cell_idx) / cell_idx if cell_idx > 0 else 0
            print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"cell={elapsed:>4.0f}s  ETA={eta:>5.0f}s "
                  f"({cell_idx:>3}/{total_cells})", flush=True)

        # Per-seed summary + cumulative analysis
        K_base = K_means_all[seed][0.0]
        plateau_K = K_means_all[seed][0.050]
        sr = (K_base - plateau_K) / K_base * 100 if K_base > 0 else 0
        mark = " +" if sr > 0 else (" -" if sr < 0 else "")
        print(f"\n  seed {seed} 완료:  K_base={K_base:.4f}  K(σ=.05)={plateau_K:.4f}  "
              f"SR={sr:+.3f}%{mark}", flush=True)
        report_combined(K_means_all, f"(seed {seed} 후)")
        report_flips(hists, K_means_all)
        print(flush=True)

    print(f"\n━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"총 소요: {time.time() - t_start:.0f}s ({(time.time() - t_start) / 60:.1f} 분)")


if __name__ == "__main__":
    main(sys.argv)
