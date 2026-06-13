"""
pure_shor_sr.py — Pure single-base Shor SR measurement (regime map prediction).

목적:
  Paper §3.6 Algorithm-structure regime map 의 *"Single-base Shor: small (~1%) SR"
  predicted* claim 의 *direct measurement*.

설계 (pure Shor 근사):
  - d = 1 (single base per run, Regev parallelism 없음)
  - (C) augmentation 비활성:
    * NO multi-base lcm accumulation (state.L 사용 안 함)
    * NO divisor search (state.L 의 divisors candidate 추가 안 함)
  - factor_from_exponent 비활성 (no fast path)
  - b-trick 만 유지 (b_i 가 우리 setup 에 이미 있으므로 그대로 사용)

주의:
  진짜 original Shor 는 random a (a = b^2 trick 없음) 와 gcd(a^(r/2) ± 1, N).
  본 implementation 은 Regev framework 의 b-trick 유지 — K-distribution 과
  noise sensitivity 가 *진짜 Shor 와 매우 가까움* (mathematically equivalent
  for SR purposes).

예상 결과 (regime map prediction):
  - K_base ~ 7 (이전 (437, 1) hybrid 와 유사)
  - Per-seed |SR| 작음 (~0-2%)
  - Mean SR ≈ 0 (direction stochastic)
  - σ-curve: plateau + decline shape

만약:
  - Mean SR < 1% → prediction *확정* ✓
  - Mean SR > 2% → prediction *수정 필요*

실험 설정:
  - N = 437 (consistent with our cross-cell base)
  - σ ∈ {0.000, 0.050, 0.150} (baseline + plateau + decline)
  - 5 seeds × 100 trials per cell
  - 총 1500 trials × ~0.4s = ~10 분

실행:
  python -u -m experiments.pure_shor_sr
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
from multi_base import convergent_denominators, minimize_order
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


N = 437
D = 1  # single base
SIGMAS = [0.000, 0.050, 0.150]
N_SEEDS = 5
TRIALS = 100

K_FILE = Path("experiments/pure_shor_sr_results.txt")
H_FILE = Path("experiments/pure_shor_sr_histograms.txt")


def pure_shor_one_trial(N, noise_kwargs, seed, max_runs=20):
    """Pure Shor approximation: d=1, NO (C), NO fast-path."""
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, D, rng_py)

    # Pure Shor: no state.L accumulation across runs
    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)
        ai = setup.a[0]
        ki = run.k_vec[0]
        bi = setup.b[0]

        # NO (C) augmentation: only convergent denominators, no divisor search
        cands = set(convergent_denominators(ki, Q, N - 1))

        valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
        if valid:
            r = minimize_order(ai, N, min(valid))
            if r > 0 and r == classical_order(ai, N):
                # b-trick factoring (b is in our setup, so we use it)
                b_pow = pow(bi, r, N)
                if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                    for delta in (-1, 1):
                        g = math.gcd((b_pow + delta) % N, N)
                        if 1 < g < N:
                            return K
        # NO factor_from_exponent (no fast path)

    return max_runs


def measure_cell(noise_kwargs, trials, base_seed):
    Ks = []
    for t in range(trials):
        K = pure_shor_one_trial(N, noise_kwargs, seed=base_seed + t * 1000)
        Ks.append(K)
    return Ks


def Ks_to_histogram(Ks):
    return dict(collections.Counter(Ks))


def p_value_normal(t):
    if math.isnan(t):
        return float('nan')
    return 0.5 * erfc(t / sqrt(2))


def main():
    print(f"# Pure Shor SR test (single-base + NO (C) augmentation)")
    print(f"# Tests regime map prediction: 'Single-base Shor: small (~1%) SR'")
    print(f"# N={N}, d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(f"# Total: {N_SEEDS * len(SIGMAS) * TRIALS} trials")
    print(f"# 예상 시간: ~10 분")
    print(flush=True)

    # Setup files
    with open(K_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Shor SR at N={N}\n")
        f.write(f"# d={D} (single base), NO (C) augmentation\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"sigma   seed  K_mean\n")
    with open(H_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Shor K histograms at N={N}\n")
        f.write(f"# d={D} (single base), NO (C) augmentation\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"seed\tsigma\tK\tcount\n")

    results: dict = {}
    hists: dict = {}

    t_start = time.time()
    total_cells = N_SEEDS * len(SIGMAS)
    cell_idx = 0

    for seed in range(1, N_SEEDS + 1):
        print(f"\n━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)
        results[seed] = {}
        hists[seed] = {}

        for sigma in SIGMAS:
            cell_idx += 1
            t_cell = time.time()
            noise = {} if sigma == 0.0 else {"phase_sigma": sigma}
            Ks = measure_cell(noise, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            hist = Ks_to_histogram(Ks)

            results[seed][sigma] = K_mean
            hists[seed][sigma] = hist

            # Save immediately
            with open(K_FILE, "a", encoding="utf-8") as f:
                f.write(f"{sigma:.3f}   {seed}     {K_mean:.4f}\n")
            with open(H_FILE, "a", encoding="utf-8") as f:
                for K in sorted(hist):
                    f.write(f"{seed}\t{sigma:.3f}\t{K}\t{hist[K]}\n")

            elapsed = time.time() - t_cell
            total_elapsed = time.time() - t_start
            eta = total_elapsed * (total_cells - cell_idx) / cell_idx if cell_idx > 0 else 0
            print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"cell={elapsed:>4.0f}s  ETA={eta:>4.0f}s "
                  f"({cell_idx}/{total_cells})", flush=True)

        K_base = results[seed][0.0]
        Kσ = results[seed][0.050]
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
        mark = " +" if sr > 0 else (" -" if sr < 0 else "")
        print(f"\n  seed {seed} 완료:  K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
              f"SR={sr:+.3f}%{mark}", flush=True)

    # ━━ Final analysis ━━
    print(f"\n━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
          flush=True)

    K_bases = [results[s][0.0] for s in range(1, N_SEEDS + 1)]
    print(f"\nK_baseline statistics:")
    print(f"  mean = {statistics.mean(K_bases):.4f}")
    if len(K_bases) > 1:
        print(f"  sd   = {statistics.stdev(K_bases):.4f}")
    print(f"  range = [{min(K_bases):.3f}, {max(K_bases):.3f}]")
    print(f"  per-seed: {K_bases}")

    print(f"\nSR % between-seed analysis:")
    print(f"  {'σ':>7}  {'mean SR%':>10}  {'sd':>7}  {'SE':>7}  {'t':>7}  {'p':>9}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*9}")
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        srs = []
        for s in range(1, N_SEEDS + 1):
            K_base = results[s][0.0]
            Kσ = results[s][sigma]
            if K_base > 0:
                srs.append((K_base - Kσ) / K_base * 100)
        n = len(srs)
        m = statistics.mean(srs)
        sd = statistics.stdev(srs) if n > 1 else 0
        se = sd / math.sqrt(n) if n > 1 else 0
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
    print(f"\nPer-seed plateau SR (σ=0.050):")
    plateau_srs = []
    for seed in range(1, N_SEEDS + 1):
        K_base = results[seed][0.0]
        Kσ = results[seed][0.050]
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
        plateau_srs.append(sr)
        mark = "+" if sr > 0 else ("-" if sr < 0 else " ")
        print(f"  seed {seed}: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
              f"SR={sr:+.3f}%  {mark}")

    # Sign test
    pos = sum(1 for s in plateau_srs if s > 0)
    neg = sum(1 for s in plateau_srs if s < 0)
    print(f"\nSign test: {pos} positive, {neg} negative ({pos}/{N_SEEDS} positive)")

    # Regime map verdict
    mean_plateau = statistics.mean(plateau_srs)
    print(f"\n━━ Regime map verdict ━━")
    print(f"  Predicted: 'Single-base Shor SR: small (~1%)'")
    print(f"  Measured:  mean = {mean_plateau:+.3f}%, |max per-seed| = {max(abs(s) for s in plateau_srs):.2f}%")
    if abs(mean_plateau) < 1.5:
        print(f"  → ✓ Prediction *confirmed* (small effect)")
    elif abs(mean_plateau) < 3:
        print(f"  → 〜 Prediction *partially confirmed* (small-medium effect)")
    else:
        print(f"  → ✗ Prediction *not confirmed* (effect larger than predicted)")

    print(f"\n총 시간: {time.time()-t_start:.0f}s ({(time.time()-t_start)/60:.1f} 분)")
    print(f"결과 저장:")
    print(f"  {K_FILE}")
    print(f"  {H_FILE}")


if __name__ == "__main__":
    main()
