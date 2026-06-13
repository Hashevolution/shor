"""
sigma_scan_N1147_d2_extended.py — (1147, 2) cell 정밀 재측정.

이전 (sigma_scan_general 1147 2 5 100 compact) 결과:
  seed 1: K_base=2.80, SR=+1.43% (K=2→K=1, classical)
  seed 2: K_base=3.25, SR=0.00%  (long jump)
  seed 3: K_base=3.39, SR=+9.44% (high-K rescue ★)
  seed 4: K_base=2.92, SR=+8.56% (high-K rescue ★)
  seed 5: K_base=2.24, SR=-2.68% (K=1→K=2, classical neg)
  Mean SR = +3.35%, p=0.082 (marginal)

본 실험 — 결정적 재확인:
  - 흥미로운 seeds 2, 3 + 새 seed 6 (high-K rescue 보편성 검증)
  - σ 값 7개 (보다 정밀한 σ-curve mapping)

설정:
  - N=1147, d=2
  - seeds [2, 3, 6]
  - σ ∈ {0.000, 0.025, 0.050, 0.075, 0.100, 0.150, 0.200} (7 값)
  - 100 trials per cell
  - 총 3 seeds × 7 σ × 100 = 2100 trials × ~3.4s = ~2 시간

목적:
  - seeds 2, 3 의 high-K rescue 재현 확인
  - seed 6 (새): high-K rescue 가 또 나타나는지 확인
  - 더 정밀한 σ-curve (plateau, decline 영역 정밀화)

실행:
  python -u -m experiments.sigma_scan_N1147_d2_extended

저장:
  experiments/sigma_scan_N1147_d2_extended_results.txt
  experiments/sigma_scan_N1147_d2_extended_histograms.txt
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


N = 1147
D = 2
SEEDS = [2, 3, 6]                                # 흥미 2,3 + 새 6
SIGMAS = [0.000, 0.025, 0.050, 0.075, 0.100, 0.150, 0.200]  # 7 값
TRIALS = 100

K_FILE = Path("experiments/sigma_scan_N1147_d2_extended_results.txt")
H_FILE = Path("experiments/sigma_scan_N1147_d2_extended_histograms.txt")


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


def init_files():
    if not K_FILE.exists():
        with open(K_FILE, "w", encoding="utf-8") as f:
            f.write(f"# extended sigma scan at N={N} d={D}\n")
            f.write(f"# seeds = {SEEDS}, sigmas = {SIGMAS}\n")
            f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"sigma   seed  K_mean\n")
    if not H_FILE.exists():
        with open(H_FILE, "w", encoding="utf-8") as f:
            f.write(f"# K histograms at N={N} d={D} (extended)\n")
            f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"seed\tsigma\tK\tcount\n")


def read_existing():
    K_means: dict = {}
    hists: dict = {}
    if K_FILE.exists():
        with open(K_FILE, encoding="utf-8", errors="replace") as f:
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
                    K_means.setdefault(seed, {})[sigma] = K_mean
                except ValueError:
                    continue
    if H_FILE.exists():
        with open(H_FILE, encoding="utf-8", errors="replace") as f:
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
    return K_means, hists


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


def report(K_means_all, hists, label=""):
    seeds = sorted(K_means_all)
    n = len(seeds)
    if n == 0:
        return
    print(f"\n  {'='*70}")
    print(f"  결합 분석 {label} (n_seeds = {n})")
    print(f"  {'='*70}")

    K_bases = [K_means_all[s].get(0.0, 0) for s in seeds if 0.0 in K_means_all[s]]
    if K_bases:
        print(f"  K_base: mean={statistics.mean(K_bases):.4f}  "
              f"range=[{min(K_bases):.3f}, {max(K_bases):.3f}]")

    # Per-sigma table
    print(f"\n  Between-seed:")
    print(f"  {'σ':>7}  {'mean SR%':>10}  {'sd':>7}  {'SE':>7}  {'t':>7}  {'p':>9}")
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        srs = []
        for s in seeds:
            if 0.0 in K_means_all[s] and sigma in K_means_all[s]:
                K_base = K_means_all[s][0.0]
                Kσ = K_means_all[s][sigma]
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

    # Per-seed plateau SR + flip
    plateau = 0.050
    print(f"\n  Per-seed (σ={plateau}):")
    for s in seeds:
        if 0.0 not in K_means_all[s] or plateau not in K_means_all[s]:
            continue
        K_base = K_means_all[s][0.0]
        Kσ = K_means_all[s][plateau]
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
        mark = "+" if sr > 0 else ("-" if sr < 0 else " ")
        flip_str = ""
        if (s in hists and 0.0 in hists[s] and plateau in hists[s]):
            flip = identify_dominant_flip(hists[s][0.0], hists[s][plateau])
            if flip:
                kf, kt, mag = flip
                flip_str = f"  flip K={kf}→K={kt} mag={mag}"
                if kf >= 5 or kt >= 5:
                    flip_str += " ★ high-K"
        print(f"    seed {s:>2}: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
              f"SR={sr:+.3f}%  {mark}{flip_str}")


def main():
    print(f"# Extended sigma scan at N={N} d={D}")
    print(f"# Seeds: {SEEDS}")
    print(f"# σ values (7): {SIGMAS}")
    print(f"# trials: {TRIALS}")
    print(f"# Total: {len(SEEDS)} × {TRIALS} × {len(SIGMAS)} = "
          f"{len(SEEDS) * TRIALS * len(SIGMAS)} trials")

    Q_size = 1 << (2 * (N - 1).bit_length())
    print(f"# Q = {Q_size:,}, 예상 per-trial ~3.4s, 총 ~2 시간\n", flush=True)

    init_files()
    K_means_all, hists = read_existing()
    print(f"# 기존 데이터 seeds: {sorted(K_means_all)}\n", flush=True)

    to_run = [s for s in SEEDS if s not in K_means_all]
    if not to_run:
        print(f"# 모든 seeds 완료. 결합 분석만 표시.")
        report(K_means_all, hists, "(완료)")
        return
    print(f"# 신규 seeds: {to_run}\n", flush=True)

    t_start = time.time()
    total_cells = len(to_run) * len(SIGMAS)
    cell_idx = 0

    for seed in to_run:
        print(f"\n━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)
        for sigma in SIGMAS:
            cell_idx += 1
            t_cell = time.time()
            noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}
            Ks = measure_cell(N, D, noise_kwargs, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            hist = Ks_to_histogram(Ks)

            with open(K_FILE, "a", encoding="utf-8") as f:
                f.write(f"{sigma:.3f}   {seed}     {K_mean:.4f}\n")
            with open(H_FILE, "a", encoding="utf-8") as f:
                for K in sorted(hist):
                    f.write(f"{seed}\t{sigma:.3f}\t{K}\t{hist[K]}\n")

            K_means_all.setdefault(seed, {})[sigma] = K_mean
            hists.setdefault(seed, {})[sigma] = hist

            elapsed = time.time() - t_cell
            total_elapsed = time.time() - t_start
            eta = total_elapsed * (total_cells - cell_idx) / cell_idx if cell_idx > 0 else 0
            print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"cell={elapsed:.0f}s  ETA={eta:.0f}s "
                  f"({cell_idx}/{total_cells})", flush=True)

        K_base = K_means_all[seed][0.0]
        Kσ = K_means_all[seed][0.050]
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
        mark = " +" if sr > 0 else (" -" if sr < 0 else "")
        print(f"\n  seed {seed} 완료: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
              f"SR={sr:+.3f}%{mark}", flush=True)
        report(K_means_all, hists, f"(seed {seed} 후)")

    print(f"\n━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"총 소요: {time.time() - t_start:.0f}s "
          f"({(time.time() - t_start) / 60:.1f} 분)")


if __name__ == "__main__":
    main()
