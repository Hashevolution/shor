"""
pure_regev_sr.py — Pure Regev-like SR measurement (regime map prediction).

목적:
  Paper §3.6 Algorithm-structure regime map 의 *"Multi-base Regev (LLL):
  negative SR" predicted* claim 의 *direct measurement*.

설계:
  Real LLL Regev 의 *full implementation* 은 매우 복잡 (lattice 구성, LLL 환원,
  short vector extraction). 본 script 는 *AND-logic proxy* 를 사용:

  Regev LLL 의 *fragility* 구조:
    - 모든 d 좌표가 *clean* (noise-free 수준) 이어야 lattice 환원 성공
    - 한 좌표라도 *corrupted* → lattice 오염 → factor 못 찾음

  AND-logic proxy:
    - d 개 좌표 모두 ord(a_i) 회수 성공해야 run 성공 ((C) 의 OR-logic 의 *반대*)
    - 한 좌표라도 실패 → run 실패
    - Noise 가 어떤 좌표든 perturb 하면 → 전체 실패
    → LLL fragility 의 *구조적 proxy*

  추가:
    - NO (C) augmentation
    - NO divisor search
    - b-trick 으로 최종 factoring

주의 (honest framing):
  본 proxy 는 *실제 LLL Regev 가 아님*. Lattice fragility 의 *구조적* representation.
  *완전한 LLL 구현* 은 future work (sympy 의 LLL + lattice 구성).
  하지만 *regime map prediction* (negative SR) 의 *core 메커니즘* 은 그대로 검증.

예상:
  - d=4 AND-logic 의 baseline 성공률 ~50-70% (모든 coord 동시 성공 어려움)
  - K_base ~2-5
  - Noise → 더 자주 실패 → K_mean 증가 → ★ negative SR 예상

실험:
  - N = 437, d = 4 (matches hybrid (437, 4))
  - σ ∈ {0.000, 0.050, 0.150}
  - 5 seeds × 100 trials
  - 총 1500 trials × ~0.5s = ~15분

실행:
  python -u -m experiments.pure_regev_sr
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
D = 4  # multi-base
SIGMAS = [0.000, 0.050, 0.150]
N_SEEDS = 5
TRIALS = 100

K_FILE = Path("experiments/pure_regev_sr_results.txt")
H_FILE = Path("experiments/pure_regev_sr_histograms.txt")


def pure_regev_one_trial(N, d, noise_kwargs, seed, max_runs=20):
    """Pure Regev-like (AND-logic proxy for LLL fragility).

    Run 성공 조건:
      - 모든 d 좌표가 ord(a_i) 회수 성공
      - 적어도 하나의 좌표가 b-trick 의 nontrivial sqrt 충족
    """
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)

        # AND-logic: all d coords must succeed in r recovery
        rs = []
        all_succeed = True
        for ai, ki in zip(setup.a, run.k_vec):
            # NO (C) augmentation (no divisor search)
            cands = set(convergent_denominators(ki, Q, N - 1))

            valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
            if not valid:
                all_succeed = False
                break
            r = minimize_order(ai, N, min(valid))
            if r == 0 or r != classical_order(ai, N):
                all_succeed = False
                break
            rs.append(r)

        if not all_succeed:
            continue  # Any failed coord → entire run failed

        # All coords succeeded. Try b-trick on any of them.
        for ai, bi, r in zip(setup.a, setup.b, rs):
            b_pow = pow(bi, r, N)
            if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                for delta in (-1, 1):
                    g = math.gcd((b_pow + delta) % N, N)
                    if 1 < g < N:
                        return K  # ✓ Factor found

    return max_runs


def measure_cell(noise_kwargs, trials, base_seed):
    Ks = []
    for t in range(trials):
        K = pure_regev_one_trial(N, D, noise_kwargs, seed=base_seed + t * 1000)
        Ks.append(K)
    return Ks


def Ks_to_histogram(Ks):
    return dict(collections.Counter(Ks))


def p_value_normal(t):
    if math.isnan(t):
        return float('nan')
    return 0.5 * erfc(t / sqrt(2))


def main():
    print(f"# Pure Regev-like SR test (AND-logic proxy for LLL fragility)")
    print(f"# Tests regime map prediction: 'Multi-base Regev (LLL): negative SR'")
    print(f"# N={N}, d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(f"# Total: {N_SEEDS * len(SIGMAS) * TRIALS} trials")
    print(f"# 예상 시간: ~15 분")
    print(flush=True)

    with open(K_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev-like SR at N={N}\n")
        f.write(f"# d={D}, AND-logic proxy for LLL fragility\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"sigma   seed  K_mean\n")
    with open(H_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev-like K histograms at N={N}\n")
        f.write(f"# d={D}, AND-logic proxy\n")
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

    # Final analysis
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

    pos = sum(1 for s in plateau_srs if s > 0)
    neg = sum(1 for s in plateau_srs if s < 0)
    print(f"\nSign test: {pos} positive, {neg} negative ({pos}/{N_SEEDS} positive)")

    mean_plateau = statistics.mean(plateau_srs)
    print(f"\n━━ Regime map verdict ━━")
    print(f"  Predicted: 'Multi-base Regev (LLL): negative SR (LLL fragile)'")
    print(f"  Measured:  mean = {mean_plateau:+.3f}%, |max per-seed| = {max(abs(s) for s in plateau_srs):.2f}%")
    if mean_plateau < -0.5:
        print(f"  → ✓ Prediction *confirmed* (negative SR, LLL fragility)")
    elif mean_plateau > 0.5:
        print(f"  → ✗ Prediction *not confirmed* (positive SR — Regev resilient?)")
    else:
        print(f"  → 〜 Prediction *partially* (neutral, fragility not clear)")

    print(f"\n총 시간: {time.time()-t_start:.0f}s ({(time.time()-t_start)/60:.1f} 분)")
    print(f"결과 저장:")
    print(f"  {K_FILE}")
    print(f"  {H_FILE}")


if __name__ == "__main__":
    main()
