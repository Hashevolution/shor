"""
pure_regev_sr.py — Pure Regev with *actual LLL* implementation.

목적:
  Paper §3.6 Algorithm-structure regime map 의 *"Multi-base Regev (LLL):
  negative SR" predicted* claim 의 *faithful* measurement.

설계 (REAL LLL Regev):
  - d=4 (multi-base measurement)
  - Lattice construction from k_vec
  - REAL LLL reduction (Gram-Schmidt 기반 self-contained 구현)
  - Short vector → b-trick factoring

Lattice 구성 (Regev-style):
  Given measurement k_vec = (k_1, ..., k_d):
  B = [[S*Q,  0,  0, ..., 0],
       [S*k_1, 1, 0, ..., 0],
       [S*k_2, 0, 1, ..., 0],
       ...,
       [S*k_d, 0, 0, ..., 1]]

  Short vector z = (z_0, z_1, ..., z_d) satisfies:
    z_0 * Q + Σ z_i * k_i ≈ 0
  → Σ z_i * k_i / Q ≈ -z_0
  → 작은 z_i 들이 *integer relation* 표현
  → b = Π b_i^z_i mod N
  → if b² = 1 (mod N) AND b ∉ {1, N-1}: factor 추출

NOTE:
  - LLL implementation = Lenstra-Lenstra-Lovász 1982 algorithm
  - Gram-Schmidt orthogonalization 기반
  - Size reduction + Lovász condition
  - δ = 0.75 (standard)

예상:
  - Without noise: 단일 측정 도 factor 충분 가능 → K_base ~1-2
  - With noise: lattice basis 오염 → LLL output 일부 trial 에서 wrong vector
  - Negative SR 예상 (regime map prediction)

실험:
  - N = 437, d = 4
  - σ ∈ {0.000, 0.050, 0.150}
  - 5 seeds × 50 trials (LLL 느리므로 100 → 50)
  - 총 750 trials × ~1-2s = ~15-30 분

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
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


N = 437
D = 4
SIGMAS = [0.000, 0.050, 0.150]
N_SEEDS = 5
TRIALS = 50
LATTICE_SCALE = 1000  # S 값 (large for LLL precision)

K_FILE = Path("experiments/pure_regev_sr_results.txt")
H_FILE = Path("experiments/pure_regev_sr_histograms.txt")


# ────────────────────────────────────────────────────────────
# Self-contained LLL implementation (Lenstra-Lenstra-Lovász 1982)
# ────────────────────────────────────────────────────────────

def lll_reduce(B, delta=0.75, max_iters=1000):
    """LLL reduction. B = list of lists (basis vectors, integer entries).

    Returns reduced basis (modified in place).
    """
    n = len(B)
    B = [list(row) for row in B]

    def gram_schmidt():
        """Compute Gram-Schmidt orthogonalization and μ coefficients."""
        B_star = [list(row) for row in B]
        mu = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i):
                inner = sum(B[i][k] * B_star[j][k] for k in range(len(B[i])))
                norm_sq = sum(B_star[j][k] ** 2 for k in range(len(B_star[j])))
                if norm_sq > 0:
                    mu[i][j] = inner / norm_sq
                    for k in range(len(B_star[i])):
                        B_star[i][k] -= mu[i][j] * B_star[j][k]
        return B_star, mu

    B_star, mu = gram_schmidt()

    k = 1
    iters = 0
    while k < n and iters < max_iters:
        # Size reduction
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > 0.5:
                q = round(mu[k][j])
                if q != 0:
                    for i in range(len(B[k])):
                        B[k][i] -= q * B[j][i]
                    B_star, mu = gram_schmidt()

        # Lovász condition
        norm_k = sum(B_star[k][i] ** 2 for i in range(len(B_star[k])))
        norm_km1 = sum(B_star[k - 1][i] ** 2 for i in range(len(B_star[k - 1])))

        if norm_k >= (delta - mu[k][k - 1] ** 2) * norm_km1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            B_star, mu = gram_schmidt()
            k = max(k - 1, 1)
        iters += 1

    return B


# ────────────────────────────────────────────────────────────
# Regev lattice factoring
# ────────────────────────────────────────────────────────────

def b_pow_mod(base, exponent, modulus):
    """Compute base^exponent mod modulus for any integer exponent."""
    if exponent >= 0:
        return pow(base, exponent, modulus)
    # Negative exponent: use modular inverse
    inv = pow(base, -exponent, modulus)
    return pow(inv, -1, modulus)


def try_regev_lll_factor(N, Q, b_values, k_vec):
    """Build lattice from k_vec, run LLL, try b-trick.

    Returns True if factor found.
    """
    d = len(k_vec)
    S = LATTICE_SCALE

    # Build (d+1) x (d+1) lattice:
    #   B[0] = (S*Q, 0, 0, ..., 0)
    #   B[i+1] = (S*k_i, 0, ..., 1 at pos i+1, ..., 0)
    B = [[0] * (d + 1) for _ in range(d + 1)]
    B[0][0] = S * Q
    for i in range(d):
        B[i + 1][0] = S * int(k_vec[i])
        B[i + 1][i + 1] = 1

    # LLL reduce
    try:
        B_reduced = lll_reduce(B)
    except Exception:
        return False

    # Try each reduced vector as candidate
    for row in B_reduced:
        # Extract z = (z_1, ..., z_d) from positions 1..d
        z = [int(row[i + 1]) for i in range(d)]

        # Skip trivial vectors
        if all(zi == 0 for zi in z):
            continue
        if any(abs(zi) > 10000 for zi in z):  # too large
            continue

        # Compute b = Π b_i^z_i mod N
        try:
            b = 1
            for i in range(d):
                if z[i] != 0:
                    b = (b * b_pow_mod(b_values[i], z[i], N)) % N

            if b == 0 or b == 1 or b == N - 1:
                continue

            # Check b² ≡ 1 (mod N) → nontrivial sqrt
            if (b * b) % N == 1:
                for delta in (-1, 1):
                    g = math.gcd((b + delta) % N, N)
                    if 1 < g < N:
                        return True  # ✓ Factor found
        except (ValueError, OverflowError):
            continue

    return False


def pure_regev_one_trial(N, d, noise_kwargs, seed, max_runs=20):
    """Pure Regev with real LLL — K = runs to factor."""
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)

        # Try LLL on current measurement
        if try_regev_lll_factor(N, Q, setup.b, run.k_vec):
            return K

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
    print(f"# Pure Regev SR test (REAL LLL implementation)")
    print(f"# Tests regime map prediction: 'Multi-base Regev (LLL): negative SR'")
    print(f"# N={N}, d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(f"# Lattice: (d+1)x(d+1) with scale {LATTICE_SCALE}")
    print(f"# Total: {N_SEEDS * len(SIGMAS) * TRIALS} trials")
    print(f"# 예상 시간: ~15-30 분 (LLL Python 구현 느림)")
    print(flush=True)

    with open(K_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev (real LLL) SR at N={N}\n")
        f.write(f"# d={D}, lattice scale={LATTICE_SCALE}\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"sigma   seed  K_mean\n")
    with open(H_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev (real LLL) K histograms at N={N}\n")
        f.write(f"# d={D}, lattice scale={LATTICE_SCALE}\n")
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
        print(f"  → 〜 Prediction *partially* (neutral, fragility ambiguous)")

    print(f"\n총 시간: {time.time()-t_start:.0f}s ({(time.time()-t_start)/60:.1f} 분)")
    print(f"결과 저장:")
    print(f"  {K_FILE}")
    print(f"  {H_FILE}")


if __name__ == "__main__":
    main()
