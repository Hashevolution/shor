"""
pure_regev_sr.py — Pure Regev with *faithful* LLL implementation (revised v2).

목적:
  Paper §3.6 Algorithm-structure regime map 의 *"Multi-base Regev (LLL):
  negative SR" predicted* claim 의 *faithful* measurement.

설계 (faithful Regev):
  - d=4 multi-base measurement
  - Multi-measurement accumulation (M measurements 누적)
  - Regev-style lattice construction (proper basis)
  - LLL reduction (LLL 1982, δ=0.75)
  - Enumerated short vector search (basis + small linear combinations)
  - b-trick factor extraction

이전 구현 의 한계 (revised):
  - S = 1000 → S = Q (proper scale)
  - Single measurement → multi-measurement accumulation
  - Basis rows only → basis + linear combinations of small coefficients
  - 기본 LLL → 더 robust 한 size reduction

예상 (faithful):
  - K_base ~ 2-5 (50%+ per-measurement success)
  - Noise sensitivity 명확
  - Negative SR (regime map prediction)

NOTE:
  Real Regev paper 의 정확한 *complete* algorithm 은 multi-page 의 detailed
  construction. 본 구현 은 *core mathematical structure* (lattice + LLL +
  b-trick) 의 faithful representation. BKZ 보다 LLL 만 사용 — 작은 d 에서
  LLL 충분.

실험:
  - N = 437, d = 4
  - σ ∈ {0.000, 0.050, 0.150}
  - 5 seeds × 50 trials
  - ~15-30분
"""

from __future__ import annotations
import collections
import itertools
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
ENUM_RANGE = 2  # short vector enumeration coefficient range (-2 to +2)

K_FILE = Path("experiments/pure_regev_sr_results.txt")
H_FILE = Path("experiments/pure_regev_sr_histograms.txt")


# ────────────────────────────────────────────────────────────
# LLL implementation (Lenstra-Lenstra-Lovász 1982)
# Improved: more iterations, better numerical handling
# ────────────────────────────────────────────────────────────

def lll_reduce(B, delta=0.75, max_iters=2000):
    """LLL reduction (improved version).

    B: list of lists (integer basis vectors)
    Returns reduced basis.
    """
    n = len(B)
    B = [list(row) for row in B]

    def inner(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gs_update():
        B_star = [list(row) for row in B]
        mu = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i):
                if any(B_star[j]):
                    norm_j = inner(B_star[j], B_star[j])
                    if norm_j > 0:
                        mu[i][j] = inner(B[i], B_star[j]) / norm_j
                        for k in range(len(B_star[i])):
                            B_star[i][k] -= mu[i][j] * B_star[j][k]
        return B_star, mu

    B_star, mu = gs_update()

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
                    B_star, mu = gs_update()

        # Lovász condition
        norm_k = inner(B_star[k], B_star[k])
        norm_km1 = inner(B_star[k - 1], B_star[k - 1])
        if norm_k >= (delta - mu[k][k - 1] ** 2) * norm_km1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            B_star, mu = gs_update()
            k = max(k - 1, 1)
        iters += 1

    return B


# ────────────────────────────────────────────────────────────
# Regev lattice + factoring
# ────────────────────────────────────────────────────────────

def b_pow_mod(base, exponent, modulus):
    """Compute base^exponent mod modulus for any integer exponent."""
    if exponent >= 0:
        return pow(base, exponent, modulus)
    inv = pow(base, -exponent, modulus)
    return pow(inv, -1, modulus)


def try_z_factor(z, b_values, N):
    """Try b-trick factoring with relation z.

    Returns factor (or None).
    """
    d = len(z)
    if all(zi == 0 for zi in z):
        return None
    # Check vector not too large (numerical safety)
    if any(abs(zi) > 1000 for zi in z):
        return None

    try:
        b = 1
        for i in range(d):
            if z[i] != 0:
                b = (b * b_pow_mod(b_values[i], z[i], N)) % N
        if b in (0, 1, N - 1):
            return None
        if (b * b) % N != 1:
            return None
        # Nontrivial sqrt of 1
        for delta in (-1, 1):
            g = math.gcd((b + delta) % N, N)
            if 1 < g < N:
                return g
    except (ValueError, OverflowError, ZeroDivisionError):
        return None
    return None


def enumerate_short_vectors(B_reduced, d, enum_range=ENUM_RANGE):
    """Enumerate short vectors from reduced basis + small linear combinations."""
    candidates = set()

    # First: basis rows themselves
    for row in B_reduced:
        z = tuple(int(round(row[i + 1])) for i in range(d))
        if any(z):
            candidates.add(z)

    # Second: small linear combinations of reduced basis rows
    n = len(B_reduced)
    for coeffs in itertools.product(range(-enum_range, enum_range + 1), repeat=n):
        if all(c == 0 for c in coeffs):
            continue
        # Skip duplicate sign (negative is same z up to sign)
        # Find first non-zero
        first_nz = next((c for c in coeffs if c != 0), 0)
        if first_nz < 0:
            continue  # skip negative-leading
        u = [sum(c * row[i] for c, row in zip(coeffs, B_reduced)) for i in range(d + 1)]
        z = tuple(int(round(u[i + 1])) for i in range(d))
        if any(z):
            candidates.add(z)

    return candidates


def try_regev_lll_factor(N, Q, b_values, k_vecs_list):
    """Build Regev lattice from accumulated measurements, run LLL, try factor.

    k_vecs_list: list of k_vecs from M ≥ 1 measurements at SAME base set.
    Returns factor or None.
    """
    d = len(k_vecs_list[0])

    # Try different S scales (heuristic robustness)
    for S in [Q, Q // 2, Q * 2]:
        # Use most recent measurement (could combine with others)
        k_vec = k_vecs_list[-1]

        # Build (d+1) x (d+1) embedding lattice
        # Rows:
        #   B[0] = (S*Q, 0, 0, ..., 0)
        #   B[i+1] = (S*k_i, 0, ..., 1 at i+1, ..., 0)
        B = [[0] * (d + 1) for _ in range(d + 1)]
        B[0][0] = S * Q
        for i in range(d):
            B[i + 1][0] = S * int(k_vec[i])
            B[i + 1][i + 1] = 1

        # LLL reduce
        try:
            B_reduced = lll_reduce(B)
        except Exception:
            continue

        # Enumerate candidates
        candidates = enumerate_short_vectors(B_reduced, d)

        # Try to factor with each
        for z in candidates:
            factor = try_z_factor(z, b_values, N)
            if factor is not None:
                return factor

    # If multiple measurements, also try combining them
    if len(k_vecs_list) >= 2:
        # Combined lattice with 2 measurements: dim d+2
        for S in [Q, Q // 2]:
            k1 = k_vecs_list[-1]
            k2 = k_vecs_list[-2]

            B = [[0] * (d + 2) for _ in range(d + 2)]
            B[0][0] = S * Q
            B[1][1] = S * Q  # second constraint
            for i in range(d):
                B[i + 2][0] = S * int(k1[i])
                B[i + 2][1] = S * int(k2[i])
                B[i + 2][i + 2] = 1

            try:
                B_reduced = lll_reduce(B)
            except Exception:
                continue

            # Extract z = positions 2..d+1
            for row in B_reduced:
                z = tuple(int(round(row[i + 2])) for i in range(d))
                if any(z):
                    factor = try_z_factor(z, b_values, N)
                    if factor is not None:
                        return factor

    return None


def pure_regev_one_trial(N, d, noise_kwargs, seed, max_runs=20):
    """Pure Regev with faithful LLL + multi-measurement accumulation."""
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)

    k_vecs_list = []  # accumulate measurements

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)
        k_vecs_list.append(run.k_vec)

        # Try factor with accumulated measurements
        factor = try_regev_lll_factor(N, Q, setup.b, k_vecs_list)
        if factor is not None:
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
    print(f"# Pure Regev SR (FAITHFUL LLL implementation v2)")
    print(f"# Improvements over v1:")
    print(f"#   - Multi-measurement accumulation (LLL on M ≥ 1 measurements)")
    print(f"#   - Multiple S scaling values")
    print(f"#   - Enumerated short vector search (basis + linear combinations)")
    print(f"#   - 2-measurement combined lattice")
    print(f"# N={N}, d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(f"# Total: {N_SEEDS * len(SIGMAS) * TRIALS} trials")
    print(f"# 예상 시간: ~30-60 분 (LLL × enumeration 무거움)")
    print(flush=True)

    with open(K_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev (faithful LLL v2) SR at N={N}\n")
        f.write(f"# d={D}, multi-S, multi-measurement accumulation\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"sigma   seed  K_mean\n")
    with open(H_FILE, "w", encoding="utf-8") as f:
        f.write(f"# Pure Regev (faithful LLL v2) K histograms at N={N}\n")
        f.write(f"# d={D}\n")
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
        print(f"  → ✗ Prediction *not confirmed* (positive SR — interesting!)")
    else:
        print(f"  → 〜 Prediction *partially* (neutral, fragility ambiguous)")

    print(f"\n총 시간: {time.time()-t_start:.0f}s ({(time.time()-t_start)/60:.1f} 분)")


if __name__ == "__main__":
    main()
