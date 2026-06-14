"""
shor_sigma_curve_model.py — Shor σ-curve closed-form 검증 (fixed (a, b) setup).

목적:
  Shor (pure single-base, b-trick) 의 phase noise σ-curve 가 closed-form 으로
  정확히 예측되는지 검증.

  Closed form (phase noise 의 noise-averaged 분포):
    각 amplitude x 에 e^{iε_x}, ε_x ~ N(0, σ²) iid:
      E_ε[|FFT(a·e^{iε})_k|²] = (1-e^{-σ²})/Q + e^{-σ²} · P_0(k)

    Per-run Bernoulli success prob:
      p(σ) = (1-e^{-σ²}) · ρ + e^{-σ²} · p_0 = ρ + (p_0 - ρ) · exp(-σ²)
    여기서 p_0 = Σ_k P_0(k)·I(k,a,b), ρ = |S_a|/Q (success set density).

    Truncated geometric K:
      E[K] = (1 - (1-p)^M) / p,  M = max_runs.

설계:
  - Fixed (a, b) per cell — NOT pure_shor_sr.py 의 (a, b) 변동 protocol.
  - 5 random (a, b) pair 를 sample, 각각이 nontrivial b-trick 조건 만족:
    ord(a) ≥ 2 AND b^{ord(a)} ∉ {1, N-1} AND (b^{ord(a)})² ≡ 1.
  - 각 setup: p_0, ρ Monte Carlo 측정.
  - 각 setup × σ ∈ {0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500}:
    K_meas = pure-Shor K-loop 500 trials.
  - K_pred(σ) = (1 - (1-p_σ)^M) / p_σ.

Reproduction:
  python -u -m experiments.shor_sigma_curve_model
"""
from __future__ import annotations
import math
import random
import statistics
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import convergent_denominators, minimize_order
from shor import simulate_period_finding
from noise import simulate_period_finding_noisy


N = 437
SIGMAS = [0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500]
N_SETUPS = 5
MAX_RUNS = 20

N_P0_SAMPLES = 2000
N_RHO_SAMPLES = 2000
N_K_TRIALS = 500

RESULTS_FILE = Path("experiments/shor_sigma_curve_model_results.txt")


def success_indicator(k: int, a: int, b: int, N_: int, Q: int) -> bool:
    cands = set(convergent_denominators(k, Q, N_ - 1))
    valid = [d_ for d_ in cands if d_ > 0 and pow(a, d_, N_) == 1]
    if not valid:
        return False
    r = minimize_order(a, N_, min(valid))
    if r <= 0 or r != classical_order(a, N_):
        return False
    b_pow = pow(b, r, N_)
    if b_pow in (1, N_ - 1):
        return False
    if (b_pow * b_pow) % N_ != 1:
        return False
    for delta in (-1, 1):
        g = math.gcd((b_pow + delta) % N_, N_)
        if 1 < g < N_:
            return True
    return False


def measure_p0(a: int, b: int, N_: int, n: int, rng: np.random.Generator) -> float:
    success = 0
    for _ in range(n):
        m = simulate_period_finding(a, N_, rng=rng)
        if success_indicator(m.k, a, b, N_, m.Q):
            success += 1
    return success / n


def measure_rho(a: int, b: int, N_: int, Q: int, n: int, rng: np.random.Generator) -> float:
    success = 0
    for _ in range(n):
        k = int(rng.integers(0, Q))
        if success_indicator(k, a, b, N_, Q):
            success += 1
    return success / n


def measure_K(a: int, b: int, N_: int, sigma: float, n_trials: int,
              rng: np.random.Generator) -> tuple[float, float]:
    """Pure Shor K-loop with fixed (a, b). Returns (mean K, sd K)."""
    Ks = []
    for _ in range(n_trials):
        K = MAX_RUNS
        for k_run in range(1, MAX_RUNS + 1):
            if sigma > 0:
                m = simulate_period_finding_noisy(a, N_, rng=rng, phase_sigma=sigma)
            else:
                m = simulate_period_finding(a, N_, rng=rng)
            if success_indicator(m.k, a, b, N_, m.Q):
                K = k_run
                break
        Ks.append(K)
    mean_K = statistics.mean(Ks)
    sd_K = statistics.stdev(Ks) if len(Ks) >= 2 else 0.0
    return mean_K, sd_K


def predicted_K(p: float, M: int = MAX_RUNS) -> float:
    if p <= 0:
        return float(M)
    q = 1.0 - p
    return (1.0 - q ** M) / p


def find_valid_setup(N_: int, rng: random.Random) -> tuple[int, int, int]:
    """Filter (a, b) pairs with nontrivial b-trick condition."""
    for _ in range(10000):
        b = rng.randrange(2, N_)
        if math.gcd(b, N_) != 1:
            continue
        a = (b * b) % N_
        if math.gcd(a, N_) != 1:
            continue
        r = classical_order(a, N_)
        if r < 2:
            continue
        b_pow = pow(b, r, N_)
        if b_pow in (1, N_ - 1):
            continue
        if (b_pow * b_pow) % N_ != 1:
            continue
        return a, b, r
    raise RuntimeError("Could not find valid setup")


def main():
    t0 = time.time()
    lines = []
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    header = (
        f"# Shor σ-curve closed-form 검증 (fixed (a, b) setup)\n"
        f"# Model: p(σ) = ρ + (p_0 - ρ) · exp(-σ²)\n"
        f"# E[K] = (1 - (1-p)^M) / p,  M = {MAX_RUNS}\n"
        f"# N = {N}, Q = {Q}, σ ∈ {SIGMAS}\n"
        f"# Setups: {N_SETUPS} (a, b) pairs with nontrivial b-trick\n"
        f"# MC: p_0 = {N_P0_SAMPLES} samples, ρ = {N_RHO_SAMPLES} samples\n"
        f"# K trials per (setup, σ): {N_K_TRIALS}\n\n"
    )
    print(header)
    lines.append(header)

    all_pred = []
    all_meas = []

    rng_setup = random.Random(202606)
    setups = []
    for s in range(N_SETUPS):
        a, b, r = find_valid_setup(N, rng_setup)
        setups.append((a, b, r))

    for idx, (a, b, r) in enumerate(setups, 1):
        seed_hdr = f"## setup {idx}: a={a}, b={b}, r={r}\n"
        print(seed_hdr)
        lines.append(seed_hdr)

        rng_p0 = np.random.default_rng(idx * 991 + 17)
        rng_rho = np.random.default_rng(idx * 991 + 23)
        t_mc = time.time()
        p0 = measure_p0(a, b, N, N_P0_SAMPLES, rng_p0)
        rho = measure_rho(a, b, N, Q, N_RHO_SAMPLES, rng_rho)
        mc_time = time.time() - t_mc
        se_p0 = math.sqrt(max(p0 * (1 - p0), 1e-12) / N_P0_SAMPLES)
        se_rho = math.sqrt(max(rho * (1 - rho), 1e-12) / N_RHO_SAMPLES)
        meas_line = (
            f"  p_0 = {p0:.4f} ± {se_p0:.4f}  ρ = {rho:.4f} ± {se_rho:.4f}  "
            f"diff = {rho - p0:+.4f}  ({mc_time:.0f}s)\n"
        )
        print(meas_line, end="")
        lines.append(meas_line)

        K_pred_0 = predicted_K(p0)

        for sigma in SIGMAS:
            decay = math.exp(-sigma * sigma)
            p_sigma = rho + (p0 - rho) * decay
            K_pred = predicted_K(p_sigma)

            rng_k = np.random.default_rng(idx * 65537 + int(sigma * 1e6))
            t_k = time.time()
            K_meas, K_sd = measure_K(a, b, N, sigma, N_K_TRIALS, rng_k)
            k_time = time.time() - t_k

            sr_pred = (K_pred_0 - K_pred) / K_pred_0 * 100 if K_pred_0 > 0 else 0.0
            K_meas_0 = None  # set below for σ=0
            if sigma == 0.0:
                K_meas_0 = K_meas
                sr_meas = 0.0
                pred_baseline_diff = K_meas - K_pred
            else:
                # Use K_meas at σ=0 from this setup (will be set above)
                pass

            row = (
                f"  σ={sigma:.3f}: p_σ={p_sigma:.4f}  "
                f"K_pred={K_pred:.3f}  K_meas={K_meas:.3f}±{K_sd:.2f}  "
                f"diff={K_meas - K_pred:+.3f}  ({k_time:.0f}s)\n"
            )
            print(row, end="")
            lines.append(row)
            all_pred.append(K_pred)
            all_meas.append(K_meas)
        lines.append("\n")
        print()

    if len(all_pred) >= 2:
        mean_meas = statistics.mean(all_meas)
        ss_res = sum((m - p) ** 2 for m, p in zip(all_meas, all_pred))
        ss_tot = sum((m - mean_meas) ** 2 for m in all_meas)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(all_meas))
        agg = (
            f"\n# Aggregate fit (all setups × σ): "
            f"R² = {r2:+.4f}, RMSE = {rmse:.3f}, n = {len(all_meas)}\n"
        )
        print(agg, end="")
        lines.append(agg)

    elapsed = time.time() - t0
    footer = f"# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
