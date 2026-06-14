"""
qpe_isolated_sigma.py — QPE isolated σ-curve closed-form (internal consistency check).

목적:
  Shor 의 quantum kernel = Quantum Phase Estimation (QPE).
  Shor (with b-trick) σ-curve closed-form 이 fit 됨 (R²=+0.95).
  QPE 를 *isolated* (no b-trick) 로 떼어내서 동일 closed form 이 적용되는지 검증.

Success criterion (QPE / original Shor 1994):
  매 measurement k 에 대해, k/Q 의 convergent 중 어떤 분모 d 가
    pow(a, d, N) == 1 AND d == classical_order(a, N)
  를 만족하면 success. (no b-trick gating.)
  → success = "this measurement recovered the period r_a".

기대 결과:
  Closed form: p(σ) = ρ + (p_0 - ρ)·exp(-σ²), 동일.
  K(σ) = (1-(1-p)^M)/p, M=20.
  Aggregate R² > 0.9.

  Shor (with b-trick) 와의 차이:
  - p_0 더 큼 (b-trick gating 제거 → K_base 작음).
  - ρ 도 더 큼 (분자 success indicator 가 더 자주 만족).
  - 절대값 다르지만 같은 functional form.

설계:
  - N = 437.
  - 5 (a, N) setups with diverse r_a ∈ {3, 11, 99, 198, ...}.
  - σ ∈ {0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500}.
  - 500 trials per (setup, σ).
  - p_0 measured via 2000 noise-free runs.
  - ρ measured via 2000 uniform k samples.

Reproduction:
  python -u -m experiments.qpe_isolated_sigma
"""
from __future__ import annotations
import math
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
MAX_RUNS = 20
N_P0_SAMPLES = 2000
N_RHO_SAMPLES = 2000
N_K_TRIALS = 500

RESULTS_FILE = Path("experiments/qpe_isolated_sigma_results.txt")


def find_setups(N_: int, target_count: int = 5) -> list[tuple[int, int]]:
    """Pick (a, r_a) pairs with diverse r."""
    setups = {}  # r -> a
    for a in range(2, N_):
        if math.gcd(a, N_) != 1:
            continue
        r = classical_order(a, N_)
        if r < 2:
            continue
        if r not in setups:
            setups[r] = a
        if len(setups) >= 20:
            break
    # pick diverse r values
    rs = sorted(setups.keys())
    # spread across small + mid + large
    chosen = []
    step = max(1, len(rs) // target_count)
    for i in range(0, len(rs), step):
        chosen.append((setups[rs[i]], rs[i]))
        if len(chosen) >= target_count:
            break
    return chosen[:target_count]


def qpe_success_indicator(k: int, a: int, N_: int, Q: int, r_a: int) -> bool:
    """QPE success: convergent of k/Q gives period r_a."""
    cands = set(convergent_denominators(k, Q, N_ - 1))
    valid = [d_ for d_ in cands if d_ > 0 and pow(a, d_, N_) == 1]
    if not valid:
        return False
    r = minimize_order(a, N_, min(valid))
    return r > 0 and r == r_a


def measure_p0(a: int, N_: int, r_a: int, n: int, rng: np.random.Generator) -> float:
    success = 0
    for _ in range(n):
        m = simulate_period_finding(a, N_, rng=rng)
        if qpe_success_indicator(m.k, a, N_, m.Q, r_a):
            success += 1
    return success / n


def measure_rho(a: int, N_: int, r_a: int, Q: int, n: int, rng: np.random.Generator) -> float:
    success = 0
    for _ in range(n):
        k = int(rng.integers(0, Q))
        if qpe_success_indicator(k, a, N_, Q, r_a):
            success += 1
    return success / n


def measure_K(a: int, N_: int, r_a: int, sigma: float, n_trials: int,
              rng: np.random.Generator) -> tuple[float, float]:
    Ks = []
    for _ in range(n_trials):
        K = MAX_RUNS
        for k_run in range(1, MAX_RUNS + 1):
            if sigma > 0:
                m = simulate_period_finding_noisy(a, N_, rng=rng, phase_sigma=sigma)
            else:
                m = simulate_period_finding(a, N_, rng=rng)
            if qpe_success_indicator(m.k, a, N_, m.Q, r_a):
                K = k_run
                break
        Ks.append(K)
    return statistics.mean(Ks), statistics.stdev(Ks) if len(Ks) >= 2 else 0.0


def predicted_K(p: float, M: int = MAX_RUNS) -> float:
    if p <= 0:
        return float(M)
    return (1.0 - (1.0 - p) ** M) / p


def main():
    t0 = time.time()
    lines = []
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    setups = find_setups(N, 5)

    header = (
        f"# QPE isolated σ-curve closed-form (internal consistency)\n"
        f"# Model: p(σ) = ρ + (p_0 - ρ) · exp(-σ²),  E[K] = (1-(1-p)^M)/p, M={MAX_RUNS}\n"
        f"# Success: convergent of k/Q yields r_a (no b-trick gating).\n"
        f"# N = {N}, Q = {Q}, σ ∈ {SIGMAS}\n"
        f"# Setups: {len(setups)} (a, r_a) pairs: {[(a, r) for a, r in setups]}\n"
        f"# MC: p_0 = {N_P0_SAMPLES}, ρ = {N_RHO_SAMPLES}, K trials = {N_K_TRIALS}\n\n"
    )
    print(header)
    lines.append(header)

    all_pred = []
    all_meas = []

    for idx, (a, r_a) in enumerate(setups, 1):
        sec = f"## setup {idx}: a={a}, r_a={r_a}\n"
        print(sec)
        lines.append(sec)

        rng_p0 = np.random.default_rng(idx * 991 + 17)
        rng_rho = np.random.default_rng(idx * 991 + 23)
        t_mc = time.time()
        p0 = measure_p0(a, N, r_a, N_P0_SAMPLES, rng_p0)
        rho = measure_rho(a, N, r_a, Q, N_RHO_SAMPLES, rng_rho)
        mc_t = time.time() - t_mc
        se_p0 = math.sqrt(max(p0 * (1 - p0), 1e-12) / N_P0_SAMPLES)
        se_rho = math.sqrt(max(rho * (1 - rho), 1e-12) / N_RHO_SAMPLES)
        ml = (
            f"  p_0 = {p0:.4f} ± {se_p0:.4f}  ρ = {rho:.4f} ± {se_rho:.4f}  "
            f"Δ = {p0 - rho:+.4f}  ({mc_t:.0f}s)\n"
        )
        print(ml, end="")
        lines.append(ml)

        K_pred_0 = predicted_K(p0)

        for sigma in SIGMAS:
            decay = math.exp(-sigma * sigma)
            p_sigma = rho + (p0 - rho) * decay
            K_pred = predicted_K(p_sigma)

            rng_k = np.random.default_rng(idx * 65537 + int(sigma * 1e6))
            t_k = time.time()
            K_meas, K_sd = measure_K(a, N, r_a, sigma, N_K_TRIALS, rng_k)
            kt = time.time() - t_k

            row = (
                f"  σ={sigma:.3f}: p_σ={p_sigma:.4f}  "
                f"K_pred={K_pred:.3f}  K_meas={K_meas:.3f}±{K_sd:.2f}  "
                f"diff={K_meas - K_pred:+.3f}  ({kt:.0f}s)\n"
            )
            print(row, end="")
            lines.append(row)
            all_pred.append(K_pred)
            all_meas.append(K_meas)
        lines.append("\n")
        print()

    if len(all_pred) >= 2:
        mm = statistics.mean(all_meas)
        ss_res = sum((m - p) ** 2 for m, p in zip(all_meas, all_pred))
        ss_tot = sum((m - mm) ** 2 for m in all_meas)
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
