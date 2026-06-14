"""
ym_multi_denom_fit.py — YM M_1,frac / Δ_ver,frac closed-form fit:
                       σ-scan emerges multi-denom at higher σ.

진단 발견 (2026-06-14):
  σ=0 sim 에서는 P_0(y) 가 peak 집중 → 어느 r_0 정의 (FIRST_VALID, BEST_CONV)
  를 써도 *single r_a* 만 나옴. Multi-denominator 분포는 *σ smearing* 의 결과.

본 script: screen 단계 건너뛰고 직접 σ-scan, multi-denom emergence + closed
form fit 검증.

설계:
  - 3 random Shor setups at N=437.
  - σ ∈ {0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500, 0.700, 1.000}.
  - 각 σ에서 BEST_CONV r_0 분포 측정 (3000 samples).
  - Multi-denom emergence: σ 따라 n_denom 증가 추적.
  - Closed form predicts m_σ(r_0) for ALL valid d ≤ N-1.
  - M_1,frac, Δ_ver,frac 측정 vs prediction (with r* shift handling).

실행:
  python -u -m experiments.ym_multi_denom_fit          # 표준 (~20 min)
  python -u -m experiments.ym_multi_denom_fit --fast   # 빠름 (~10 min)
"""
from __future__ import annotations
import math
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import convergent_denominators
from shor import simulate_period_finding
from noise import simulate_period_finding_noisy


N = 437
SIGMAS = [0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500, 0.700, 1.000]

FAST = "--fast" in sys.argv
N_SETUPS = 3
N_UNIFORM_SAMPLES = 5000 if not FAST else 2500
N_SIGMA_SAMPLES = 2000 if not FAST else 1000

RESULTS_FILE = Path("experiments/ym_multi_denom_fit_results.txt")


def find_best_convergent(k, a, N_, Q):
    """YM Eq.1: best convergent denominator (largest ≤ N-1); valid iff a^d ≡ 1."""
    cands = convergent_denominators(k, Q, N_ - 1)
    if not cands:
        return None
    d_best = cands[-1]
    if d_best > 0 and pow(a, d_best, N_) == 1:
        return d_best
    return None


def measure_r0_mass(a, N_, Q, sigma, n_samples, rng):
    """{r_0 → prob} mass under BEST_CONV definition."""
    mass = {}
    for _ in range(n_samples):
        if sigma > 0:
            m = simulate_period_finding_noisy(a, N_, rng=rng, phase_sigma=sigma)
        else:
            m = simulate_period_finding(a, N_, rng=rng)
        r_0 = find_best_convergent(m.k, a, N_, Q)
        if r_0 is not None:
            mass[r_0] = mass.get(r_0, 0.0) + 1.0 / n_samples
    return mass


def measure_S_frac(a, N_, Q, n_samples, rng):
    """{r_0 → |S_{r_0}|/Q} under BEST_CONV definition."""
    S = {}
    for _ in range(n_samples):
        k = int(rng.integers(0, Q))
        r_0 = find_best_convergent(k, a, N_, Q)
        if r_0 is not None:
            S[r_0] = S.get(r_0, 0.0) + 1.0 / n_samples
    return S


def predict_mass_per_r0(u, m_0, S_frac, all_r0):
    """Predict m_σ(r_0) for each r_0 in all_r0."""
    return {
        r_0: (1 - u) * S_frac.get(r_0, 0.0) + u * m_0.get(r_0, 0.0)
        for r_0 in all_r0
    }


def features_from_mass(mass):
    """Return (M_1_frac, Delta_ver_frac, sorted_pairs)."""
    if not mass:
        return 0.0, 0.0, []
    sorted_m = sorted(mass.items(), key=lambda kv: -kv[1])
    M_ver = sum(m for _, m in sorted_m)
    if M_ver == 0:
        return 0.0, 0.0, sorted_m
    M_1 = sorted_m[0][1]
    M_2 = sorted_m[1][1] if len(sorted_m) >= 2 else 0.0
    return M_1 / M_ver, (M_1 - M_2) / M_ver, sorted_m


def main():
    t_global = time.time()
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    log = []

    def out(*args):
        msg = " ".join(str(a) for a in args)
        print(msg, flush=True)
        log.append(msg + "\n")

    def hr(c="─", w=72):
        return c * w

    out(hr("="))
    out(f"║ YM Multi-Denom σ-scan (BEST_CONV definition, multi-denom emerges)")
    out(f"║ N={N}, Q={Q}, mode={'FAST' if FAST else 'STANDARD'}")
    out(f"║ Strategy: skip σ=0 screen, directly σ-scan to observe emergence")
    out(hr("="))
    out("")

    rng_setup = random.Random(202616)
    setups = []
    while len(setups) < N_SETUPS:
        a = rng_setup.randrange(2, N)
        if math.gcd(a, N) != 1 or any(s[0] == a for s in setups):
            continue
        r = classical_order(a, N)
        if r < 4:
            continue
        setups.append((a, r))

    # Per-setup aggregates
    all_m1_meas, all_m1_pred = [], []
    all_d_meas, all_d_pred = [], []

    for setup_idx, (a, r_a) in enumerate(setups, 1):
        out(hr("─"))
        out(f"### Setup {setup_idx}: a={a}, ord(a)={r_a}")
        out(hr("─"))

        # All valid d ≤ N-1: multiples of r_a within [1, N-1]
        all_valid_d = list(range(r_a, N, r_a))
        out(f"  All valid d ≤ N-1 (multiples of r_a): {all_valid_d[:8]}"
            f"{'...' if len(all_valid_d) > 8 else ''} (n={len(all_valid_d)})")

        # Noise-free histogram → m_0(r_0)
        out(f"  [Measuring noise-free m_0(r_0)... {N_SIGMA_SAMPLES} samples]")
        rng_0 = np.random.default_rng(setup_idx * 991)
        t_0 = time.time()
        m_0 = measure_r0_mass(a, N, Q, 0.0, N_SIGMA_SAMPLES, rng_0)
        out(f"  → m_0: {dict(sorted(m_0.items(), key=lambda kv: -kv[1])[:5])} "
            f"({time.time()-t_0:.0f}s)")

        # Uniform → |S_{r_0}|/Q
        out(f"  [Measuring |S_r0|/Q via {N_UNIFORM_SAMPLES} uniform samples...]")
        rng_u = np.random.default_rng(setup_idx * 991 + 47)
        t_u = time.time()
        S_frac = measure_S_frac(a, N, Q, N_UNIFORM_SAMPLES, rng_u)
        top_S = sorted(S_frac.items(), key=lambda kv: -kv[1])[:5]
        out(f"  → |S_ver|/Q={sum(S_frac.values()):.4f}, top: "
            + ", ".join(f"r={d}:{s*100:.2f}%" for d, s in top_S)
            + f" ({time.time()-t_u:.0f}s)")
        out("")

        # Header
        out(f"  {'σ':>6}  {'n_denom':>7}  {'M_1,frac':>9}  {'(pred)':>9}  {'diff':>7}  "
            f"{'r*meas':>7}  {'r*pred':>7}  {'Δ_ver':>9}  {'(pred)':>9}")
        out("  " + hr("─", 92))

        for sigma in SIGMAS:
            u = math.exp(-sigma * sigma)
            t_s = time.time()
            rng_sig = np.random.default_rng(
                setup_idx * 991 + int(sigma * 1e7) + 101
            )
            mass_s = measure_r0_mass(a, N, Q, sigma, N_SIGMA_SAMPLES, rng_sig)
            M_1_meas, D_meas, sorted_meas = features_from_mass(mass_s)
            r_star_meas = sorted_meas[0][0] if sorted_meas else 0

            # Closed-form prediction over ALL valid d
            mass_pred = predict_mass_per_r0(u, m_0, S_frac, all_valid_d)
            M_1_pred, D_pred, sorted_pred = features_from_mass(mass_pred)
            r_star_pred = sorted_pred[0][0] if sorted_pred else 0

            n_denom = len(mass_s)
            dt = time.time() - t_s
            out(f"  {sigma:>6.3f}  {n_denom:>7}  {M_1_meas:>9.4f}  {M_1_pred:>9.4f}  "
                f"{M_1_meas-M_1_pred:>+7.4f}  {r_star_meas:>7}  {r_star_pred:>7}  "
                f"{D_meas:>9.4f}  {D_pred:>9.4f}  ({dt:.0f}s)")

            all_m1_meas.append(M_1_meas)
            all_m1_pred.append(M_1_pred)
            all_d_meas.append(D_meas)
            all_d_pred.append(D_pred)
        out("")

    # Aggregate R²
    out(hr("="))
    out("║ Aggregate R² summary (BEST_CONV closed-form fit)")
    out(hr("="))
    out("")

    def r2(meas, pred):
        if len(meas) < 2:
            return 0.0, 0.0
        mm = statistics.mean(meas)
        ss_res = sum((m - p) ** 2 for p, m in zip(pred, meas))
        ss_tot = sum((m - mm) ** 2 for m in meas)
        return (1 - ss_res / ss_tot if ss_tot > 0 else 1.0,
                math.sqrt(ss_res / len(meas)))

    m1_r, m1_rmse = r2(all_m1_meas, all_m1_pred)
    d_r, d_rmse = r2(all_d_meas, all_d_pred)
    out(f"  Feature                  R²         RMSE     n")
    out("  " + hr("─", 50))
    out(f"  M_1,frac                {m1_r:>+7.4f}    {m1_rmse:.4f}   "
        f"{len(all_m1_meas)}")
    out(f"  Δ_ver,frac              {d_r:>+7.4f}    {d_rmse:.4f}   "
        f"{len(all_d_meas)}")
    out("")

    # Diagnose
    m1_meas_range = max(all_m1_meas) - min(all_m1_meas) if all_m1_meas else 0
    if m1_meas_range < 0.05:
        out("  Note: M_1,frac varies < 5% across σ — single-r_0 dominates")
        out("        throughout. Closed form holds trivially.")
    elif m1_r > 0.9:
        out("  ✓ Closed-form rational fit verified (R² > 0.9). YM features")
        out("    analytically derive from our framework.")
    elif m1_r > 0.5:
        out("  ~ Moderate fit. r* shift or finite-sample noise contributes.")
    else:
        out("  ✗ Poor fit. Investigate r* shift or framework assumption.")

    elapsed = time.time() - t_global
    out("")
    out(f"  Total elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    out(f"  Results saved to: {RESULTS_FILE}")

    RESULTS_FILE.write_text("".join(log), encoding="utf-8")


if __name__ == "__main__":
    main()
