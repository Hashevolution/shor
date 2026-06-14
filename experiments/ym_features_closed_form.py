"""
ym_features_closed_form.py — Yang-Markidis features의 closed form 검증.

목적:
  Yang-Markidis (arXiv:2605.16074, ICS Workshops '26) 의 4 features 가 우리
  universal form `P_σ(y) = (1-u)/Q + u·P_0(y)` (u = e^{-σ²}) 로부터 *analytical
  derive* 됨을 sim 검증.

Derivations:
  Define q_σ(y) := P_σ(y) - 1/Q. Then q_σ = u·q_0.

  (1) A_peak(σ): autocorrelation peak strength
        A(ℓ) = Σ_y q(y)·q(y+ℓ)
        A_σ(ℓ) = u² · A_0(ℓ)
        A_peak = max_{ℓ≠0} A(ℓ)/A(0) — *u²-scale cancels in ratio*.
        → **A_peak is σ-invariant** (under pure phase noise).

  (2) H_norm(σ): normalized entropy
        H_σ = -Σ P_σ(y) log P_σ(y)
        No clean closed form (non-linear in P). But monotone in σ:
          σ=0:  H_norm = H_0 / log Q
          σ→∞: H_norm → 1 (uniform).

  (3) M_{1,frac}(σ): dominant verified mass fraction
        m_σ(r_0) = (1-u)·|S_{r_0}|/Q + u·m_0(r_0)
        M_ver(σ) = (1-u)·|S_ver|/Q + u·M_ver(0)
        M_1(σ)   = (1-u)·|S_{r*}|/Q + u·m_0(r*)         (leading r* invariant)
        M_1,frac(σ) = M_1(σ)/M_ver(σ)
        → **rational function in u**.

  (4) Δ_ver,frac(σ): verified margin fraction
        M_2(σ) = (1-u)·|S_{r**}|/Q + u·m_0(r**)         (second leading r**)
        Δ_ver,frac(σ) = (M_1(σ) - M_2(σ)) / M_ver(σ)
        → also **rational function in u**.

검증 plan:
  - Pure Shor d=1 setup at N=437.
  - 3 setups × 8 σ levels (= 24 cells).
  - Per cell: 1500 noise-free samples to estimate P_σ histogram.
  - Compute 4 features from histogram + Yang-Markidis Table 2 definitions.
  - Compare against closed-form predictions.

Reproduction:
  python -u -m experiments.ym_features_closed_form
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
N_SAMPLES = 2000
N_SETUPS = 3

RESULTS_FILE = Path("experiments/ym_features_closed_form_results.txt")


def find_valid_setup(N_, rng):
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


def candidate_denominator_for_k(k, a, N_, Q):
    """Yang-Markidis's r_0(y) for k: leading convergent denominator that passes modular verification."""
    cands = list(convergent_denominators(k, Q, N_ - 1))
    # Pick smallest valid (matches paper's continued-fraction step picking k/Q best convergent)
    for d_ in cands:
        if d_ > 0 and pow(a, d_, N_) == 1:
            r = minimize_order(a, N_, d_)
            if r > 0:
                return r
    return None  # no verified candidate


def measure_histogram(a, N_, sigma, n_samples, rng):
    """Sample n_samples measurements with phase noise σ. Return frequency dict {k: count}."""
    hist = {}
    for _ in range(n_samples):
        if sigma > 0:
            m = simulate_period_finding_noisy(a, N_, rng=rng, phase_sigma=sigma)
        else:
            m = simulate_period_finding(a, N_, rng=rng)
        hist[m.k] = hist.get(m.k, 0) + 1
    return hist


def compute_p_dist(hist, n_samples):
    """Empirical P(y) from histogram, full Q-sized dict (zeros omitted)."""
    return {k: c / n_samples for k, c in hist.items()}


def compute_apeak(p_dist, Q):
    """A_peak = max_{ℓ≠0} A(ℓ)/A(0).  A(ℓ) = Σ q(y)·q(y+ℓ), q = p - 1/Q."""
    p_array = np.zeros(Q)
    for y, p in p_dist.items():
        p_array[y] = p
    u = 1.0 / Q
    q_array = p_array - u
    # Autocorrelation via FFT
    qf = np.fft.fft(q_array)
    acorr = np.real(np.fft.ifft(qf * np.conj(qf)))
    A0 = acorr[0]
    if A0 == 0:
        return 0.0
    A_peak = float(np.max(acorr[1:]) / A0)
    return A_peak


def compute_hnorm(p_dist, Q):
    """H_norm = -Σ p log p / log Q. Only nonzero p contribute."""
    h = 0.0
    for p in p_dist.values():
        if p > 0:
            h -= p * math.log(p)
    return h / math.log(Q)


def compute_verified_masses(p_dist, a, N_, Q):
    """For each measurement y, find its r_0(y) = candidate denominator. Aggregate masses m(r_0)."""
    m = {}  # r_0 -> mass
    for y, p in p_dist.items():
        r_0 = candidate_denominator_for_k(y, a, N_, Q)
        if r_0 is None:
            continue
        m[r_0] = m.get(r_0, 0.0) + p
    return m


def compute_m1_frac_and_margin(m_dict):
    """M_1,frac and Δ_ver,frac from verified-mass dict."""
    if not m_dict:
        return 0.0, 0.0
    sorted_m = sorted(m_dict.values(), reverse=True)
    M_ver = sum(sorted_m)
    if M_ver == 0:
        return 0.0, 0.0
    M_1 = sorted_m[0]
    M_2 = sorted_m[1] if len(sorted_m) > 1 else 0.0
    M_1_frac = M_1 / M_ver
    delta_ver_frac = (M_1 - M_2) / M_ver
    return M_1_frac, delta_ver_frac


def main():
    t0 = time.time()
    lines = []
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    header = (
        f"# Yang-Markidis features closed form 검증\n"
        f"# Under our universal form P_σ(y) = (1-u)/Q + u·P_0(y), u = exp(-σ²)\n"
        f"# Predictions:\n"
        f"#  A_peak(σ) = A_peak(0)            [σ-invariant — q_σ = u·q_0, ratio cancels]\n"
        f"#  H_norm(σ): monotone, no closed form (asymptotes 0→1)\n"
        f"#  M_1,frac(σ): rational function in u\n"
        f"#  Δ_ver,frac(σ): rational function in u\n"
        f"# N={N}, Q={Q}, σ ∈ {SIGMAS}\n"
        f"# {N_SETUPS} setups × {N_SAMPLES} samples per cell\n\n"
    )
    print(header)
    lines.append(header)

    # Aggregate per-feature R² across (setup, σ)
    apeak_meas = []
    apeak_pred = []
    m1_meas = []
    m1_pred = []
    margin_meas = []
    margin_pred = []
    hnorm_data = []  # only for reporting trend

    rng_setup = random.Random(202606)
    setups = [find_valid_setup(N, rng_setup) for _ in range(N_SETUPS)]

    for idx, (a, b, r) in enumerate(setups, 1):
        sec = f"## setup {idx}: a={a}, b={b}, r={r}\n"
        print(sec)
        lines.append(sec)

        # Noise-free histogram (σ=0) → derives all closed-form parameters
        rng_0 = np.random.default_rng(idx * 991)
        hist_0 = measure_histogram(a, N, 0.0, N_SAMPLES, rng_0)
        P_0 = compute_p_dist(hist_0, N_SAMPLES)
        A_peak_0 = compute_apeak(P_0, Q)
        H_norm_0 = compute_hnorm(P_0, Q)
        m_dict_0 = compute_verified_masses(P_0, a, N, Q)
        M_1_frac_0, delta_0 = compute_m1_frac_and_margin(m_dict_0)
        # Universal form predicts: m_∞(r_0) = (1-u)·|S_{r_0}|/Q at u=0.
        # |S_{r_0}| = Q · (uniform success density for that r_0). Estimate from uniform samples.
        # For simplicity here: take noise-free histogram support sizes (each y with positive p_0 sees its own r_0).
        # This is an approximation; for full rigor, separately sample uniform k and aggregate.
        baseline = (
            f"  baseline (σ=0): A_peak={A_peak_0:.4f}  H_norm={H_norm_0:.4f}  "
            f"M_1,frac={M_1_frac_0:.4f}  Δ_ver,frac={delta_0:.4f}\n"
        )
        print(baseline, end="")
        lines.append(baseline)

        # Estimate |S_{r_0}|/Q for each r_0 by sampling uniform k (5000 samples)
        rng_uniform = np.random.default_rng(idx * 991 + 1)
        uniform_counts = {}
        n_uniform = 5000
        for _ in range(n_uniform):
            k = int(rng_uniform.integers(0, Q))
            r_0 = candidate_denominator_for_k(k, a, N, Q)
            if r_0 is not None:
                uniform_counts[r_0] = uniform_counts.get(r_0, 0) + 1
        S_frac = {r_0: c / n_uniform for r_0, c in uniform_counts.items()}  # |S_{r_0}|/Q
        S_ver_frac = sum(S_frac.values())  # |S_ver|/Q

        # Identify r* (leading) and r** (second leading) from noise-free
        sorted_r0 = sorted(m_dict_0.items(), key=lambda kv: -kv[1])
        if not sorted_r0:
            continue
        r_star = sorted_r0[0][0]
        m0_rstar = sorted_r0[0][1]
        S_rstar = S_frac.get(r_star, 0.0)
        if len(sorted_r0) >= 2:
            r_second = sorted_r0[1][0]
            m0_rsecond = sorted_r0[1][1]
            S_rsecond = S_frac.get(r_second, 0.0)
        else:
            r_second = None
            m0_rsecond = 0.0
            S_rsecond = 0.0
        params_line = (
            f"  closed-form params: r*={r_star} (m_0={m0_rstar:.3f}, |S|/Q={S_rstar:.4f}), "
            f"r**={r_second} (m_0={m0_rsecond:.3f}, |S|/Q={S_rsecond:.4f}), "
            f"|S_ver|/Q={S_ver_frac:.4f}, M_ver(0)={sum(m_dict_0.values()):.4f}\n"
        )
        print(params_line, end="")
        lines.append(params_line)

        # Now scan σ values
        for sigma in SIGMAS:
            u = math.exp(-sigma * sigma)
            t_c = time.time()
            rng_s = np.random.default_rng(idx * 991 + int(sigma * 1e6) + 100)
            hist_s = measure_histogram(a, N, sigma, N_SAMPLES, rng_s)
            P_s = compute_p_dist(hist_s, N_SAMPLES)
            A_peak_s = compute_apeak(P_s, Q)
            H_norm_s = compute_hnorm(P_s, Q)
            m_dict_s = compute_verified_masses(P_s, a, N, Q)
            M_1_frac_s, delta_s = compute_m1_frac_and_margin(m_dict_s)

            # Closed-form predictions
            A_peak_pred = A_peak_0  # σ-invariant
            M_ver_pred = (1 - u) * S_ver_frac + u * sum(m_dict_0.values())
            M_1_pred = (1 - u) * S_rstar + u * m0_rstar
            M_2_pred = (1 - u) * S_rsecond + u * m0_rsecond
            M_1_frac_pred = M_1_pred / M_ver_pred if M_ver_pred > 0 else 0.0
            delta_pred = (M_1_pred - M_2_pred) / M_ver_pred if M_ver_pred > 0 else 0.0

            ct = time.time() - t_c
            row = (
                f"  σ={sigma:.3f}: A_peak={A_peak_s:.4f} (pred {A_peak_pred:.4f}), "
                f"H_norm={H_norm_s:.4f}, M_1,frac={M_1_frac_s:.4f} (pred {M_1_frac_pred:.4f}), "
                f"Δ={delta_s:.4f} (pred {delta_pred:.4f}) ({ct:.0f}s)\n"
            )
            print(row, end="")
            lines.append(row)

            apeak_meas.append(A_peak_s)
            apeak_pred.append(A_peak_pred)
            m1_meas.append(M_1_frac_s)
            m1_pred.append(M_1_frac_pred)
            margin_meas.append(delta_s)
            margin_pred.append(delta_pred)
            hnorm_data.append((sigma, H_norm_s))
        lines.append("\n")
        print()

    # Compute R² per feature
    def r2(meas, pred):
        if len(meas) < 2:
            return 0.0, 0.0
        mm = statistics.mean(meas)
        ss_res = sum((m - p) ** 2 for p, m in zip(pred, meas))
        ss_tot = sum((m - mm) ** 2 for m in meas)
        r2_val = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(meas))
        return r2_val, rmse

    summary = "\n## Closed-form fit summary\n\n"
    summary += "| Feature | R² | RMSE | n |\n|---|---:|---:|---:|\n"
    apeak_r2, apeak_rmse = r2(apeak_meas, apeak_pred)
    m1_r2, m1_rmse = r2(m1_meas, m1_pred)
    margin_r2, margin_rmse = r2(margin_meas, margin_pred)
    summary += f"| A_peak (σ-invariant pred) | {apeak_r2:+.4f} | {apeak_rmse:.4f} | {len(apeak_meas)} |\n"
    summary += f"| M_1,frac (rational in u)  | {m1_r2:+.4f} | {m1_rmse:.4f} | {len(m1_meas)} |\n"
    summary += f"| Δ_ver,frac (rational in u)| {margin_r2:+.4f} | {margin_rmse:.4f} | {len(margin_meas)} |\n"
    print(summary)
    lines.append(summary)

    # H_norm trend (sanity, no closed form)
    print("\n## H_norm trend (no closed form, expected monotone increase)\n")
    lines.append("\n## H_norm trend (no closed form, expected monotone increase)\n")
    for setup_idx in range(N_SETUPS):
        sigma_h = [hnorm_data[setup_idx * len(SIGMAS) + i] for i in range(len(SIGMAS))]
        line = f"  setup {setup_idx+1}: " + ", ".join(
            f"σ={s:.3f}→H={h:.4f}" for s, h in sigma_h
        ) + "\n"
        print(line, end="")
        lines.append(line)

    elapsed = time.time() - t0
    footer = f"\n# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
