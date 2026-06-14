"""
multi_noise_closed_form.py — Universal closed form 검증 across 4 noise models.

목적:
  v0.3.0 closed form `p(σ) = ρ + (p_0 - ρ)·exp(-σ²)` 의 universal form 일반화:
    p(noise) = (1 - ε(noise)) · p_0 + ε(noise) · g_∞(noise)
  를 4 noise models 에서 검증.

각 noise model의 ε와 g_∞:
  - **Phase σ** (이미 v0.3.0):  ε = 1 - exp(-σ²),    g_∞ = ρ = uniform success
  - **Depolarizing p**:          ε = p,                g_∞ = ρ = uniform success
  - **Bias zero p**:             ε = p,                g_∞ = I(k=0, a, b)
  - **Amp damping γ**:           ε = 1 - exp(-γQ_eff), g_∞ ≈ ρ_damped (weighted)

paper v0.2.1 §3.3 Theorem 3 와의 정합:
  `g_M(η) = (1 - η)·g_0 + η·g_unif_M` 가 위 form 의 정확한 사례.
  본 v0.3.0 framework는 그 generalization.

설계:
  - Pure Shor d=1 setup 사용 (shor_sigma_curve_model.py 와 동일 5 setups).
  - 각 noise model 별 ε scan: depol/bias p ∈ {0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9},
    amp_damp γ ∈ {0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2}.
  - 200 trials per (setup, noise model, ε value).
  - p_0, ρ 측정은 phase noise 실험과 공유 (cache 또는 re-measure).
  - 각 noise model 별 closed form fit + per-model R².

Reproduction:
  python -u -m experiments.multi_noise_closed_form
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
MAX_RUNS = 20
N_P0_SAMPLES = 1500
N_RHO_SAMPLES = 1500
N_K_TRIALS = 200

NOISE_LEVELS = {
    "depolarizing": [0.000, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90],
    "bias_zero":    [0.000, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90],
    "amplitude_damp": [0.0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2],
}

RESULTS_FILE = Path("experiments/multi_noise_closed_form_results.txt")


def success_indicator(k, a, b, N_, Q):
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


def measure_p0(a, b, N_, n, rng):
    success = 0
    for _ in range(n):
        m = simulate_period_finding(a, N_, rng=rng)
        if success_indicator(m.k, a, b, N_, m.Q):
            success += 1
    return success / n


def measure_rho(a, b, N_, Q, n, rng):
    success = 0
    for _ in range(n):
        k = int(rng.integers(0, Q))
        if success_indicator(k, a, b, N_, Q):
            success += 1
    return success / n


def success_at_k_zero(a, b, N_, Q):
    """I(k=0, a, b) — for bias_zero noise's g_∞."""
    return 1.0 if success_indicator(0, a, b, N_, Q) else 0.0


def measure_K(a, b, N_, noise_kwargs, n_trials, rng):
    Ks = []
    for _ in range(n_trials):
        K = MAX_RUNS
        for k_run in range(1, MAX_RUNS + 1):
            if noise_kwargs:
                m = simulate_period_finding_noisy(a, N_, rng=rng, **noise_kwargs)
            else:
                m = simulate_period_finding(a, N_, rng=rng)
            if success_indicator(m.k, a, b, N_, m.Q):
                K = k_run
                break
        Ks.append(K)
    return statistics.mean(Ks), statistics.stdev(Ks) if len(Ks) >= 2 else 0.0


def predicted_K(p, M=MAX_RUNS):
    if p <= 0:
        return float(M)
    q = 1.0 - p
    return (1.0 - q ** M) / p


def closed_form_p(noise_model, level, p_0, rho, I_zero):
    """Universal closed form p(ε) = (1-ε)·p_0 + ε·g_∞."""
    if noise_model == "depolarizing":
        eps = level
        g_inf = rho
    elif noise_model == "bias_zero":
        eps = level
        g_inf = I_zero
    elif noise_model == "amplitude_damp":
        # amp damping: amp[x] *= exp(-γ·x). Effective coherence loss after Q steps.
        # For small γ·Q, ε ≈ 1 - exp(-γ·Q/2). Empirically fit.
        # We use ε = 1 - exp(-γ · Q) as approximation; g_inf = rho (peak destruction → uniform-ish).
        eps = 1.0 - math.exp(-level * 100.0)  # γ·100 heuristic
        g_inf = rho
    else:
        raise ValueError(noise_model)
    return (1 - eps) * p_0 + eps * g_inf


def main():
    t0 = time.time()
    lines = []
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    header = (
        f"# Multi-noise closed-form universality verification\n"
        f"# Model: p(noise) = (1-ε(noise))·p_0 + ε(noise)·g_∞\n"
        f"# Noise levels: {NOISE_LEVELS}\n"
        f"# N={N}, Q={Q}, MC samples: p_0={N_P0_SAMPLES}, ρ={N_RHO_SAMPLES}\n"
        f"# {N_K_TRIALS} trials per (setup, noise, ε)\n\n"
    )
    print(header)
    lines.append(header)

    # Pick 3 setups (smaller than phase noise's 5 to save time, same seed scheme)
    rng_setup = random.Random(202606)
    setups = []
    for s in range(3):
        a, b, r = find_valid_setup(N, rng_setup)
        setups.append((a, b, r))

    # Per-noise-model aggregate
    model_results = {}  # model -> list of (K_pred, K_meas)

    for idx, (a, b, r) in enumerate(setups, 1):
        sec = f"## setup {idx}: a={a}, b={b}, r={r}\n"
        print(sec)
        lines.append(sec)

        # Measure p_0, rho, I_zero
        rng_p0 = np.random.default_rng(idx * 991 + 17)
        rng_rho = np.random.default_rng(idx * 991 + 23)
        p0 = measure_p0(a, b, N, N_P0_SAMPLES, rng_p0)
        rho = measure_rho(a, b, N, Q, N_RHO_SAMPLES, rng_rho)
        I_zero = success_at_k_zero(a, b, N, Q)
        meas_line = (
            f"  p_0 = {p0:.4f}, ρ = {rho:.4f}, I(k=0) = {I_zero:.0f}\n"
        )
        print(meas_line, end="")
        lines.append(meas_line)

        K_pred_baseline = predicted_K(p0)

        for noise_model, levels in NOISE_LEVELS.items():
            section = f"\n  ── {noise_model} ──\n"
            print(section, end="")
            lines.append(section)
            if noise_model not in model_results:
                model_results[noise_model] = []
            for level in levels:
                noise_kwargs = {} if level == 0.0 else {noise_model: level}
                rng_k = np.random.default_rng(
                    idx * 65537 + hash(noise_model) % 10000 + int(level * 1e6)
                )
                t_c = time.time()
                K_meas, K_sd = measure_K(a, b, N, noise_kwargs, N_K_TRIALS, rng_k)
                ct = time.time() - t_c

                p_pred = closed_form_p(noise_model, level, p0, rho, I_zero)
                K_pred = predicted_K(p_pred)

                row = (
                    f"  level={level:.4f}: p_pred={p_pred:.4f}  "
                    f"K_pred={K_pred:.3f}  K_meas={K_meas:.3f}±{K_sd:.2f}  "
                    f"diff={K_meas - K_pred:+.3f}  ({ct:.0f}s)\n"
                )
                print(row, end="")
                lines.append(row)
                model_results[noise_model].append((K_pred, K_meas))
        lines.append("\n")

    # Aggregate R² per noise model
    print("\n## Per-noise-model fit summary\n")
    lines.append("\n## Per-noise-model fit summary\n")
    tbl = "| noise model | n | R² | RMSE |\n|---|---:|---:|---:|\n"
    print(tbl, end="")
    lines.append(tbl)
    for model, pairs in model_results.items():
        preds = [p for p, m in pairs]
        meass = [m for p, m in pairs]
        mm = statistics.mean(meass)
        ss_res = sum((m - p) ** 2 for p, m in zip(preds, meass))
        ss_tot = sum((m - mm) ** 2 for m in meass)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(meass))
        row = f"| {model} | {len(pairs)} | {r2:+.4f} | {rmse:.3f} |\n"
        print(row, end="")
        lines.append(row)

    elapsed = time.time() - t0
    footer = f"\n# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
