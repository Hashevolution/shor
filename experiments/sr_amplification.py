"""
sr_amplification.py — Engineered sub-functional baseline + noise SR amplification test.

목적:
  Mechanism amplification 가능성 *직접* 검증.
    - 정상 hybrid: K_base ≈ 1.92, SR = +0.144% (작음, our finding)
    - Thinned hybrid: 알고리즘을 *deliberately sub-functional* 로 만듦
      → K_base 크게 증가 (예: 5-10)
      → noise 가 *기능 복원* 가능성 → SR 증폭 가능

설계 (thinned hybrid):
  - convergent candidate 를 *smallest 1 개* 로 제한 (정상은 모든 convergent)
  - (C) divisor search 비활성 (정상은 state.L > 1 면 divisor pool 확장)
  - factor_from_exponent 비활성 (정상은 fast-path 추가)
  → 알고리즘 의 자원 제거 → "stuck" 영역 모방

실험:
  - N = 437, d = 4 (정상 (437, 4) 와 같은 cell)
  - σ ∈ {0.000, 0.050, 0.150} (3 값)
  - 3 seeds × 100 trials
  - 총 900 trials @ 0.25s = ~4분

비교:
  - 같은 seeds 에서 normal hybrid SR% (reference)
  - thinned hybrid SR%
  - 증폭 비율: thinned / normal

예상:
  - 만약 thinned baseline K_base 가 크고 (~5+)
  - 그 위에서 SR 가 ENAQT-style 큰 (예: 10-50%)
  - → amplification *확인* — paper §3.6 격상

실행:
  python -u -m experiments.sr_amplification
"""

from __future__ import annotations
import collections
import math
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    factor_from_exponent,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


N = 437
D = 4
SIGMAS = [0.000, 0.050, 0.150]
N_SEEDS = 3
TRIALS = 100


def hybrid_one_trial_thinned(N, d, noise_kwargs, seed, max_runs=20):
    """Thinned hybrid: smallest convergent only + (C) 비활성.

    Sub-functional baseline 만들기 위해:
      1. convergent candidates = {smallest 1개}
      2. divisor search 비활성 (state.L 누적은 유지하지만 candidate 추가 안 함)
      3. factor_from_exponent 비활성
    """
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)
    state = MultiBaseState()

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)
        for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
            # Thinned: only smallest convergent denominator
            all_cands = sorted(convergent_denominators(ki, Q, N - 1))
            cands = {all_cands[0]} if all_cands else set()

            # NO divisor search (state.L 의 divisors 추가 안 함)
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
        # NO factor_from_exponent (fast-path 비활성)

    return max_runs


def hybrid_one_trial_normal(N, d, noise_kwargs, seed, max_runs=20):
    """정상 hybrid (reference)."""
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


def measure(trial_fn, N, d, noise_kwargs, trials, base_seed):
    Ks = []
    for t in range(trials):
        K = trial_fn(N, d, noise_kwargs, seed=base_seed + t * 1000)
        Ks.append(K)
    return Ks


def main():
    print(f"# SR Amplification Test")
    print(f"# Thinned hybrid (smallest convergent only, no divisor search)")
    print(f"# vs Normal hybrid")
    print(f"# N={N} d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(flush=True)

    t_start = time.time()
    save_path = Path("experiments/sr_amplification_results.txt")
    with open(save_path, "w", encoding="utf-8") as f:
        f.write(f"# SR amplification at N={N} d={D}\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"mode    seed  sigma   K_mean\n")

    # Normal results
    normal: dict[int, dict[float, float]] = {}
    thinned: dict[int, dict[float, float]] = {}

    for seed in range(1, N_SEEDS + 1):
        print(f"━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

        # Normal
        normal[seed] = {}
        for sigma in SIGMAS:
            noise = {} if sigma == 0.0 else {"phase_sigma": sigma}
            t_cell = time.time()
            Ks = measure(hybrid_one_trial_normal, N, D, noise, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            normal[seed][sigma] = K_mean
            with open(save_path, "a", encoding="utf-8") as f:
                f.write(f"normal  {seed}    {sigma:.3f}  {K_mean:.4f}\n")
            print(f"  normal  σ={sigma:.3f}  K={K_mean:.4f}  ({time.time()-t_cell:.0f}s)",
                  flush=True)

        # Thinned
        thinned[seed] = {}
        for sigma in SIGMAS:
            noise = {} if sigma == 0.0 else {"phase_sigma": sigma}
            t_cell = time.time()
            Ks = measure(hybrid_one_trial_thinned, N, D, noise, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            thinned[seed][sigma] = K_mean
            with open(save_path, "a", encoding="utf-8") as f:
                f.write(f"thinned {seed}    {sigma:.3f}  {K_mean:.4f}\n")
            print(f"  thinned σ={sigma:.3f}  K={K_mean:.4f}  ({time.time()-t_cell:.0f}s)",
                  flush=True)

        # Per-seed SR
        for sigma in SIGMAS:
            if sigma == 0.0:
                continue
            n_base = normal[seed][0.0]
            n_sig = normal[seed][sigma]
            t_base = thinned[seed][0.0]
            t_sig = thinned[seed][sigma]
            n_sr = (n_base - n_sig) / n_base * 100 if n_base > 0 else 0
            t_sr = (t_base - t_sig) / t_base * 100 if t_base > 0 else 0
            amp = t_sr / n_sr if abs(n_sr) > 0.01 else float('inf')
            print(f"  σ={sigma:.3f}  normal SR={n_sr:+.2f}%  thinned SR={t_sr:+.2f}%  "
                  f"amp={amp:+.1f}x", flush=True)
        print(flush=True)

    # Final analysis
    print(f"━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

    # K_base comparison
    n_K_bases = [normal[s][0.0] for s in range(1, N_SEEDS + 1)]
    t_K_bases = [thinned[s][0.0] for s in range(1, N_SEEDS + 1)]
    print(f"\nK_baseline (σ=0):")
    print(f"  normal:  mean={statistics.mean(n_K_bases):.4f}  vals={n_K_bases}")
    print(f"  thinned: mean={statistics.mean(t_K_bases):.4f}  vals={t_K_bases}")
    print(f"  → thinned 의 sub-functional 정도: {statistics.mean(t_K_bases) / statistics.mean(n_K_bases):.2f}x")

    # SR comparison
    print(f"\nSR % comparison:")
    print(f"  {'σ':>7}  {'normal mean':>12}  {'thinned mean':>13}  {'amplification':>14}")
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        n_srs = []
        t_srs = []
        for s in range(1, N_SEEDS + 1):
            n_base = normal[s][0.0]
            n_sig = normal[s][sigma]
            t_base = thinned[s][0.0]
            t_sig = thinned[s][sigma]
            if n_base > 0:
                n_srs.append((n_base - n_sig) / n_base * 100)
            if t_base > 0:
                t_srs.append((t_base - t_sig) / t_base * 100)
        n_mean = statistics.mean(n_srs) if n_srs else 0
        t_mean = statistics.mean(t_srs) if t_srs else 0
        amp = t_mean / n_mean if abs(n_mean) > 0.01 else float('inf')
        marker = ""
        if abs(t_mean) > 5:
            marker = " ★★ ENAQT-style"
        elif abs(t_mean) > 2 * abs(n_mean) and abs(n_mean) > 0.1:
            marker = " ★ amplified"
        print(f"  {sigma:>7.3f}  {n_mean:>+11.3f}%  {t_mean:>+12.3f}%  {amp:>+13.1f}x{marker}")

    print(f"\n총 시간: {time.time()-t_start:.0f}s")
    print(f"결과 저장: {save_path}")


if __name__ == "__main__":
    main()
