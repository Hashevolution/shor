"""
sr_mild_amplification.py — Mild thinned hybrid SR test.

이전 실험 (sr_amplification.py) 의 thinned 가 *과도하게* sub-functional
(K_base=19.81, 거의 max_runs) → noise 효과 0.

본 실험 의 mild thinned:
  - ALL convergent denominators 유지 (이전: smallest 1개만)
  - (C) divisor search 비활성 (state.L 누적 안 함)
  - factor_from_exponent 비활성 (end-of-run fast path 안 함)
  → 알고리즘 의 augmentation 만 제거
  → K_base 예상: 5-15 (functional 하지만 slower)
  → Borderline trial 유지 → noise effect 가능 ★

실험 설계:
  - N = 437, d = 4 (정상 cell)
  - σ ∈ {0.000, 0.050, 0.150} (3 값)
  - 3 seeds × 100 trials
  - Normal hybrid (reference) + Mild thinned 둘 다

예상 결과 시나리오:
  - K_base normal ~2, mild thinned ~5-10 (확인 필요)
  - 만약 mild thinned 의 SR > normal 의 SR → amplification ★
  - 만약 mild thinned 의 SR ~ normal → mechanism 이 base set 결정
  - 만약 mild thinned 도 SR 0 → noise 가 augmentation 통해 작동했음을 시사

시간 추정 (N=437 d=4):
  - 정상 hybrid: ~0.3s/trial × 900 trials = ~5분
  - Mild thinned: K_base 5-10 면 ~2s/trial × 900 = ~30분
  - 총 ~35-40분

실행:
  python -u -m experiments.sr_mild_amplification
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


def hybrid_one_trial_mild_thinned(N, d, noise_kwargs, seed, max_runs=20):
    """Mild thinned: ALL convergents but NO (C) augmentation.

    제거된 것:
      - state.L 의 divisors() candidate 추가 (NO divisor search)
      - factor_from_exponent end-of-run check (NO fast path)

    유지된 것:
      - ALL convergent denominators per measurement
      - b-trick (ord 발견 시 b_i^ord 검사)
      - state.L 업데이트 (다음 측정에 영향 없음)
    """
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)
    state = MultiBaseState()

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs)
        for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
            # ALL convergent denominators (mild thinned: 유지)
            cands = set(convergent_denominators(ki, Q, N - 1))

            # NO divisor search (mild thinned: 제거)
            # if state.L > 1:
            #     cands.update(divisors(state.L))

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

        # NO factor_from_exponent (mild thinned: 제거)
        # if state.L > 1:
        #     rng_f = random.Random(seed)
        #     res = factor_from_exponent(N, state.L, rng_f, max_attempts=5)
        #     if res and 1 < res.factor < N:
        #         return K

    return max_runs


def hybrid_one_trial_normal(N, d, noise_kwargs, seed, max_runs=20):
    """Reference normal hybrid (with all (C) augmentations)."""
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
    print(f"# SR *Mild* Amplification Test")
    print(f"# Mild thinned hybrid (ALL convergents, NO (C) augmentation)")
    print(f"# vs Normal hybrid (full hybrid)")
    print(f"# N={N} d={D}, σ ∈ {SIGMAS}, {N_SEEDS} seeds × {TRIALS} trials")
    print(flush=True)

    t_start = time.time()
    save_path = Path("experiments/sr_mild_amplification_results.txt")
    with open(save_path, "w", encoding="utf-8") as f:
        f.write(f"# SR mild amplification at N={N} d={D}\n")
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"mode    seed  sigma   K_mean\n")

    normal: dict[int, dict[float, float]] = {}
    mild: dict[int, dict[float, float]] = {}

    for seed in range(1, N_SEEDS + 1):
        print(f"━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

        # Normal
        normal[seed] = {}
        for sigma in SIGMAS:
            noise = {} if sigma == 0.0 else {"phase_sigma": sigma}
            t_cell = time.time()
            Ks = measure(hybrid_one_trial_normal, N, D, noise, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            normal[seed][sigma] = K_mean
            with open(save_path, "a", encoding="utf-8") as f:
                f.write(f"normal  {seed}    {sigma:.3f}  {K_mean:.4f}\n")
            print(f"  normal       σ={sigma:.3f}  K={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"({time.time()-t_cell:.0f}s)", flush=True)

        # Mild thinned
        mild[seed] = {}
        for sigma in SIGMAS:
            noise = {} if sigma == 0.0 else {"phase_sigma": sigma}
            t_cell = time.time()
            Ks = measure(hybrid_one_trial_mild_thinned, N, D, noise, TRIALS, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            mild[seed][sigma] = K_mean
            with open(save_path, "a", encoding="utf-8") as f:
                f.write(f"mild    {seed}    {sigma:.3f}  {K_mean:.4f}\n")
            print(f"  mild thinned σ={sigma:.3f}  K={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"({time.time()-t_cell:.0f}s)", flush=True)

        # Per-seed SR comparison
        for sigma in SIGMAS:
            if sigma == 0.0:
                continue
            n_base = normal[seed][0.0]
            n_sig = normal[seed][sigma]
            m_base = mild[seed][0.0]
            m_sig = mild[seed][sigma]
            n_sr = (n_base - n_sig) / n_base * 100 if n_base > 0 else 0
            m_sr = (m_base - m_sig) / m_base * 100 if m_base > 0 else 0
            amp = m_sr / n_sr if abs(n_sr) > 0.01 else float('inf')
            print(f"  σ={sigma:.3f}  normal SR={n_sr:+.2f}%  "
                  f"mild SR={m_sr:+.2f}%  amp={amp:+.2f}x", flush=True)
        print(flush=True)

    # Final analysis
    print(f"━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

    n_K_bases = [normal[s][0.0] for s in range(1, N_SEEDS + 1)]
    m_K_bases = [mild[s][0.0] for s in range(1, N_SEEDS + 1)]
    print(f"\nK_baseline (σ=0):")
    print(f"  normal:       mean={statistics.mean(n_K_bases):.4f}  vals={n_K_bases}")
    print(f"  mild thinned: mean={statistics.mean(m_K_bases):.4f}  vals={m_K_bases}")
    print(f"  → sub-functional 정도: {statistics.mean(m_K_bases) / statistics.mean(n_K_bases):.2f}x")

    print(f"\nSR % comparison:")
    print(f"  {'σ':>7}  {'normal mean':>12}  {'mild mean':>10}  {'amplification':>14}")
    for sigma in SIGMAS:
        if sigma == 0.0:
            continue
        n_srs, m_srs = [], []
        for s in range(1, N_SEEDS + 1):
            if normal[s][0.0] > 0:
                n_srs.append((normal[s][0.0] - normal[s][sigma]) / normal[s][0.0] * 100)
            if mild[s][0.0] > 0:
                m_srs.append((mild[s][0.0] - mild[s][sigma]) / mild[s][0.0] * 100)
        n_mean = statistics.mean(n_srs) if n_srs else 0
        m_mean = statistics.mean(m_srs) if m_srs else 0
        amp = m_mean / n_mean if abs(n_mean) > 0.01 else float('inf')
        marker = ""
        if abs(m_mean) > 5:
            marker = " ★★ ENAQT-style"
        elif abs(m_mean) > 2 * abs(n_mean) and abs(n_mean) > 0.1:
            marker = " ★ amplified"
        elif abs(m_mean) < 0.2 and m_K_bases[0] > 15:
            marker = " (over-thinned)"
        print(f"  {sigma:>7.3f}  {n_mean:>+11.3f}%  {m_mean:>+9.3f}%  {amp:>+13.2f}x{marker}")

    # Direction comparison
    print(f"\nDirection (σ=0.050):")
    n_dirs = collections.Counter()
    m_dirs = collections.Counter()
    for s in range(1, N_SEEDS + 1):
        n_sr = (normal[s][0.0] - normal[s][0.050]) / normal[s][0.0] * 100 if normal[s][0.0] > 0 else 0
        m_sr = (mild[s][0.0] - mild[s][0.050]) / mild[s][0.0] * 100 if mild[s][0.0] > 0 else 0
        n_dirs["+" if n_sr > 0 else ("-" if n_sr < 0 else "0")] += 1
        m_dirs["+" if m_sr > 0 else ("-" if m_sr < 0 else "0")] += 1
    print(f"  normal:       +{n_dirs['+']}  -{n_dirs['-']}  0{n_dirs['0']}")
    print(f"  mild thinned: +{m_dirs['+']}  -{m_dirs['-']}  0{m_dirs['0']}")

    print(f"\n총 시간: {time.time()-t_start:.0f}s")
    print(f"결과 저장: {save_path}")


if __name__ == "__main__":
    main()
