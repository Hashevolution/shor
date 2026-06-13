"""
N=1147 d=1 confirm — *실시간 진행 출력* verbose 버전.

매 trial 마다 진행 표시.
- 5 trial 마다: 진행률 + 현재 평균 K
- 매 seed 종료시: SR 결과

실행:
    python -u -m experiments.confirm_verbose
또는 N, d, trials 조절:
    python -u -m experiments.confirm_verbose 1147 1 200
"""

from __future__ import annotations
import math
import random
import statistics
import sys
import time

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState, convergent_denominators, divisors, minimize_order,
    factor_from_exponent,
)
from experiments.rv_filter_lll import regev_setup_bases, simulate_regev_run


def hybrid_one_trial(N: int, d: int, noise_kwargs: dict, seed: int, max_runs: int = 20):
    """한 trial 의 hybrid 실행. K (성공까지 runs) 반환."""
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    setup = regev_setup_bases(N, d, rng_py)
    state = MultiBaseState()

    for K in range(1, max_runs + 1):
        run = simulate_regev_run(
            setup.a, N, Q, rng_np, noise_kwargs=noise_kwargs,
        )
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


def measure_verbose(N: int, d: int, noise_kwargs: dict, trials: int, base_seed: int):
    """실시간 진행 표시하며 측정."""
    Ks = []
    t_start = time.time()
    for t in range(trials):
        K = hybrid_one_trial(N, d, noise_kwargs, seed=base_seed + t * 1000)
        Ks.append(K)

        # 매 5 trial 마다 진행 표시
        if (t + 1) % 5 == 0 or (t + 1) == trials:
            elapsed = time.time() - t_start
            rate = (t + 1) / elapsed if elapsed > 0 else 0
            remaining_s = (trials - t - 1) / rate if rate > 0 else 0
            mean_K = sum(Ks) / len(Ks)
            print(f"      [trial {t+1:>4}/{trials}] K_mean={mean_K:.3f}  "
                  f"elapsed={elapsed:>5.0f}s  eta={remaining_s:>4.0f}s",
                  flush=True)
    return sum(Ks) / len(Ks)


def main(argv):
    N = int(argv[1]) if len(argv) > 1 else 1147
    d = int(argv[2]) if len(argv) > 2 else 1
    trials = int(argv[3]) if len(argv) > 3 else 100
    n_seeds = int(argv[4]) if len(argv) > 4 else 3

    print(f"# Verbose confirm: N={N}, d={d}, trials={trials}, seeds={n_seeds}")
    print(f"# 각 5 trial 마다 진행 출력\n", flush=True)

    srs = []
    K0s = []
    for seed in range(1, n_seeds + 1):
        print(f"━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)
        print(f"  [σ=0   ] 측정 중...", flush=True)
        K0 = measure_verbose(N, d, {}, trials, seed)
        print(f"  → K(σ=0) = {K0:.3f}", flush=True)

        print(f"  [σ=.05 ] 측정 중...", flush=True)
        K1 = measure_verbose(N, d, {"phase_sigma": 0.05}, trials, seed)
        print(f"  → K(σ=.05) = {K1:.3f}", flush=True)

        sr = (K0 - K1) / K0 * 100 if K0 > 0 else 0
        marker = " ★" if sr > 2 else (" ↓" if sr > 0.5 else (" ↑" if sr < -0.5 else ""))
        print(f"  ▶ seed {seed} SR = {sr:+.2f}%{marker}\n", flush=True)
        srs.append(sr)
        K0s.append(K0)

    print(f"\n━━ 최종 결과 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)
    mean_sr = statistics.mean(srs)
    sd_sr = statistics.stdev(srs) if len(srs) > 1 else 0
    print(f"  Mean SR     = {mean_sr:+.2f}%  (sd={sd_sr:.2f})")
    print(f"  Mean K_base = {statistics.mean(K0s):.3f}")
    n_positive = sum(1 for sr in srs if sr > 0)
    print(f"  Positive    : {n_positive}/{n_seeds}")
    print(f"  All SRs     : {srs}")


if __name__ == "__main__":
    main(sys.argv)
