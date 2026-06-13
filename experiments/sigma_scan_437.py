"""
σ scan at N=437 d=4 — V3 anchor cell 의 정밀 σ → SR curve mapping.

목적:
  V3 의 single robust cell (N=437 d=4, SR ≈ +0.91% at σ=0.05) 에서
  σ → K(σ) curve 를 multi-seed paired analysis 로 정밀 측정.
  σ_opt 위치 + curve shape + per-σ statistical significance 확정.

설계:
  - σ grid: 12 값 (0, 0.005, ..., 0.200) — σ_opt ≈ 0.01-0.05 근처 dense
  - paired: 각 seed × trial_index 에 대해 같은 base_seed → 같은 base set 사용
            σ만 바뀌므로 K[σ=0,t] 와 K[σ=x,t] 는 paired observation
  - 분석:
    * per-seed: mean K(σ) — between-seed sd 추정
    * paired (pooled): 각 trial 의 K0-Kσ → 작은 SE 달성
  - 실시간 progress + ETA + 저장

실행:
    python -u -m experiments.sigma_scan_437                # 3 seeds × 200 trials × 12 σ
    python -u -m experiments.sigma_scan_437 5              # 5 seeds × 200 trials
    python -u -m experiments.sigma_scan_437 3 300          # 3 seeds × 300 trials
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


N = 437
D = 4

# σ grid: σ_opt 추정 영역 (0.01-0.05) 근처 dense, 양쪽 tail 도 포함
SIGMA_GRID = [
    0.000,   # baseline
    0.005,
    0.010,   # ★ V3 σ scan 의 σ_opt 후보
    0.015,
    0.020,
    0.025,
    0.035,
    0.050,   # ★ V3 의 robust σ
    0.075,
    0.100,
    0.150,
    0.200,
]


def hybrid_one_trial(N: int, d: int, noise_kwargs: dict, seed: int, max_runs: int = 20):
    """한 trial 의 hybrid 실행. K (성공까지 runs) 반환."""
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


def measure_cell(N, d, noise_kwargs, trials, base_seed):
    """한 (seed, σ) cell 측정. 각 trial 의 K 를 list 로 반환."""
    Ks = []
    for t in range(trials):
        trial_seed = base_seed + t * 1000  # paired: 같은 (base_seed, t) → 같은 base set
        K = hybrid_one_trial(N, d, noise_kwargs, seed=trial_seed)
        Ks.append(K)
    return Ks


def fmt_row(label, mean_sr, se_sr, t_val, per_seed):
    """한 줄 표 format."""
    marker = ""
    if t_val is not None and not math.isnan(t_val):
        if abs(t_val) > 3:
            marker = " ★★★"
        elif abs(t_val) > 2:
            marker = " ★★"
        elif abs(t_val) > 1.5:
            marker = " ★"
    seeds_str = " ".join(f"{s:+5.2f}" for s in per_seed)
    return (f"  {label:>7}  {mean_sr:>+8.3f}%  {se_sr:>6.3f}  "
            f"{t_val:>+6.2f}  [{seeds_str}]{marker}")


def main(argv):
    n_seeds = int(argv[1]) if len(argv) > 1 else 3
    trials = int(argv[2]) if len(argv) > 2 else 200

    total_cells = n_seeds * len(SIGMA_GRID)
    total_trials = total_cells * trials

    print(f"# σ scan at N={N} d={D} — V3 anchor cell precise mapping")
    print(f"# {n_seeds} seeds × {trials} trials × {len(SIGMA_GRID)} σ values")
    print(f"# Total: {total_trials} trials ({total_cells} cells)")
    print(f"# σ grid: {SIGMA_GRID}")
    print(f"# 예상: cell 당 ~{trials * 0.7:.0f}s (N=437 작음) → 전체 ~{total_cells * trials * 0.7 / 60:.0f} 분")
    print(flush=True)

    # results[seed][sigma] = list of per-trial Ks
    results: dict = {seed: {} for seed in range(1, n_seeds + 1)}

    t_global = time.time()
    cell_idx = 0

    for seed in range(1, n_seeds + 1):
        print(f"━━ seed {seed}/{n_seeds} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

        for sigma in SIGMA_GRID:
            cell_idx += 1
            t_cell = time.time()
            noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}

            Ks = measure_cell(N, D, noise_kwargs, trials, base_seed=seed)
            results[seed][sigma] = Ks

            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            cell_elapsed = time.time() - t_cell
            global_elapsed = time.time() - t_global
            eta = global_elapsed * (total_cells - cell_idx) / cell_idx

            # 진행 표시
            print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"cell={cell_elapsed:>4.0f}s  total ETA={eta:>5.0f}s "
                  f"({cell_idx:>3}/{total_cells})",
                  flush=True)

        # seed 종료시 quick summary
        K_base = statistics.mean(results[seed][0.0])
        print(f"\n  seed {seed} per-σ SR (K_base = {K_base:.4f}):", flush=True)
        for sigma in SIGMA_GRID:
            if sigma == 0.0:
                continue
            Kσ = statistics.mean(results[seed][sigma])
            sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
            mark = " ↓" if sr > 0.3 else (" ↑" if sr < -0.3 else "")
            print(f"    σ={sigma:.3f}  Kσ={Kσ:.4f}  SR={sr:+.3f}%{mark}", flush=True)
        print(flush=True)

    # ━━━ 최종 분석 (paired + between-seed) ━━━
    print(f"━━ FINAL σ → SR curve ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", flush=True)

    # K_baseline pooled across seeds
    K_base_per_seed = [statistics.mean(results[s][0.0]) for s in range(1, n_seeds + 1)]
    K_base_pooled = statistics.mean(K_base_per_seed)
    K_base_se = (statistics.stdev(K_base_per_seed) / math.sqrt(n_seeds)
                 if n_seeds > 1 else 0)
    print(f"  K_baseline (σ=0) pooled = {K_base_pooled:.4f}  "
          f"per-seed: {K_base_per_seed}\n", flush=True)

    # Between-seed analysis (per-seed mean SR)
    print(f"  ─── Between-seed analysis ───")
    print(f"  {'σ':>7}  {'mean SR%':>9}  {'SE':>6}  {'t':>6}  {'per-seed SR%':<35}")
    print(f"  {'─'*7}  {'─'*9}  {'─'*6}  {'─'*6}  {'─'*35}")

    between_seed_results = {}
    for sigma in SIGMA_GRID:
        if sigma == 0.0:
            continue
        per_seed_srs = []
        for seed in range(1, n_seeds + 1):
            K_base = statistics.mean(results[seed][0.0])
            Kσ = statistics.mean(results[seed][sigma])
            sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
            per_seed_srs.append(sr)

        mean_sr = statistics.mean(per_seed_srs)
        sd_sr = statistics.stdev(per_seed_srs) if n_seeds > 1 else 0
        se_sr = sd_sr / math.sqrt(n_seeds) if n_seeds > 1 else 0
        t_val = mean_sr / se_sr if se_sr > 0 else float('nan')
        between_seed_results[sigma] = (mean_sr, sd_sr, se_sr, t_val, per_seed_srs)

        print(fmt_row(f"{sigma:.3f}", mean_sr, se_sr, t_val, per_seed_srs))

    # Paired analysis (pooled per-trial diffs)
    print(f"\n  ─── Paired analysis (pooled across seeds) ───")
    print(f"  {'σ':>7}  {'mean SR%':>9}  {'SE':>6}  {'t':>6}  {'n_paired':<10}")
    print(f"  {'─'*7}  {'─'*9}  {'─'*6}  {'─'*6}  {'─'*10}")

    paired_results = {}
    for sigma in SIGMA_GRID:
        if sigma == 0.0:
            continue
        all_diffs = []
        for seed in range(1, n_seeds + 1):
            K0_list = results[seed][0.0]
            Kσ_list = results[seed][sigma]
            for k0, ks in zip(K0_list, Kσ_list):
                all_diffs.append(k0 - ks)

        n_paired = len(all_diffs)
        mean_diff = statistics.mean(all_diffs)
        sd_diff = statistics.stdev(all_diffs) if n_paired > 1 else 0
        se_diff = sd_diff / math.sqrt(n_paired) if n_paired > 1 else 0

        sr_paired = mean_diff / K_base_pooled * 100 if K_base_pooled > 0 else 0
        sr_paired_se = se_diff / K_base_pooled * 100 if K_base_pooled > 0 else 0
        t_paired = sr_paired / sr_paired_se if sr_paired_se > 0 else float('nan')
        paired_results[sigma] = (sr_paired, sd_diff, se_diff, t_paired, n_paired)

        marker = ""
        if not math.isnan(t_paired):
            if abs(t_paired) > 3:
                marker = " ★★★"
            elif abs(t_paired) > 2:
                marker = " ★★"
            elif abs(t_paired) > 1.5:
                marker = " ★"
        print(f"  {sigma:>7.3f}  {sr_paired:>+8.3f}%  {sr_paired_se:>6.3f}  "
              f"{t_paired:>+6.2f}  n={n_paired:<8}{marker}",
              flush=True)

    # σ_opt 추정 (paired analysis 기준)
    pos_sigma = [(s, paired_results[s][0]) for s in paired_results
                 if paired_results[s][0] > 0]
    if pos_sigma:
        best_sigma, best_sr = max(pos_sigma, key=lambda x: x[1])
        best_t = paired_results[best_sigma][3]
        print(f"\n  σ_opt (max SR, paired) = {best_sigma:.3f}  "
              f"→ SR = {best_sr:+.3f}%  (t = {best_t:+.2f})")
    else:
        print(f"\n  σ_opt: 모든 σ 에서 SR ≤ 0 — robust 영역 미검출")

    # 저장
    save_path = "experiments/sigma_scan_437_d4_results.txt"
    with open(save_path, "w") as f:
        f.write(f"# σ scan at N={N} d={D}\n")
        f.write(f"# {n_seeds} seeds × {trials} trials × {len(SIGMA_GRID)} σ\n")
        f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Total time: {time.time() - t_global:.0f}s\n\n")
        f.write(f"# K_baseline pooled = {K_base_pooled:.4f}\n\n")
        f.write(f"# Per-cell K_mean\n")
        f.write(f"sigma   seed  K_mean\n")
        for seed in range(1, n_seeds + 1):
            for sigma in SIGMA_GRID:
                K_mean = statistics.mean(results[seed][sigma])
                f.write(f"{sigma:.3f}   {seed}     {K_mean:.4f}\n")
        f.write(f"\n# Paired analysis (pooled per-trial diffs)\n")
        f.write(f"sigma   SR%      SE      t       n_paired\n")
        for sigma in SIGMA_GRID:
            if sigma == 0.0:
                continue
            sr, sd_d, se, t, n = paired_results[sigma]
            f.write(f"{sigma:.3f}   {sr:+.3f}  {se:.3f}  {t:+.2f}  {n}\n")

    print(f"\n  결과 저장: {save_path}")
    print(f"  총 소요: {time.time() - t_global:.0f}s ({(time.time() - t_global) / 60:.1f} 분)")


if __name__ == "__main__":
    main(sys.argv)
