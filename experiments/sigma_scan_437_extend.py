"""
sigma_scan_437_extend.py — N=437 d=4 σ scan 추가 seed 진행 + 결합 분석.

목적:
  기존 seeds 1-3 결과 (sigma_scan_437_d4_results.txt) 를 *유지* 한 채
  새 seeds 추가 → 통계 검정력 강화 → SR finding *결정적* 통과/탈락.

설계:
  - 기존 K_means 자동 로드 (seeds 1-3 재실행 안 함)
  - 새 seeds 4-13 (또는 사용자 지정 범위) 추가 실행
  - 각 seed 완료시 *즉시* 저장 (중단되어도 결과 보존)
  - 매 seed 완료시 *현재까지* 결합 통계 표시 (early stopping 판단 가능)
  - 종료시 최종 between-seed + sign test + 정직한 conclusion

실행:
  python -u -m experiments.sigma_scan_437_extend                  # 기본: seeds 4-13 (10 추가, ~90분)
  python -u -m experiments.sigma_scan_437_extend 4 15             # seeds 4-15 (12 추가)
  python -u -m experiments.sigma_scan_437_extend 4 20             # seeds 4-20 (17 추가, ~2.5시간)
  python -u -m experiments.sigma_scan_437_extend 4 13 300         # 300 trials per seed (더 정밀)
"""

from __future__ import annotations
import math
import statistics
import sys
import time
from math import comb, erfc, sqrt
from pathlib import Path

from experiments.sigma_scan_437 import measure_cell, SIGMA_GRID, N, D


BASE_FILE = Path("experiments/sigma_scan_437_d4_results.txt")
EXTEND_FILE = Path("experiments/sigma_scan_437_d4_extended.txt")


def read_K_means_from_files():
    """기존 + 확장 파일에서 K_means 로드. K_means[seed][sigma] = K_mean."""
    K_means: dict[int, dict[float, float]] = {}
    for path in [BASE_FILE, EXTEND_FILE]:
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("sigma"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    sigma = float(parts[0])
                    seed = int(parts[1])
                    K_mean = float(parts[2])
                except ValueError:
                    continue
                K_means.setdefault(seed, {})[sigma] = K_mean
    return K_means


def init_extend_file():
    """확장 파일 헤더 생성 (없을 때만)."""
    if EXTEND_FILE.exists():
        return
    with open(EXTEND_FILE, "w", encoding="utf-8") as f:
        f.write(f"# σ scan extension at N={N} d={D}\n")
        f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Per-cell K_mean (append-only)\n\n")
        f.write(f"sigma   seed  K_mean\n")


def append_seed_K_means(seed: int, K_means: dict[float, float]):
    """한 seed 의 결과를 확장 파일에 즉시 append."""
    with open(EXTEND_FILE, "a", encoding="utf-8") as f:
        for sigma in sorted(K_means):
            f.write(f"{sigma:.3f}   {seed}     {K_means[sigma]:.4f}\n")


def p_value_normal(t: float) -> float:
    """t → 1-sided p (normal approximation)."""
    if math.isnan(t):
        return float('nan')
    return 0.5 * erfc(t / sqrt(2))


def report_combined_analysis(all_K_means: dict, label: str = ""):
    """현재 데이터로 between-seed 분석 출력."""
    all_seeds = sorted(all_K_means.keys())
    n_total = len(all_seeds)
    if n_total == 0:
        print("  데이터 없음")
        return

    K_bases = [all_K_means[s].get(0.0, 0.0) for s in all_seeds if 0.0 in all_K_means[s]]
    K_base_mean = statistics.mean(K_bases) if K_bases else 0
    K_base_sd = statistics.stdev(K_bases) if len(K_bases) > 1 else 0

    print(f"\n  {'='*70}")
    print(f"  결합 분석 {label} (n_seeds = {n_total})")
    print(f"  {'='*70}")
    print(f"  K_baseline: mean={K_base_mean:.4f}  sd={K_base_sd:.4f}  "
          f"range=[{min(K_bases):.3f}, {max(K_bases):.3f}]")

    # Per-σ between-seed analysis
    print(f"\n  Between-seed table:")
    print(f"  {'σ':>7}  {'mean SR%':>10}  {'sd':>7}  {'SE':>7}  {'t':>7}  "
          f"{'p (1-sd)':>9}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*9}")

    sigma_stats = {}
    for sigma in SIGMA_GRID:
        if sigma == 0.0:
            continue
        per_seed_srs = []
        for seed in all_seeds:
            if 0.0 not in all_K_means[seed] or sigma not in all_K_means[seed]:
                continue
            K_base = all_K_means[seed][0.0]
            Kσ = all_K_means[seed][sigma]
            if K_base > 0:
                per_seed_srs.append((K_base - Kσ) / K_base * 100)
        if not per_seed_srs:
            continue
        n_s = len(per_seed_srs)
        mean_sr = statistics.mean(per_seed_srs)
        sd_sr = statistics.stdev(per_seed_srs) if n_s > 1 else 0
        se_sr = sd_sr / math.sqrt(n_s) if n_s > 1 else 0
        t_val = mean_sr / se_sr if se_sr > 0 else float('nan')
        p_val = p_value_normal(t_val)
        sigma_stats[sigma] = (mean_sr, sd_sr, se_sr, t_val, p_val, per_seed_srs)

        marker = ""
        if not math.isnan(t_val):
            if abs(t_val) > 3:
                marker = " ★★★"
            elif abs(t_val) > 2:
                marker = " ★★"
            elif abs(t_val) > 1.5:
                marker = " ★"
        print(f"  {sigma:>7.3f}  {mean_sr:>+9.3f}%  {sd_sr:>7.3f}  "
              f"{se_sr:>7.3f}  {t_val:>+7.2f}  {p_val:>9.4f}{marker}")

    # Sign test at σ=0.050 (most representative)
    plateau_sigma = 0.050
    pos = sum(1 for s in all_seeds
              if 0.0 in all_K_means[s] and plateau_sigma in all_K_means[s]
              and all_K_means[s][0.0] > all_K_means[s][plateau_sigma])
    neg = sum(1 for s in all_seeds
              if 0.0 in all_K_means[s] and plateau_sigma in all_K_means[s]
              and all_K_means[s][0.0] < all_K_means[s][plateau_sigma])
    zero_dir = n_total - pos - neg

    print(f"\n  Sign test at σ={plateau_sigma} (plateau center):")
    print(f"    positive direction (noise helps): {pos}/{n_total}")
    print(f"    negative direction (noise hurts): {neg}/{n_total}")
    print(f"    zero / undefined: {zero_dir}/{n_total}")
    if n_total > 0 and (pos + neg) > 0:
        n_pn = pos + neg
        p_sign = sum(comb(n_pn, k) / 2**n_pn for k in range(pos, n_pn + 1))
        print(f"    P(≥{pos}/{n_pn} positive | null 50/50) = {p_sign:.4f}")

    # Per-seed plateau SR list
    print(f"\n  Per-seed plateau SR (σ=0.050):")
    sr_list = []
    for seed in all_seeds:
        if 0.0 in all_K_means[seed] and plateau_sigma in all_K_means[seed]:
            K_base = all_K_means[seed][0.0]
            Kσ = all_K_means[seed][plateau_sigma]
            sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
            sr_list.append(sr)
            mark = "+" if sr > 0 else ("-" if sr < 0 else " ")
            print(f"    seed {seed:>2}: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
                  f"SR={sr:+.3f}%  {mark}")

    if len(sr_list) >= 2:
        print(f"\n  Plateau SR distribution: "
              f"mean={statistics.mean(sr_list):+.3f}%, "
              f"sd={statistics.stdev(sr_list):.3f}%, "
              f"range=[{min(sr_list):+.2f}, {max(sr_list):+.2f}]")


def main(argv):
    start_seed = int(argv[1]) if len(argv) > 1 else 4
    end_seed = int(argv[2]) if len(argv) > 2 else 13
    trials = int(argv[3]) if len(argv) > 3 else 200

    # 기존 데이터 로드
    existing = read_K_means_from_files()
    print(f"# σ scan extension at N={N} d={D}")
    print(f"# 기존 seeds: {sorted(existing.keys())} ({len(existing)} 개)")

    requested = list(range(start_seed, end_seed + 1))
    to_run = [s for s in requested if s not in existing]
    if not to_run:
        print(f"# 모든 seeds {start_seed}..{end_seed} 이미 계산됨. 결합 분석만 출력.")
        report_combined_analysis(existing, label="(기존 데이터만)")
        return

    print(f"# 실행할 새 seeds: {to_run}")
    print(f"# 각 seed: {trials} trials × {len(SIGMA_GRID)} σ values")
    print(f"# 예상: ~9 분/seed → 총 ~{len(to_run) * 9} 분")
    print(flush=True)

    init_extend_file()

    t_global = time.time()
    total_cells = len(to_run) * len(SIGMA_GRID)
    cell_idx = 0

    for seed in to_run:
        print(f"━━ seed {seed} (start) ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
              flush=True)
        seed_K_means: dict[float, float] = {}
        t_seed = time.time()

        for sigma in SIGMA_GRID:
            cell_idx += 1
            t_cell = time.time()
            noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}
            Ks = measure_cell(N, D, noise_kwargs, trials, base_seed=seed)
            K_mean = statistics.mean(Ks)
            K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0
            seed_K_means[sigma] = K_mean

            elapsed = time.time() - t_cell
            global_elapsed = time.time() - t_global
            eta = global_elapsed * (total_cells - cell_idx) / cell_idx

            print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                  f"cell={elapsed:>4.0f}s  ETA={eta:>5.0f}s "
                  f"({cell_idx:>3}/{total_cells})", flush=True)

        # 즉시 저장 (중단 대비)
        append_seed_K_means(seed, seed_K_means)

        # Per-seed summary
        K_base = seed_K_means[0.0]
        plateau_K = seed_K_means[0.050]
        sr = (K_base - plateau_K) / K_base * 100 if K_base > 0 else 0
        seed_time = time.time() - t_seed
        mark = " +" if sr > 0 else (" -" if sr < 0 else "")
        print(f"\n  seed {seed} 완료 ({seed_time:.0f}s):  "
              f"K_base={K_base:.4f}  K(σ=.05)={plateau_K:.4f}  "
              f"SR={sr:+.3f}%{mark}", flush=True)

        # 실시간 결합 분석 (early stopping 판단)
        current = read_K_means_from_files()
        report_combined_analysis(current,
                                 label=f"(seed {seed} 추가 후)")
        print(flush=True)

    print(f"\n━━ FINAL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
          flush=True)
    print(f"총 소요: {time.time() - t_global:.0f}s "
          f"({(time.time() - t_global) / 60:.1f} 분)")
    print(f"확장 파일: {EXTEND_FILE}")
    print(f"  → 다음 실행은 자동으로 이 결과를 이어서 사용")


if __name__ == "__main__":
    main(sys.argv)
