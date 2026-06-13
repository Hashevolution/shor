"""
sigma_scan_437_extend.py — N=437 d=4 σ scan 추가 seed + K-histogram + borderline 검증.

목적:
  1. 추가 seeds 진행 → between-seed statistical 검정력 강화
  2. K-distribution histogram 저장 → borderline K-boundary 가설 직접 검증
     - positive SR seed: 어느 K 에서 어디로 trial 이 이동?
     - negative SR seed: 다른 K-boundary 에서 flip?

설계:
  - 기존 seeds 1-3 K_means 자동 로드 (BASE_FILE)
  - 기존 seeds 의 histogram 부족 시 backfill (σ=0, 0.050 만, ~14분)
  - 새 seeds 의 K_mean + histogram 동시 저장
  - 매 cell / seed 즉시 파일 append (중단 보존)
  - 매 seed 완료 시 실시간 결합 분석 출력 (early stopping 판단)
  - 종료시 borderline K-boundary 검증 표

실행:
  python -u -m experiments.sigma_scan_437_extend                  # 기본: seeds 4-13 (10 추가)
  python -u -m experiments.sigma_scan_437_extend 4 15             # seeds 4-15 (12 추가)
  python -u -m experiments.sigma_scan_437_extend 4 20             # seeds 4-20 (17 추가, ~2.5시간)
"""

from __future__ import annotations
import collections
import math
import statistics
import sys
import time
from math import comb, erfc, sqrt
from pathlib import Path

from experiments.sigma_scan_437 import measure_cell, SIGMA_GRID, N, D


BASE_FILE = Path("experiments/sigma_scan_437_d4_results.txt")
EXTEND_FILE = Path("experiments/sigma_scan_437_d4_extended.txt")
HIST_FILE = Path("experiments/sigma_scan_437_d4_histograms.txt")

# Histogram backfill 시 핵심 σ (baseline + plateau center)
BACKFILL_SIGMAS = [0.0, 0.050]


# ────────────────────────────────────────────────────────────
# I/O
# ────────────────────────────────────────────────────────────

def read_K_means():
    """Base + extend 파일에서 K_means 로드."""
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


def read_histograms():
    """K histograms 로드. hist[seed][sigma] = {K: count}."""
    hists: dict[int, dict[float, dict[int, int]]] = {}
    if not HIST_FILE.exists():
        return hists
    with open(HIST_FILE, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("seed"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                seed = int(parts[0])
                sigma = float(parts[1])
                K = int(parts[2])
                count = int(parts[3])
            except ValueError:
                continue
            hists.setdefault(seed, {}).setdefault(sigma, {})[K] = count
    return hists


def init_files():
    """확장 + histogram 파일 헤더 (없을 때만)."""
    if not EXTEND_FILE.exists():
        with open(EXTEND_FILE, "w", encoding="utf-8") as f:
            f.write(f"# σ scan extension at N={N} d={D}\n")
            f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Per-cell K_mean (append-only)\n\n")
            f.write(f"sigma   seed  K_mean\n")
    if not HIST_FILE.exists():
        with open(HIST_FILE, "w", encoding="utf-8") as f:
            f.write(f"# K histograms for σ scan at N={N} d={D}\n")
            f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Borderline K-boundary 가설 검증용\n\n")
            f.write(f"seed\tsigma\tK\tcount\n")


def append_K_mean(seed, sigma, K_mean):
    with open(EXTEND_FILE, "a", encoding="utf-8") as f:
        f.write(f"{sigma:.3f}   {seed}     {K_mean:.4f}\n")


def append_histogram(seed, sigma, hist):
    with open(HIST_FILE, "a", encoding="utf-8") as f:
        for K in sorted(hist):
            f.write(f"{seed}\t{sigma:.3f}\t{K}\t{hist[K]}\n")


# ────────────────────────────────────────────────────────────
# 측정
# ────────────────────────────────────────────────────────────

def Ks_to_histogram(Ks):
    return dict(collections.Counter(Ks))


def measure_with_hist(N, d, noise_kwargs, trials, base_seed):
    """measure_cell wrapper: K_mean + histogram 동시 반환."""
    Ks = measure_cell(N, d, noise_kwargs, trials, base_seed)
    return statistics.mean(Ks), Ks_to_histogram(Ks), Ks


# ────────────────────────────────────────────────────────────
# 분석
# ────────────────────────────────────────────────────────────

def p_value_normal(t):
    if math.isnan(t):
        return float('nan')
    return 0.5 * erfc(t / sqrt(2))


def report_combined_analysis(K_means_all, label=""):
    """Between-seed + sign test 결합 분석."""
    all_seeds = sorted(K_means_all.keys())
    n_total = len(all_seeds)
    if n_total == 0:
        return

    K_bases = [K_means_all[s].get(0.0, 0.0) for s in all_seeds if 0.0 in K_means_all[s]]
    K_base_mean = statistics.mean(K_bases) if K_bases else 0
    K_base_sd = statistics.stdev(K_bases) if len(K_bases) > 1 else 0

    print(f"\n  {'='*72}")
    print(f"  결합 분석 {label} (n_seeds = {n_total})")
    print(f"  {'='*72}")
    print(f"  K_baseline: mean={K_base_mean:.4f}  sd={K_base_sd:.4f}  "
          f"range=[{min(K_bases):.3f}, {max(K_bases):.3f}]")

    print(f"\n  Between-seed table:")
    print(f"  {'σ':>7}  {'mean SR%':>10}  {'sd':>7}  {'SE':>7}  {'t':>7}  "
          f"{'p (1-sd)':>9}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*9}")

    for sigma in SIGMA_GRID:
        if sigma == 0.0:
            continue
        per_seed_srs = []
        for seed in all_seeds:
            if 0.0 not in K_means_all[seed] or sigma not in K_means_all[seed]:
                continue
            K_base = K_means_all[seed][0.0]
            Kσ = K_means_all[seed][sigma]
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

    # Sign test at σ=0.050
    plateau_sigma = 0.050
    pos = sum(1 for s in all_seeds
              if 0.0 in K_means_all[s] and plateau_sigma in K_means_all[s]
              and K_means_all[s][0.0] > K_means_all[s][plateau_sigma])
    neg = sum(1 for s in all_seeds
              if 0.0 in K_means_all[s] and plateau_sigma in K_means_all[s]
              and K_means_all[s][0.0] < K_means_all[s][plateau_sigma])
    print(f"\n  Sign test at σ={plateau_sigma}: "
          f"positive={pos}/{n_total}, negative={neg}/{n_total}")
    if (pos + neg) > 0:
        n_pn = pos + neg
        p_sign = sum(comb(n_pn, k) / 2**n_pn for k in range(pos, n_pn + 1))
        print(f"    P(≥{pos}/{n_pn} positive | null 50/50) = {p_sign:.4f}")

    print(f"\n  Per-seed plateau SR (σ=0.050):")
    sr_list = []
    for seed in all_seeds:
        if 0.0 in K_means_all[seed] and plateau_sigma in K_means_all[seed]:
            K_base = K_means_all[seed][0.0]
            Kσ = K_means_all[seed][plateau_sigma]
            sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0
            sr_list.append(sr)
            mark = "+" if sr > 0 else ("-" if sr < 0 else " ")
            print(f"    seed {seed:>2}: K_base={K_base:.4f}  K(σ=.05)={Kσ:.4f}  "
                  f"SR={sr:+.3f}%  {mark}")


def identify_dominant_flip(h0, h_sigma):
    """diff hist 에서 주요 flip (K_from → K_to) 식별."""
    K_range = sorted(set(h0.keys()) | set(h_sigma.keys()))
    diffs = {K: h_sigma.get(K, 0) - h0.get(K, 0) for K in K_range}
    gains = [(K, d) for K, d in diffs.items() if d > 0]
    losses = [(K, -d) for K, d in diffs.items() if d < 0]
    if not gains or not losses:
        return None
    K_to, gain_mag = max(gains, key=lambda x: x[1])
    K_from, loss_mag = max(losses, key=lambda x: x[1])
    return K_from, K_to, min(gain_mag, loss_mag)


def report_histogram_analysis(hists, K_means_all):
    """Borderline K-boundary 가설 검증."""
    print(f"\n  {'='*72}")
    print(f"  Histogram analysis — borderline K-boundary 가설 검증")
    print(f"  {'='*72}")

    eligible_seeds = sorted(
        s for s in hists
        if 0.0 in hists[s] and 0.050 in hists[s]
        and 0.0 in K_means_all.get(s, {}) and 0.050 in K_means_all.get(s, {})
    )
    if not eligible_seeds:
        print(f"  분석 가능 seed 없음 (histogram 부족).")
        return

    print(f"  분석 대상 seeds: {eligible_seeds}\n")

    # 방향별 분류
    positive_seeds = []
    negative_seeds = []
    for seed in eligible_seeds:
        K_base = K_means_all[seed][0.0]
        K_sigma = K_means_all[seed][0.050]
        sr = (K_base - K_sigma) / K_base * 100 if K_base > 0 else 0
        if sr > 0.01:
            positive_seeds.append(seed)
        elif sr < -0.01:
            negative_seeds.append(seed)

    print(f"  방향별 seed 분류 (σ=0.050 plateau):")
    print(f"    Positive SR (noise helps): {positive_seeds} ({len(positive_seeds)})")
    print(f"    Negative SR (noise hurts): {negative_seeds} ({len(negative_seeds)})")

    # Per-seed flip 분석
    print(f"\n  Per-seed K-flip:")
    print(f"  {'seed':>5} {'K_base':>7} {'SR%':>7} {'dir':>4} "
          f"{'flip':>15} {'magnitude':>10}")
    print(f"  {'-'*5} {'-'*7} {'-'*7} {'-'*4} {'-'*15} {'-'*10}")

    pos_flips = []  # (K_from, K_to) for positive seeds
    neg_flips = []  # (K_from, K_to) for negative seeds

    for seed in eligible_seeds:
        h0 = hists[seed][0.0]
        h_s = hists[seed][0.050]
        K_base = K_means_all[seed][0.0]
        K_sigma = K_means_all[seed][0.050]
        sr = (K_base - K_sigma) / K_base * 100 if K_base > 0 else 0
        direction = "+" if sr > 0 else ("-" if sr < 0 else "=")

        flip = identify_dominant_flip(h0, h_s)
        if flip is None:
            print(f"  {seed:>5} {K_base:>7.3f} {sr:>+6.2f}% {direction:>4} "
                  f"{'no change':>15} {'-':>10}")
            continue

        K_from, K_to, mag = flip
        flip_str = f"K={K_from}→K={K_to}"
        print(f"  {seed:>5} {K_base:>7.3f} {sr:>+6.2f}% {direction:>4} "
              f"{flip_str:>15} {mag:>10}")

        if sr > 0.01:
            pos_flips.append((K_from, K_to))
        elif sr < -0.01:
            neg_flips.append((K_from, K_to))

    # 가설 검증
    print(f"\n  ─── 가설 검증 ───")
    if pos_flips:
        from_Ks = collections.Counter(K_from for K_from, _ in pos_flips)
        to_Ks = collections.Counter(K_to for _, K_to in pos_flips)
        print(f"  Positive seeds 주요 flip: ")
        print(f"    K_from 분포: {dict(from_Ks)}")
        print(f"    K_to 분포:   {dict(to_Ks)}")
        avg_from = statistics.mean(f for f, _ in pos_flips)
        avg_to = statistics.mean(t for _, t in pos_flips)
        print(f"    평균: K={avg_from:.1f} → K={avg_to:.1f}")
        if avg_from > avg_to:
            print(f"    → 가설 지지: positive seeds 는 *높은 K → 낮은 K* (성공 측면)")
    if neg_flips:
        from_Ks = collections.Counter(K_from for K_from, _ in neg_flips)
        to_Ks = collections.Counter(K_to for _, K_to in neg_flips)
        print(f"  Negative seeds 주요 flip: ")
        print(f"    K_from 분포: {dict(from_Ks)}")
        print(f"    K_to 분포:   {dict(to_Ks)}")
        avg_from = statistics.mean(f for f, _ in neg_flips)
        avg_to = statistics.mean(t for _, t in neg_flips)
        print(f"    평균: K={avg_from:.1f} → K={avg_to:.1f}")
        if avg_from < avg_to:
            print(f"    → 가설 지지: negative seeds 는 *낮은 K → 높은 K* (실패 측면)")

    if pos_flips and neg_flips:
        avg_pos_from = statistics.mean(f for f, _ in pos_flips)
        avg_neg_from = statistics.mean(f for f, _ in neg_flips)
        if avg_pos_from > avg_neg_from:
            print(f"\n  ★ 가설 A+B 의 핵심 prediction 만족: ")
            print(f"    Positive flip 의 K_from ({avg_pos_from:.1f}) > "
                  f"Negative flip 의 K_from ({avg_neg_from:.1f})")
            print(f"    → noise effect direction 은 base set 의 borderline K-bin 결정")


# ────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────

def main(argv):
    start_seed = int(argv[1]) if len(argv) > 1 else 4
    end_seed = int(argv[2]) if len(argv) > 2 else 13
    trials = int(argv[3]) if len(argv) > 3 else 200

    K_means_all = read_K_means()
    hists = read_histograms()
    existing_seeds = sorted(K_means_all.keys())

    print(f"# σ scan extension at N={N} d={D}  (Histogram + Borderline 검증)")
    print(f"# 기존 K_means seeds: {existing_seeds}")
    print(f"# 기존 histograms: {sorted(hists.keys())}")

    init_files()

    # Backfill 필요 판단
    backfill_seeds = []
    for seed in existing_seeds:
        if (seed not in hists
                or 0.0 not in hists.get(seed, {})
                or 0.050 not in hists.get(seed, {})):
            backfill_seeds.append(seed)

    if backfill_seeds:
        print(f"# Histogram backfill 필요: {backfill_seeds}")
        print(f"#   σ ∈ {BACKFILL_SIGMAS} 만 측정 (~{len(backfill_seeds) * 2 * 70 / 60:.0f}분 추가)")

    requested = list(range(start_seed, end_seed + 1))
    to_run = [s for s in requested if s not in existing_seeds]
    if to_run:
        print(f"# 신규 seeds: {to_run}")
        print(f"#   예상: ~9분/seed × {len(to_run)} = ~{len(to_run) * 9}분")

    if not backfill_seeds and not to_run:
        print(f"# 작업 없음. 현재 데이터 분석만 출력.")
        report_combined_analysis(K_means_all, "(현재까지)")
        if hists:
            report_histogram_analysis(hists, K_means_all)
        return

    print(flush=True)
    t_global = time.time()

    # ─── BACKFILL ───
    if backfill_seeds:
        print(f"━━ BACKFILL: 기존 seeds histogram 보강 ━━━━━━━━━━━━━━━━━━━━",
              flush=True)
        for seed in backfill_seeds:
            for sigma in BACKFILL_SIGMAS:
                t_cell = time.time()
                noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}
                K_mean, hist, _ = measure_with_hist(
                    N, D, noise_kwargs, trials, base_seed=seed)
                append_histogram(seed, sigma, hist)
                hists.setdefault(seed, {})[sigma] = hist
                elapsed = time.time() - t_cell
                # K_mean 비교 (deterministic check)
                ref = K_means_all.get(seed, {}).get(sigma, K_mean)
                match = "✓" if abs(K_mean - ref) < 0.001 else f"≠{ref:.4f}"
                print(f"  seed {seed}  σ={sigma:.3f}  K_mean={K_mean:.4f} {match}  "
                      f"hist={dict(sorted(hist.items()))}  ({elapsed:.0f}s)",
                      flush=True)
        print(flush=True)

        # Backfill 완료 후 즉시 1차 histogram 분석
        print(f"━━ Backfill 후 첫 borderline 분석 (3 seeds) ━━━", flush=True)
        report_histogram_analysis(hists, K_means_all)
        print(flush=True)

    # ─── NEW SEEDS ───
    if to_run:
        print(f"━━ NEW SEEDS ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
              flush=True)
        total_cells = len(to_run) * len(SIGMA_GRID)
        cell_idx = 0
        t_new_start = time.time()

        for seed in to_run:
            print(f"\n━━ seed {seed} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                  flush=True)
            t_seed = time.time()

            for sigma in SIGMA_GRID:
                cell_idx += 1
                t_cell = time.time()
                noise_kwargs = {} if sigma == 0.0 else {"phase_sigma": sigma}
                K_mean, hist, Ks = measure_with_hist(
                    N, D, noise_kwargs, trials, base_seed=seed)
                K_sd = statistics.stdev(Ks) if len(Ks) > 1 else 0

                append_K_mean(seed, sigma, K_mean)
                append_histogram(seed, sigma, hist)
                K_means_all.setdefault(seed, {})[sigma] = K_mean
                hists.setdefault(seed, {})[sigma] = hist

                elapsed = time.time() - t_cell
                new_elapsed = time.time() - t_new_start
                eta = (new_elapsed * (total_cells - cell_idx) / cell_idx
                       if cell_idx > 0 else 0)
                print(f"  σ={sigma:.3f}  K_mean={K_mean:.4f} (sd={K_sd:.2f})  "
                      f"cell={elapsed:>4.0f}s  ETA={eta:>5.0f}s "
                      f"({cell_idx:>3}/{total_cells})", flush=True)

            # Per-seed summary
            K_base = K_means_all[seed][0.0]
            plateau_K = K_means_all[seed][0.050]
            sr = (K_base - plateau_K) / K_base * 100 if K_base > 0 else 0
            seed_time = time.time() - t_seed
            mark = " +" if sr > 0 else (" -" if sr < 0 else "")
            print(f"\n  seed {seed} 완료 ({seed_time:.0f}s):  "
                  f"K_base={K_base:.4f}  K(σ=.05)={plateau_K:.4f}  "
                  f"SR={sr:+.3f}%{mark}", flush=True)

            # 실시간 결합 분석
            report_combined_analysis(K_means_all, f"(seed {seed} 후)")
            print(flush=True)

    # ─── FINAL ANALYSIS ───
    print(f"\n━━ FINAL FULL ANALYSIS ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
          flush=True)
    report_combined_analysis(K_means_all, "(최종)")
    report_histogram_analysis(hists, K_means_all)

    print(f"\n  파일 저장:")
    print(f"    K_means: {EXTEND_FILE}")
    print(f"    Hist:    {HIST_FILE}")
    print(f"  총 소요: {time.time() - t_global:.0f}s "
          f"({(time.time() - t_global) / 60:.1f} 분)")


if __name__ == "__main__":
    main(sys.argv)
