"""
analyze_histograms.py — sigma_scan_437_d4_histograms.txt 의 모든 seed 검증.

목적:
  Seeds 1-9 (또는 현재까지) 의 K-histogram 비교로 *각 seed* 의 K=1/K=2 boundary
  flip mechanism 직접 확인. seeds 1-3 의 backfill 분석을 *모든 seeds* 로 확장.

분석:
  1. 각 seed 의 σ=0 vs σ=0.050 histogram 비교
  2. Diff histogram (gain/loss per K bin)
  3. Dominant flip K_from → K_to 식별 + magnitude
  4. K_mean 변화와 *직접* 일치 검증
  5. 전체 패턴 통계 (몇 seed 가 K=2→K=1, 몇 seed 가 K=1→K=2, etc.)

실행:
  python -u -m experiments.analyze_histograms                     # σ=0.050 비교
  python -u -m experiments.analyze_histograms 0.025               # σ=0.025 비교
  python -u -m experiments.analyze_histograms 0.150               # σ=0.150 비교 (decline 구간)
"""

from __future__ import annotations
import collections
import statistics
import sys
from pathlib import Path


HIST_FILE = Path("experiments/sigma_scan_437_d4_histograms.txt")


def read_histograms():
    hists: dict[int, dict[float, dict[int, int]]] = {}
    if not HIST_FILE.exists():
        print(f"ERROR: {HIST_FILE} 없음.")
        return hists
    with open(HIST_FILE, encoding="utf-8", errors="replace") as f:
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


def histogram_K_mean(hist):
    total_K = sum(K * count for K, count in hist.items())
    total_n = sum(hist.values())
    return total_K / total_n if total_n > 0 else 0, total_n


def main(argv):
    target_sigma = float(argv[1]) if len(argv) > 1 else 0.050

    hists = read_histograms()
    if not hists:
        return

    # 분석 가능 seed
    eligible = sorted(
        s for s in hists
        if 0.0 in hists[s] and target_sigma in hists[s]
    )
    if not eligible:
        print(f"σ={target_sigma} 데이터 있는 seed 없음.")
        return

    print(f"=" * 80)
    print(f"K-histogram 분석: σ=0 vs σ={target_sigma}  (n_seeds = {len(eligible)})")
    print(f"=" * 80)

    pos_flips = []  # (seed, K_from, K_to, magnitude)
    neg_flips = []

    for seed in eligible:
        h0 = hists[seed][0.0]
        hσ = hists[seed][target_sigma]

        K_base, n0 = histogram_K_mean(h0)
        Kσ, nσ = histogram_K_mean(hσ)
        sr = (K_base - Kσ) / K_base * 100 if K_base > 0 else 0

        # Diff histogram
        K_range = sorted(set(h0) | set(hσ))
        diffs = [(K, h0.get(K, 0), hσ.get(K, 0), hσ.get(K, 0) - h0.get(K, 0))
                 for K in K_range]

        # Dominant flip
        gains = [(K, d) for K, _, _, d in diffs if d > 0]
        losses = [(K, -d) for K, _, _, d in diffs if d < 0]
        if gains and losses:
            K_to, gain_mag = max(gains, key=lambda x: x[1])
            K_from, loss_mag = max(losses, key=lambda x: x[1])
            mag = min(gain_mag, loss_mag)
            flip_str = f"K={K_from}→K={K_to}"

            if sr > 0:
                pos_flips.append((seed, K_from, K_to, mag))
            elif sr < 0:
                neg_flips.append((seed, K_from, K_to, mag))
        else:
            K_from = K_to = -1
            mag = 0
            flip_str = "no change"

        # Per-seed 출력
        direction = "+" if sr > 0 else ("-" if sr < 0 else "=")
        print(f"\n▶ seed {seed}  K_base={K_base:.4f} (n={n0})  "
              f"K(σ={target_sigma})={Kσ:.4f} (n={nσ})  "
              f"SR={sr:+.3f}%  {direction}  flip={flip_str} (mag={mag})")

        # Histogram table
        print(f"   {'K':>3}  {'σ=0':>5}  {'σ=' + str(target_sigma):>7}  {'Δ':>5}  hint")
        for K, c0, cσ, d in diffs:
            if c0 == 0 and cσ == 0:
                continue
            hint = ""
            if d > 0:
                hint = f"+{d} ↑gain"
            elif d < 0:
                hint = f"{d} ↓loss"
            print(f"   {K:>3}  {c0:>5}  {cσ:>7}  {d:>+5}  {hint}")

        # Sanity: total K_mean diff 의 분석
        total_diff_K = sum(K * d for K, _, _, d in diffs)
        expected_diff = (Kσ - K_base) * n0
        print(f"   sanity: histogram Δ K total = {total_diff_K:+}, "
              f"K_mean Δ × n = {expected_diff:+.2f}  "
              f"({'✓ match' if abs(total_diff_K - expected_diff) < 0.5 else 'MISMATCH'})")

    # ─── 전체 패턴 통계 ───
    print(f"\n{'='*80}")
    print(f"전체 패턴 통계")
    print(f"{'='*80}")

    print(f"\nPositive SR seeds ({len(pos_flips)}):")
    for seed, K_from, K_to, mag in pos_flips:
        print(f"   seed {seed:>2}: K={K_from}→K={K_to}  mag={mag}")

    print(f"\nNegative SR seeds ({len(neg_flips)}):")
    for seed, K_from, K_to, mag in neg_flips:
        print(f"   seed {seed:>2}: K={K_from}→K={K_to}  mag={mag}")

    if pos_flips:
        from_counts = collections.Counter(K_from for _, K_from, _, _ in pos_flips)
        to_counts = collections.Counter(K_to for _, _, K_to, _ in pos_flips)
        mag_total = sum(mag for _, _, _, mag in pos_flips)
        print(f"\nPositive seeds 종합:")
        print(f"   K_from 분포: {dict(sorted(from_counts.items()))}")
        print(f"   K_to 분포:   {dict(sorted(to_counts.items()))}")
        print(f"   평균 magnitude: {mag_total / len(pos_flips):.1f}")

    if neg_flips:
        from_counts = collections.Counter(K_from for _, K_from, _, _ in neg_flips)
        to_counts = collections.Counter(K_to for _, _, K_to, _ in neg_flips)
        mag_total = sum(mag for _, _, _, mag in neg_flips)
        print(f"\nNegative seeds 종합:")
        print(f"   K_from 분포: {dict(sorted(from_counts.items()))}")
        print(f"   K_to 분포:   {dict(sorted(to_counts.items()))}")
        print(f"   평균 magnitude: {mag_total / len(neg_flips):.1f}")

    # K=1/K=2 boundary universality
    print(f"\n─── K=1/K=2 boundary universality 검증 ───")
    k12_seeds = []
    other_seeds = []
    for seed, K_from, K_to, mag in pos_flips + neg_flips:
        if {K_from, K_to} == {1, 2}:
            k12_seeds.append(seed)
        else:
            other_seeds.append((seed, K_from, K_to))

    n_total = len(pos_flips) + len(neg_flips)
    print(f"   K=1/K=2 boundary 가 dominant: {len(k12_seeds)}/{n_total} seeds")
    print(f"     seeds: {k12_seeds}")
    if other_seeds:
        print(f"   다른 boundary dominant: {len(other_seeds)}")
        for seed, K_from, K_to in other_seeds:
            print(f"     seed {seed}: K={K_from}→K={K_to}")
    else:
        print(f"   ★ 모든 seeds 가 K=1/K=2 boundary flip — universal mechanism")


if __name__ == "__main__":
    main(sys.argv)
