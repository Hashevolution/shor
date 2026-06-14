"""
method_bakeoff.py — ε×method crossover 측정 (noise-adaptive 후처리 발명의 make-or-break).

핵심 질문:
  노이즈 σ(→ε = 1 - exp(-σ²)) 를 올릴 때, 인수까지의 run 수 E[K] 의
  최소값을 주는 후처리 method 가 *바뀌는* crossover 구간이 존재하는가?

  - crossover 有 (regret > 0)  → 어떤 단일 method 도 모든 ε 에서 이기지 못함
                                 → noise-keyed selector 가 synergy 를 가짐 (발명 정당화).
  - crossover 無 (regret ≈ 0)  → 한 method (Theorem 5 상 hybrid) 가 ε-지배
                                 → "hybrid 은 ε-robust" robustness 결과 (Theorem 7 후보).

  Theorem 6 (SR no-go) 와의 정합: 이건 noise-as-resource 가 아니라 robustness/
  efficiency 측정이다. 적응 규칙은 (만약 정당화되면) per-run 관측 불가한 ε 대신
  하드웨어 보고 noise 에 키잉되어야 한다 (Yang-Markidis IBM 데이터 연계).

비교 method (CF 는 (C) 의 부분집합이라 제외; 실질 선택지만):
  1. (C) lcm only        — 좌표별 (C) → L 누적 → factor_from_exponent.
  2. Regev b-trick       — base 위수 회수 즉시 b_i^r 의 nontrivial sqrt.
  3. (C)+b-trick hybrid  — 본 paper 권고 (Theorem 5). 위 둘의 superset.

  (RV filter-then-LLL 은 rv_filter_lll.py 에서 아직 end-to-end 인수 STUB 이라
   여기 미포함. 구현되면 4번째 method 로 추가 — selector 가 RV-fixed 보다
   이득인지 확인 필요.)

공정성: 매 trial 에서 세 method 가 *동일한* 측정 시퀀스를 공유한다
  (rv_filter_lll.end_to_end_factor_comparison_v2 의 설계를 그대로 따름).

Reproduction:
  python -u -m experiments.method_bakeoff
  python -u -m experiments.method_bakeoff --N 437 --trials 30 --d 4
  python -u -m experiments.method_bakeoff --N 1147 --trials 20   # 판정 재확인
"""
from __future__ import annotations

import argparse
import math
import random
import statistics
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import (
    MultiBaseState,
    convergent_denominators,
    divisors,
    factor_from_exponent,
    minimize_order,
)
from noise import simulate_period_finding_noisy
from shor import simulate_period_finding


def regev_bases(N: int, d: int, rng: random.Random) -> tuple[list[int], list[int]]:
    """d 개 random b_i ∈ (Z/N)*, a_i = b_i² mod N. (rv_filter_lll 와 동일.)"""
    b: list[int] = []
    while len(b) < d:
        cand = rng.randrange(2, N)
        if math.gcd(cand, N) == 1:
            b.append(cand)
    a = [(bi * bi) % N for bi in b]
    return b, a


def regev_run(bases: list[int], N: int, Q: int, rng: np.random.Generator,
              noise_kwargs: dict | None) -> list[int]:
    """한 Regev run = d 개 base 의 (noisy) 측정값 k_vec."""
    out = []
    for a in bases:
        if noise_kwargs:
            m = simulate_period_finding_noisy(a, N, rng=rng, **noise_kwargs)
        else:
            m = simulate_period_finding(a, N, rng=rng)
        out.append(m.k)
    return out

# σ sweep: ε = 1 - exp(-σ²). 저노이즈 ~ 거의-완전 dephasing 전구간.
#   σ:  0.0  0.2   0.4   0.6   0.8   1.0   1.5   2.0
#   ε:  0.0  .039  .148  .302  .473  .632  .895  .982
DEFAULT_SIGMAS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]

METHOD_LABELS = [
    ("c_lcm", "(C) lcm only"),
    ("b_trick", "Regev b-trick"),
    ("hybrid", "(C)+b-trick hybrid"),
]

RESULTS_FILE = Path("experiments/method_bakeoff_results.txt")


def measure_methods(
    N: int, d: int, max_runs: int, n_trials: int, phase_sigma: float, seed: int,
) -> dict[str, list[int]]:
    """세 method 의 trial 별 K (인수까지 run 수, 실패시 max_runs) 를 측정.

    rv_filter_lll.end_to_end_factor_comparison_v2 의 로직을 그대로 따르되
    출력 대신 K 리스트를 반환한다. 세 method 는 trial 마다 동일한 측정
    시퀀스(runs)를 공유한다.
    """
    noise_kwargs = {"phase_sigma": phase_sigma} if phase_sigma > 0 else {}
    Q = 1 << (2 * max(1, (N - 1).bit_length()))
    results: dict[str, list[int]] = {"c_lcm": [], "b_trick": [], "hybrid": []}

    for t in range(n_trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)
        b_vec, a_vec = regev_bases(N, d, rng_py)

        # 동일 측정 시퀀스 (공정 비교)
        runs = [
            regev_run(a_vec, N, Q, rng_np, noise_kwargs)
            for _ in range(max_runs)
        ]

        # 방법 1: (C) lcm only
        state = MultiBaseState()
        k1 = max_runs
        for K, k_vec in enumerate(runs, start=1):
            for ai, ki in zip(a_vec, k_vec):
                cands = set(convergent_denominators(ki, Q, N - 1))
                if state.L > 1:
                    cands.update(divisors(state.L))
                valid = [x for x in cands if x > 0 and pow(ai, x, N) == 1]
                if valid:
                    r = minimize_order(ai, N, min(valid))
                    if r > 0 and r == classical_order(ai, N):
                        state.update(ai, r)
            if state.L > 1:
                res = factor_from_exponent(N, state.L, random.Random(t), max_attempts=10)
                if res and 1 < res.factor < N:
                    k1 = K
                    break
        results["c_lcm"].append(k1)

        # 방법 2: pure b-trick (L 누적 없음)
        k2 = max_runs
        seen: set[int] = set()
        for K, k_vec in enumerate(runs, start=1):
            for ai, ki, bi in zip(a_vec, k_vec, b_vec):
                if ai in seen:
                    continue
                for d_ in convergent_denominators(ki, Q, N - 1):
                    if d_ > 0 and pow(ai, d_, N) == 1:
                        r = minimize_order(ai, N, d_)
                        if r > 0 and r == classical_order(ai, N):
                            seen.add(ai)
                            b_pow = pow(bi, r, N)
                            if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                                for delta in (-1, 1):
                                    g = math.gcd((b_pow + delta) % N, N)
                                    if 1 < g < N:
                                        k2 = K
                                        break
                            break
                if k2 < max_runs:
                    break
            if k2 < max_runs:
                break
        results["b_trick"].append(k2)

        # 방법 3: hybrid — 좌표별 (C) recovery + b-trick
        state3 = MultiBaseState()
        k3 = max_runs
        for K, k_vec in enumerate(runs, start=1):
            for ai, ki, bi in zip(a_vec, k_vec, b_vec):
                cands = set(convergent_denominators(ki, Q, N - 1))
                if state3.L > 1:
                    cands.update(divisors(state3.L))
                valid = [x for x in cands if x > 0 and pow(ai, x, N) == 1]
                if valid:
                    r = minimize_order(ai, N, min(valid))
                    if r > 0 and r == classical_order(ai, N):
                        state3.update(ai, r)
                        b_pow = pow(bi, r, N)
                        if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                            for delta in (-1, 1):
                                g = math.gcd((b_pow + delta) % N, N)
                                if 1 < g < N:
                                    k3 = K
                                    break
            if k3 < max_runs:
                break
            if state3.L > 1:
                res = factor_from_exponent(N, state3.L, random.Random(t), max_attempts=10)
                if res and 1 < res.factor < N:
                    k3 = K
                    break
        results["hybrid"].append(k3)

    return results


def summarize(ks: list[int], max_runs: int) -> tuple[float, float, float]:
    """(mean K, standard error, success rate) 반환."""
    n = len(ks)
    mean = statistics.mean(ks)
    se = (statistics.pstdev(ks) / math.sqrt(n)) if n > 1 else 0.0
    success = sum(1 for k in ks if k < max_runs) / n
    return mean, se, success


def main():
    ap = argparse.ArgumentParser(description="ε×method crossover bake-off")
    ap.add_argument("--N", type=int, default=437)
    ap.add_argument("--d", type=int, default=4)
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--max-runs", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--success-floor", type=float, default=0.8,
        help="best-method 판정에 포함할 최소 success rate",
    )
    args = ap.parse_args()

    t0 = time.time()
    lines: list[str] = []

    def emit(s: str = "") -> None:
        print(s)
        lines.append(s + "\n")

    emit(f"# ε×method bake-off  (N={args.N}, d={args.d}, {args.trials} trials, "
         f"max_runs={args.max_runs})")
    emit(f"# ε = 1 - exp(-σ²).  best = min mean-K among methods with "
         f"success ≥ {args.success_floor:.0%}")
    emit(f"# regret(σ) = mean_K[hybrid] - mean_K[best].  crossover ⇔ best "
         f"method changes with σ AND gap > 2·SE")
    emit()

    # per-σ 측정
    rows = []  # (sigma, eps, {method: (mean, se, success)})
    best_per_sigma = []  # (sigma, best_key 또는 None)
    for sigma in DEFAULT_SIGMAS:
        eps = 1.0 - math.exp(-(sigma ** 2))
        t_s = time.time()
        ks = measure_methods(
            args.N, args.d, args.max_runs, args.trials, sigma, args.seed,
        )
        stats = {k: summarize(ks[k], args.max_runs) for k, _ in METHOD_LABELS}
        dt = time.time() - t_s

        emit(f"## σ={sigma:.2f}  (ε={eps:.3f})   [{dt:.0f}s]")
        emit(f"  {'method':<22} {'mean K':>8} {'±SE':>7} {'success':>9}")
        for key, label in METHOD_LABELS:
            mean, se, succ = stats[key]
            emit(f"  {label:<22} {mean:>8.2f} {se:>7.2f} {succ:>8.0%}")

        # best 판정 (success_floor 이상만)
        eligible = [(stats[k][0], k) for k, _ in METHOD_LABELS
                    if stats[k][2] >= args.success_floor]
        if eligible:
            best_mean, best_key = min(eligible)
            hyb_mean = stats["hybrid"][0]
            # tie-break: hybrid 가 best 와 사실상 동률(within SE)이면 hybrid 를 best 로.
            # alphabetical min() 이 1.30=1.30 동률에서 b_trick 을 가짜 winner 로 뽑는 것 방지.
            hyb_se = stats["hybrid"][1]
            if (best_key != "hybrid"
                    and stats["hybrid"][2] >= args.success_floor
                    and hyb_mean <= best_mean + max(hyb_se, 1e-9)):
                best_key, best_mean = "hybrid", hyb_mean
            regret = hyb_mean - best_mean
            label = dict(METHOD_LABELS)[best_key]
            emit(f"  → best: {label}  (mean K={best_mean:.2f}); "
                 f"always-hybrid regret = {regret:+.2f}")
            best_per_sigma.append((sigma, best_key, best_mean, stats[best_key][1]))
        else:
            emit("  → best: (모든 method success < floor — all-fail regime)")
            best_per_sigma.append((sigma, None, None, None))
        emit()

        rows.append((sigma, eps, stats))
        RESULTS_FILE.write_text("".join(lines), encoding="utf-8")  # 증분 저장

    # ── verdict: hybrid 의 always-on regret 을 headline 으로 ──
    # 핵심 질문은 "best method 가 바뀌나"가 아니라 "항상 hybrid 를 써서 얼마나 손해보나".
    # max regret ≈ 0 이면 noise-keyed selector 는 불필요 (synergy 0) → robustness 결과.
    emit("## Verdict — always-hybrid regret across σ")
    regrets = []  # (sigma, regret, thr)
    for (sigma, eps, stats), (s, bkey, bmean, bse) in zip(rows, best_per_sigma):
        if bkey is None:
            continue
        hyb_mean, hyb_se, _ = stats["hybrid"]
        regret = hyb_mean - bmean
        thr = 2.0 * (hyb_se + (bse or 0.0))
        regrets.append((sigma, regret, thr))
    if not regrets:
        emit("  평가가능 σ 없음 (전 구간 all-fail).")
    else:
        max_regret = max(r for _, r, _ in regrets)
        s_at_max, _, thr_at_max = max(regrets, key=lambda t: t[1])
        emit(f"  max regret = {max_regret:+.2f} at σ={s_at_max:.2f} "
             f"(2·SE there = {thr_at_max:.2f})")
        if max_regret <= thr_at_max:
            emit("  → 모든 σ 에서 hybrid 의 손해가 noise floor 안쪽. "
                 "crossover 없음 → noise-keyed selector synergy 0.")
            emit("  → 양수 해석: 'hybrid 은 measured ε∈[0,0.98] 에서 ε-dominant' "
                 "(robustness 관측; selector 발명 불필요).")
        else:
            emit("  → 어떤 σ 에서 hybrid 가 best 에 유의하게 뒤짐 → crossover 후보. "
                 "noise-keyed selector 정당화 여지 (trials 늘려 재확인 권장).")
            emit("  주의: novelty 는 'selector 착상'이 아니라 'crossover 발견 + "
                 "규칙'에 있음. RV-filter/YM-predictor 와 구분되게 좁게 framing.")

    elapsed = time.time() - t0
    emit()
    emit(f"# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
