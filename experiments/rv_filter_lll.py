"""
Ragavan-Vaikuntanathan (2023, arXiv:2310.00899) Algorithm 6.1 — filter-then-LLL
구현 stub (Phase 5 follow-up).

RV's noise-tolerance 접근:
1. m = α·d 개 noisy 샘플 수집 (일부 corrupted by independent noise distribution D).
2. **Filter**: 각 sample 부분집합에 대해 격자 환원으로 "이 sample 들이 격자 근처에 있는가"
   판정. uncorrupted = 격자 근처, corrupted = 격자 멀리.
3. uncorrupted 만 모아 표준 Regev LLL 후처리.

본 모듈은 RV 의 Algorithm 6.1 을 numpy + sympy LLL 로 단순화 구현. 후속 세션에서
정확화 + 우리 (C) 좌표별 과의 head-to-head 비교에 사용.

상태: **STUB** — 핵심 구조 작성됨, 정확한 RV 파라미터 매핑 + Regev LLL 후처리는 후속.

References:
    [RV23] Ragavan & Vaikuntanathan, arXiv:2310.00899, §6.
    [Reg23] Regev, arXiv:2308.06572, classical post-processing.
"""

from __future__ import annotations
import math
import random
from dataclasses import dataclass

import numpy as np
from sympy import Matrix, Rational


@dataclass
class RegevSample:
    """Regev 측정의 한 출력 — d-dim 벡터 ∈ {0, ..., Q-1}^d."""
    k_vec: list[int]      # 측정값
    Q: int                # 2^t
    is_corrupted: bool = False  # ground truth (시뮬용)


# ──────────────────────────────────────────────────────────────────
# RV's Algorithm 6.1 — filter-then-LLL stub
# ──────────────────────────────────────────────────────────────────

def build_rv_lattice(
    samples_E: list[RegevSample], S: int,
) -> Matrix:
    """RV §6 Algorithm 6.1 (c): 격자 Λ ⊆ Z^(d+|E|) 의 basis 행렬 H.

    H 의 (d+|E|) × (d+|E|) 형태:
        [ S·I_d     S·W   ]
        [   0     I_|E|   ]
    여기서 W ∈ Z^(d × |E|) 의 열 = 샘플 벡터 w_i = k_i (Q-스케일 안 한 정수).

    Note: RV 는 w_i ∈ R^d/Z^d (= k_i/Q) 를 사용. 정수 LLL 을 위해 Q 로 곱해 사용.
    이로 인해 S 의 의미가 RV 원본과 약간 다를 수 있음 (실용적 검증용).

    Returns: sympy Matrix, 각 *행* 이 basis 벡터 (sympy LLL convention).
    """
    if not samples_E:
        return Matrix([])
    d = len(samples_E[0].k_vec)
    e = len(samples_E)
    n = d + e

    # H 를 행 단위로 생성 (sympy LLL 은 행 = basis vector)
    rows = []
    # 첫 d 개 행: (S·e_i, 0)
    for i in range(d):
        row = [0] * n
        row[i] = S
        rows.append(row)
    # 다음 |E| 개 행: (S·w_j, e_j)
    for j, s in enumerate(samples_E):
        row = [0] * n
        for i in range(d):
            row[i] = S * s.k_vec[i]
        row[d + j] = 1
        rows.append(row)

    return Matrix(rows)


def rv_filter_round(
    samples_E: list[RegevSample], S: int,
) -> list[int]:
    """RV Algorithm 6.1 의 한 round: LLL → 짧은 벡터 → support 인덱스.

    반환: a_j ≠ 0 인 j 의 인덱스 리스트 (uncorrupted 로 판정된 E 내 인덱스).
    """
    H = build_rv_lattice(samples_E, S)
    if H.rows == 0:
        return []
    try:
        reduced = H.lll()
    except Exception:
        return []

    d = len(samples_E[0].k_vec)
    # 가장 짧은 벡터 (행) 찾기
    norms = []
    for i in range(reduced.rows):
        norm_sq = sum(int(reduced[i, j]) ** 2 for j in range(reduced.cols))
        norms.append((norm_sq, i))
    norms.sort()
    short_row = reduced.row(norms[0][1])

    # 짧은 벡터의 마지막 |E| 좌표 = a_E
    a_E = [int(short_row[0, d + j]) for j in range(len(samples_E))]
    # a_j ≠ 0 인 인덱스
    return [j for j, a in enumerate(a_E) if a != 0]


def filter_uncorrupted(
    samples: list[RegevSample], target_size: int, S: int = 100,
    max_iter: int = 30,
) -> list[RegevSample]:
    """RV Algorithm 6.1 — full filter loop.

    매 iteration: 남은 샘플에서 E 선택 → rv_filter_round → support 인덱스를 B 에 추가.
    |B| ≥ target_size 가 될 때까지 반복.

    S: 격자 scaling. RV 의 S = 2^(An/d) 와 다르게 단순 정수 (실용 검증용).
       너무 크면 LLL 느림, 너무 작으면 corrupted 와 uncorrupted 가 구분 안 됨.
    """
    B: list[RegevSample] = []
    remaining = list(samples)
    for _ in range(max_iter):
        if len(B) >= target_size or not remaining:
            break
        # E 선택: 남은 모든 샘플
        E = remaining
        support_idx = rv_filter_round(E, S)
        if not support_idx:
            # LLL 이 nontrivial short vector 못 찾음 — 보수적으로 모두 통과
            B.extend(E[:target_size - len(B)])
            break
        # support 의 샘플을 B 에 추가, 나머지는 다음 round 대기
        added = [E[j] for j in support_idx]
        B.extend(added)
        remaining = [E[j] for j in range(len(E)) if j not in support_idx]
    return B[:target_size] if len(B) > target_size else B


def regev_lll_postprocess_stub(
    samples: list[RegevSample], N: int,
) -> int | None:
    """표준 Regev LLL 후처리 (TODO: 정확한 구현).

    현재 stub: filter 통과한 샘플 수만 출력. 실제 Regev factoring 은 미구현.
    """
    if len(samples) < 2:
        return None
    # TODO: Regev 의 정확한 lattice 구성 + LLL → 인수
    return len(samples)  # placeholder


# ──────────────────────────────────────────────────────────────────
# 시뮬레이션 helper — Regev 측정 + corruption 주입
# ──────────────────────────────────────────────────────────────────

def simulate_regev_run(
    bases: list[int], N: int, Q: int, rng: np.random.Generator,
    corrupt_prob: float = 0.0, noise_kwargs: dict | None = None,
) -> RegevSample:
    """단순화 Regev run: d 개 독립 Shor-like 측정.

    corrupt_prob 확률로 통째로 uniform 으로 교체 (RV 의 "overwrite" 모델).
    """
    from shor import simulate_period_finding
    from noise import simulate_period_finding_noisy

    if rng.random() < corrupt_prob:
        # 통째로 corrupted
        k_vec = [int(rng.integers(0, Q)) for _ in bases]
        return RegevSample(k_vec=k_vec, Q=Q, is_corrupted=True)

    k_vec = []
    for a in bases:
        if noise_kwargs:
            m = simulate_period_finding_noisy(a, N, rng=rng, **noise_kwargs)
        else:
            m = simulate_period_finding(a, N, rng=rng)
        k_vec.append(m.k)
    return RegevSample(k_vec=k_vec, Q=Q, is_corrupted=False)


# ──────────────────────────────────────────────────────────────────
# 임시 데모
# ──────────────────────────────────────────────────────────────────

def demo(N: int = 437, d: int = 4, m: int = 16, corrupt_prob: float = 0.2,
         seed: int = 0):
    """Demo: m 샘플 수집 → RV filter → 살아남은 수 vs ground truth."""
    rng_np = np.random.default_rng(seed)
    rng_py = random.Random(seed)
    Q = 1 << (2 * max(1, (N - 1).bit_length()))

    samples = []
    for _ in range(m):
        bases = []
        while len(bases) < d:
            a = rng_py.randrange(2, N)
            if math.gcd(a, N) == 1:
                bases.append(a)
        s = simulate_regev_run(bases, N, Q, rng_np, corrupt_prob=corrupt_prob)
        samples.append(s)

    true_corrupted = sum(1 for s in samples if s.is_corrupted)
    print(f"N={N}, d={d}, m={m}, corrupt_prob={corrupt_prob}")
    print(f"  실제 corrupted: {true_corrupted}/{m}")

    filtered = filter_uncorrupted(samples, target_size=d)
    filtered_corrupted = sum(1 for s in filtered if s.is_corrupted)
    print(f"  filter 통과: {len(filtered)} (그 중 corrupted: {filtered_corrupted})")
    print(f"  precision: {(len(filtered) - filtered_corrupted) / max(1, len(filtered)):.1%}")


def compare_with_c(
    N: int = 437, d: int = 4, corrupt_probs: list[float] | None = None,
    n_trials: int = 50, seed: int = 0,
):
    """(C) 좌표별 vs RV filter 의 head-to-head (부분).

    각 corrupt_prob 에 대해:
    - 우리 (C) 좌표별: 노이즈 + corruption 하에서 L=λ(N) 도달 # runs.
    - RV filter alone: precision (corrupted 검출률). 전체 RV (filter+LLL) 의 인수
      산출은 Regev LLL 미구현으로 비교 불가 — partial comparison only.
    """
    from classical import classical_order
    from multi_base import (
        MultiBaseState, convergent_denominators, divisors, minimize_order,
    )

    corrupt_probs = corrupt_probs or [0.0, 0.1, 0.2, 0.3]

    def lambda_of(N):
        for p in range(2, int(N**0.5) + 1):
            if N % p == 0:
                return math.lcm(p - 1, N // p - 1)

    lam = lambda_of(N)
    Q = 1 << (2 * max(1, (N - 1).bit_length()))

    print(f"\n# (C) 좌표별 vs RV filter 부분 비교, N={N}, d={d}, {n_trials} trials")
    print(f"# 우리 (C): L=λ(N) 도달 평균 # runs")
    print(f"# RV filter: 4 samples 선택 후 precision (1.00 = 모두 uncorrupted)")
    print(f"  {'corrupt_p':>10}  {'(C) runs':>10}  {'RV precision':>14}")

    for cp in corrupt_probs:
        c_runs_total = 0
        rv_precision_sum = 0
        rv_trials_succ = 0
        for t in range(n_trials):
            rng_np = np.random.default_rng(seed + t * 1000)
            rng_py = random.Random(seed + t * 1000)

            # (C) 시뮬: 매 run = d 개 base 좌표별 측정 + (C)
            state = MultiBaseState()
            c_runs = 0
            for K in range(1, 51):
                c_runs = K
                bases = []
                while len(bases) < d:
                    a = rng_py.randrange(2, N)
                    if math.gcd(a, N) == 1:
                        bases.append(a)
                run = simulate_regev_run(bases, N, Q, rng_np, corrupt_prob=cp)
                # (C) 좌표별
                for a, k in zip(bases, run.k_vec):
                    cands = set(convergent_denominators(k, Q, N - 1))
                    if state.L > 1:
                        cands.update(divisors(state.L))
                    valid = [x for x in cands if x > 0 and pow(a, x, N) == 1]
                    if valid:
                        r = minimize_order(a, N, min(valid))
                        if r > 0 and r == classical_order(a, N):
                            state.update(a, r)
                if state.L == lam:
                    break
            c_runs_total += c_runs

            # RV filter 시뮬: 동일 양의 노이즈 샘플 수집 후 filter 정밀도
            rng_np2 = np.random.default_rng(seed + t * 1000 + 999)
            rng_py2 = random.Random(seed + t * 1000 + 999)
            samples = []
            for _ in range(16):  # m=16 고정
                bases = []
                while len(bases) < d:
                    a = rng_py2.randrange(2, N)
                    if math.gcd(a, N) == 1:
                        bases.append(a)
                s = simulate_regev_run(bases, N, Q, rng_np2, corrupt_prob=cp)
                samples.append(s)
            filtered = filter_uncorrupted(samples, target_size=d)
            if filtered:
                corrupt_in_filtered = sum(1 for s in filtered if s.is_corrupted)
                prec = 1 - corrupt_in_filtered / len(filtered)
                rv_precision_sum += prec
                rv_trials_succ += 1

        avg_c = c_runs_total / n_trials
        avg_rv = rv_precision_sum / max(1, rv_trials_succ)
        print(f"  {cp:>10.2f}  {avg_c:>10.2f}  {avg_rv:>14.2%}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "demo":
        demo()
    elif len(sys.argv) > 1 and sys.argv[1] == "compare":
        compare_with_c()
    else:
        print("Usage: python -m experiments.rv_filter_lll demo")
        print("\n본 모듈은 STUB. 후속 작업 (todo):")
        print("  1. RV Algorithm 6.1 의 정확한 격자 구성 (Section 6, Algorithm 6.1)")
        print("  2. Regev LLL 후처리의 numpy/sympy 구현 (Regev 2023 §3.2)")
        print("  3. (C) 좌표별과의 head-to-head 비교 실험")
        print("  4. well-spread 가정의 numerical verification")
