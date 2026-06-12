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

def lattice_close(samples: list[RegevSample], threshold: float) -> bool:
    """주어진 sample 부분집합이 dual lattice 에 가까운지 판정.

    RV §6 의 직관: uncorrupted 샘플 사이의 small integer combination 이 lattice
    근처에 있어야 함. corrupted 가 섞이면 멀어짐.

    구현: 샘플 행렬의 첫 행과 norm 비교. (정확한 RV 격자 구성은 TODO).
    """
    if not samples:
        return False
    Q = samples[0].Q
    d = len(samples[0].k_vec)

    # 매우 단순화: 샘플 사이의 차이 벡터의 (mod Q) norm 으로 근접도 판단
    if len(samples) < 2:
        return True
    diffs = []
    for i in range(1, len(samples)):
        diff = [
            (samples[i].k_vec[j] - samples[0].k_vec[j]) % Q
            for j in range(d)
        ]
        # signed to [-Q/2, Q/2)
        diff = [x - Q if x > Q // 2 else x for x in diff]
        diffs.append(diff)

    # sympy LLL on diff matrix
    if not diffs:
        return True
    try:
        M = Matrix(diffs)
        reduced = M.lll()
        # 짧은 벡터의 norm
        shortest = min(
            sum(int(reduced[i, j]) ** 2 for j in range(reduced.cols)) ** 0.5
            for i in range(reduced.rows)
        )
        return shortest < threshold
    except Exception:
        return True  # LLL 실패시 보수적으로 통과


def filter_uncorrupted(
    samples: list[RegevSample], target_size: int, threshold: float = 1e6,
) -> list[RegevSample]:
    """RV Algorithm 6.1: corrupted 샘플 제거.

    각 후보 sample 에 대해, 기존 B 와 함께 lattice_close 가 True 면 B 에 추가.
    |B| = target_size 가 될 때까지 반복.
    """
    B: list[RegevSample] = []
    for s in samples:
        if lattice_close(B + [s], threshold):
            B.append(s)
            if len(B) >= target_size:
                break
    return B


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


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "demo":
        demo()
    else:
        print("Usage: python -m experiments.rv_filter_lll demo")
        print("\n본 모듈은 STUB. 후속 작업 (todo):")
        print("  1. RV Algorithm 6.1 의 정확한 격자 구성 (Section 6, Algorithm 6.1)")
        print("  2. Regev LLL 후처리의 numpy/sympy 구현 (Regev 2023 §3.2)")
        print("  3. (C) 좌표별과의 head-to-head 비교 실험")
        print("  4. well-spread 가정의 numerical verification")
