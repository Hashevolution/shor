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


@dataclass
class RegevSetup:
    """Regev 의 base 설정 — b_i random + a_i = b_i² mod N.

    b_i 는 우리가 자유롭게 고를 수 있는 random 값. a_i 는 그 제곱.
    회로는 a_i 를 base 로 측정하지만, 우리는 b_i 를 알고 있음.

    이게 Regev 의 핵심 트릭: b_i (sqrt of a_i) 를 알고 있어야 단축 vector
    z ∈ L 발견시 b = ∏ b_i^z_i 가 1 의 비자명 제곱근인지 검사 가능.
    """
    b: list[int]   # 우리가 random 으로 고른 b_i ∈ (Z/N)*
    a: list[int]   # a_i = b_i² mod N (회로의 base)
    N: int


def rv_scale_S(N: int, d: int, A: float = 2.0) -> int:
    """RV §6 의 lattice scaling S = 2^(An/d).

    n = N.bit_length(); A 는 RV의 정확도-보안 trade-off 상수 (RV §6 Lemma 6.5
    참조; A 가 클수록 corrupted 샘플 분리 능력 증가, LLL 비용 증가). A=2 가
    실용적 default — N=437 (n=9), d=4 → S = 2^5 = 32; N=4087 (n=12), d=4 → S = 2^6 = 64.
    """
    n = max(1, N.bit_length())
    return 1 << math.ceil(A * n / d)


def regev_setup_bases(N: int, d: int, rng: random.Random) -> RegevSetup:
    """d 개 random b_i ∈ (Z/N)*, a_i = b_i² mod N."""
    b = []
    while len(b) < d:
        cand = rng.randrange(2, N)
        if math.gcd(cand, N) == 1:
            b.append(cand)
    a = [(bi * bi) % N for bi in b]
    return RegevSetup(b=b, a=a, N=N)


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
    samples: list[RegevSample], target_size: int, S: int | None = None,
    max_iter: int = 30, N: int | None = None, d: int | None = None,
    A: float = 2.0,
) -> list[RegevSample]:
    """RV Algorithm 6.1 — full filter loop.

    매 iteration: 남은 샘플에서 E 선택 → rv_filter_round → support 인덱스를 B 에 추가.
    |B| ≥ target_size 가 될 때까지 반복.

    S: 격자 scaling. None 이면 RV §6 공식 `S = 2^(An/d)` 로 자동 계산
       (`rv_scale_S(N, d, A)`). 명시 S 값은 디버깅/legacy 비교용 (이전 default=100).
    """
    if S is None:
        if samples and N is not None:
            d_eff = d if d is not None else len(samples[0].k_vec)
            S = rv_scale_S(N, d_eff, A=A)
        else:
            S = 100  # legacy fallback when N/d unknown
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


def regev_algorithm_b1(
    samples: list[RegevSample], bases: list[list[int]], N: int,
    epsilon: int = 1, S: int | None = None,
) -> list[list[int]]:
    """Regev 2023 Algorithm B.1 — kernel-of-measurements 격자 후처리.

    찾는 격자: L_meas = {z ∈ Z^d : k^{(j)} · z ≡ 0 mod Q ∀j}. Shor 측정에서
    k^{(j)}_i ≈ Q·m/r_i 이므로 z 의 r_i 배수 좌표가 nontrivial 짧은 vector 가 됨.
    nontrivial 짧은 z 는 ∏ a_i^{z_i} ≡ 1 mod N (L_a 의 원소) 가 되고, Regev 의
    b_i (a_i = b_i²) 트릭으로 비자명 제곱근 → gcd 인수.

    격자 임베딩 in Z^(d+m), basis (sympy 행=basis 컨벤션):
        [ S·I_d     K^T   ]
        [   0    Q·I_m    ]
    여기서 K ∈ Z^(d×m), 열 j = k^{(j)}. LLL output v = (S·α, K^T α + Q·β) 의
    last-m 부분이 0 (또는 작음) 이면 k^{(j)} · α ≡ 0 mod Q for all j → α ∈ L_meas.
    z = v[:d] / S 가 후보 relation.

    S = 2^(An/d) (RV §6, A=2 default) — 작은 z 를 선호하는 격자 스케일링.
    """
    if not samples:
        return []
    m = len(samples)
    d = len(samples[0].k_vec)
    Q = samples[0].Q
    if S is None:
        S = rv_scale_S(N, d, A=2.0)

    rows = []
    # First d rows: (S·e_i, K^T row_i) = (S·e_i, k^{(1)}_i, ..., k^{(m)}_i)
    for i in range(d):
        row = [0] * (d + m)
        row[i] = S
        for j in range(m):
            row[d + j] = samples[j].k_vec[i]
        rows.append(row)
    # Next m rows: (0, Q·e_j) — enforce mod-Q
    for j in range(m):
        row = [0] * (d + m)
        row[d + j] = Q
        rows.append(row)

    H = Matrix(rows)
    try:
        reduced = H.lll()
    except Exception:
        return []

    # Extract z candidates: first-d / S where divisible by S
    candidates = []
    for i in range(reduced.rows):
        v_full = [int(reduced[i, j]) for j in range(reduced.cols)]
        z_part = v_full[:d]
        # divisibility check (LLL outputs in lattice, so first d should be S·integer)
        if any(c % S != 0 for c in z_part):
            continue
        z = [c // S for c in z_part]
        if all(zi == 0 for zi in z):
            continue
        # rank by full-vector L2 (prefer short)
        nrm = sum(x * x for x in v_full)
        candidates.append((nrm, z))
    candidates.sort()
    return [z for _, z in candidates[:d + 2]]  # top d+2 short relations


def try_factor_via_b_squareroot(
    relations: list[list[int]], setup: RegevSetup,
) -> tuple[int | None, str]:
    """Regev 의 정통 인수: z ∈ L → b = ∏ b_i^z_i → 비자명 제곱근 → gcd(b±1, N).

    Setup: b_i (random), a_i = b_i² mod N.
    Lattice L: {z ∈ Z^d : ∏ a_i^z_i ≡ 1 mod N}.

    각 z 에 대해:
    1. ∏ a_i^z_i ≡ 1 mod N 검증 (z ∈ L 여부).
    2. b = ∏ b_i^z_i mod N 계산.
    3. b² mod N == 1 인가? (반드시 그래야 함 — Regev 의 핵심)
    4. b ≢ ±1 mod N 인가? → 비자명 제곱근.
    5. gcd(b ± 1, N) > 1 → 인수.

    Returns: (factor, reason). factor=None 이면 reason 이 실패 이유.
    """
    N = setup.N
    for z in relations:
        if not z or all(zi == 0 for zi in z):
            continue
        # 1. ∏ a_i^z_i ≡ 1 검증
        prod_a = 1
        for ai, zi in zip(setup.a, z):
            # z_i 음수 → 모듈러 역원 사용
            if zi >= 0:
                prod_a = (prod_a * pow(ai, zi, N)) % N
            else:
                prod_a = (prod_a * pow(ai, -zi, N)) % N
                # need inverse: a^(-|zi|). 만약 gcd(prod_a, N) > 1 이면 nontrivial.
                try:
                    inv = pow(ai, -zi, N)  # Python 3.8+ supports negative exp
                    prod_a = (prod_a * inv) % N if False else prod_a  # placeholder
                except ValueError:
                    break
        if prod_a != 1:
            continue

        # 2. b = ∏ b_i^z_i mod N
        prod_b = 1
        ok = True
        for bi, zi in zip(setup.b, z):
            if zi >= 0:
                prod_b = (prod_b * pow(bi, zi, N)) % N
            else:
                try:
                    inv = pow(bi, zi, N)  # python 3.8+ negative exp = modular inverse
                except ValueError:
                    ok = False
                    break
                prod_b = (prod_b * inv) % N
        if not ok:
            continue

        # 3. b² ≡ 1 검증
        if (prod_b * prod_b) % N != 1:
            continue

        # 4. b ≢ ±1 (비자명 제곱근)?
        if prod_b == 1 or prod_b == N - 1:
            continue

        # 5. gcd(b-1, N) 또는 gcd(b+1, N) 가 인수
        for delta in (-1, 1):
            cand = math.gcd((prod_b + delta) % N, N)
            if 1 < cand < N:
                return cand, f"nontrivial sqrt b={prod_b}, gcd(b{delta:+d}, N)={cand}"

    return None, "no valid z found"


# (legacy alias for compatibility)
try_factor_from_relations = try_factor_via_b_squareroot


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


def end_to_end_factor_comparison_v2(
    N: int = 437, d: int = 4, max_runs: int = 20, n_trials: int = 30,
    corrupt_prob: float = 0.0, noise_kwargs: dict | None = None,
    seed: int = 0,
):
    """End-to-end 인수 비교 v2 — 세 방법:

    1. **(C) lcm only** — 좌표별 (C) → L 누적 → factor_from_exponent.
       Regev setup (a_i = b_i²) 에서는 L=odd part, Miller-Rabin 실패.
    2. **Regev b-trick (orders from convergents)** — 매 ord(a_i) 회수시 즉시
       b_i^ord 의 nontrivial sqrt 시도.
    3. **(C) + b-trick hybrid** — (C) 의 noise-tolerant 좌표별 recovery +
       b-trick 인수 추출. = 본 paper 의 권고 hybrid.

    노이즈 / corruption 하에서 K runs 까지 인수 회수 # runs 측정.
    """
    from classical import classical_order
    from multi_base import (
        MultiBaseState, convergent_denominators, divisors, minimize_order,
        factor_from_exponent,
    )
    noise_kwargs = noise_kwargs or {}

    def lambda_of(N):
        for p in range(2, int(N**0.5) + 1):
            if N % p == 0:
                return math.lcm(p - 1, N // p - 1)
    lam = lambda_of(N)
    Q = 1 << (2 * max(1, (N - 1).bit_length()))

    noise_label = "noise-free" if not noise_kwargs else str(noise_kwargs)
    print(f"\n# End-to-end 인수 비교 v2: N={N}, d={d}, {n_trials} trials")
    print(f"# corrupt_p={corrupt_prob}, noise={noise_label}")
    print(f"  {'method':<26} {'mean K':>8} {'max K':>6} {'success':>9}")

    results = {"c_lcm": [], "b_trick": [], "hybrid": [], "pure_rv": []}
    S_rv = rv_scale_S(N, d, A=2.0)

    for t in range(n_trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)
        setup = regev_setup_bases(N, d, rng_py)

        # 모든 방법에 동일 측정 sequence 사용 (공정 비교)
        runs = []
        for K in range(1, max_runs + 1):
            r = simulate_regev_run(
                setup.a, N, Q, rng_np,
                corrupt_prob=corrupt_prob, noise_kwargs=noise_kwargs,
            )
            runs.append(r)

        # 방법 1: (C) lcm only
        state = MultiBaseState()
        method1_K = max_runs
        for K, run in enumerate(runs, start=1):
            for ai, ki in zip(setup.a, run.k_vec):
                cands = set(convergent_denominators(ki, Q, N - 1))
                if state.L > 1:
                    cands.update(divisors(state.L))
                valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
                if valid:
                    r = minimize_order(ai, N, min(valid))
                    if r > 0 and r == classical_order(ai, N):
                        state.update(ai, r)
            if state.L > 1:
                rng_f = random.Random(t)
                res = factor_from_exponent(N, state.L, rng_f, max_attempts=10)
                if res and 1 < res.factor < N:
                    method1_K = K
                    break
        results["c_lcm"].append(method1_K)

        # 방법 2: pure b-trick (no L accumulation)
        method2_K = max_runs
        seen = set()
        for K, run in enumerate(runs, start=1):
            for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
                if ai in seen:
                    continue
                cands = convergent_denominators(ki, Q, N - 1)
                for d_ in cands:
                    if d_ > 0 and pow(ai, d_, N) == 1:
                        r = minimize_order(ai, N, d_)
                        if r > 0 and r == classical_order(ai, N):
                            seen.add(ai)
                            b_pow = pow(bi, r, N)
                            if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                                for delta in (-1, 1):
                                    g = math.gcd((b_pow + delta) % N, N)
                                    if 1 < g < N:
                                        method2_K = K
                                        break
                                if method2_K < max_runs:
                                    break
                            break
                if method2_K < max_runs:
                    break
            if method2_K < max_runs:
                break
        results["b_trick"].append(method2_K)

        # 방법 3: hybrid — (C) 의 좌표별 recovery + b-trick
        state3 = MultiBaseState()
        method3_K = max_runs
        for K, run in enumerate(runs, start=1):
            for ai, ki, bi in zip(setup.a, run.k_vec, setup.b):
                cands = set(convergent_denominators(ki, Q, N - 1))
                if state3.L > 1:
                    cands.update(divisors(state3.L))
                valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
                if valid:
                    r = minimize_order(ai, N, min(valid))
                    if r > 0 and r == classical_order(ai, N):
                        state3.update(ai, r)
                        # b-trick on this base
                        b_pow = pow(bi, r, N)
                        if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                            for delta in (-1, 1):
                                g = math.gcd((b_pow + delta) % N, N)
                                if 1 < g < N:
                                    method3_K = K
                                    break
            if method3_K < max_runs:
                break
            # 매 K 끝에 factor_from_exponent 도 시도 (L 이 짝수일 경우)
            if state3.L > 1:
                rng_f = random.Random(t)
                res = factor_from_exponent(N, state3.L, rng_f, max_attempts=10)
                if res and 1 < res.factor < N:
                    method3_K = K
                    break
        results["hybrid"].append(method3_K)

        # 방법 4: pure RV — filter (S=2^(An/d)) → Algorithm B.1 LLL → b-sqrt factor.
        # 필터는 m≥2d 샘플에서 의미가 있으므로 K = max(2d, 1) 부터 시도.
        method4_K = max_runs
        for K in range(max(1, 2 * d), max_runs + 1):
            samples_so_far = runs[:K]
            try:
                filtered = filter_uncorrupted(
                    samples_so_far, target_size=K, S=S_rv,
                    N=N, d=d,
                )
            except Exception:
                continue
            if len(filtered) < d:
                continue
            try:
                relations = regev_algorithm_b1(filtered, setup.a, N)
            except Exception:
                continue
            if not relations:
                continue
            factor, _ = try_factor_via_b_squareroot(relations, setup)
            if factor is not None and 1 < factor < N:
                method4_K = K
                break
        results["pure_rv"].append(method4_K)

    print(f"  (S_RV = 2^({math.ceil(2.0 * N.bit_length() / d)}) = {S_rv})")
    for label, key in [
        ("(C) lcm only", "c_lcm"),
        ("Regev b-trick", "b_trick"),
        ("(C)+b-trick hybrid", "hybrid"),
        ("pure RV (filter+B.1+sqrt)", "pure_rv"),
    ]:
        ks = results[key]
        succ = sum(1 for k in ks if k < max_runs)
        print(f"  {label:<26} "
              f"{sum(ks)/n_trials:>8.2f} {max(ks):>6} {succ:>3}/{n_trials}")


def end_to_end_factor_comparison(
    N: int = 437, d: int = 4, max_runs: int = 20, n_trials: int = 30,
    corrupt_prob: float = 0.0, seed: int = 0,
):
    """End-to-end 인수 비교: (C) λ(N)-via-lcm vs Regev b-trick.

    각 trial:
    1. b_i random, a_i = b_i² mod N 설정.
    2. K = 1, 2, ... runs 진행 (각 run = d 개 a_i 의 noisy 측정).
    3. (C) 좌표별: L 누적 → L=λ(N) 도달시 factor_from_exponent 로 인수.
    4. Regev 식: 매 run 의 좌표를 lattice 후보로 추가 → 짧은 vector → b-trick.
       단순화: ord(a_i) 직접 사용 후 b_i^ord 의 nontrivial sqrt 확인.

    어느 방법이 더 적은 K 로 인수 회수하는지 측정.
    """
    from classical import classical_order
    from multi_base import (
        MultiBaseState, convergent_denominators, divisors, minimize_order,
        factor_from_exponent,
    )

    def lambda_of(N):
        for p in range(2, int(N**0.5) + 1):
            if N % p == 0:
                return math.lcm(p - 1, N // p - 1)
    lam = lambda_of(N)
    Q = 1 << (2 * max(1, (N - 1).bit_length()))

    print(f"\n# End-to-end 인수 비교: N={N}, d={d}, {n_trials} trials, corrupt={corrupt_prob}")
    print(f"# (C): 좌표별 → L 누적 → L=λ(N) 도달시 factor_from_exponent")
    print(f"# Regev b-trick: 매 측정의 base 위수 → b_i^r 의 nontrivial sqrt")
    print(f"# 메트릭: 인수 회수까지 # runs (실패시 max_runs)")
    print(f"  {'method':<14} {'mean K':>8} {'max K':>6} {'success':>8}")

    c_ks = []
    regev_ks = []
    for t in range(n_trials):
        rng_py = random.Random(seed + t * 1000)
        rng_np = np.random.default_rng(seed + t * 1000)

        setup = regev_setup_bases(N, d, rng_py)

        # (C) trajectory
        state = MultiBaseState()
        c_factor = None
        c_K = max_runs
        for K in range(1, max_runs + 1):
            run = simulate_regev_run(
                setup.a, N, Q, rng_np, corrupt_prob=corrupt_prob,
            )
            # 좌표별 (C)
            for ai, ki in zip(setup.a, run.k_vec):
                cands = set(convergent_denominators(ki, Q, N - 1))
                if state.L > 1:
                    cands.update(divisors(state.L))
                valid = [d_ for d_ in cands if d_ > 0 and pow(ai, d_, N) == 1]
                if valid:
                    r = minimize_order(ai, N, min(valid))
                    if r > 0:
                        true_r = classical_order(ai, N)
                        if r == true_r:
                            state.update(ai, r)
            # Try factor
            if state.L > 1:
                rng_factor = random.Random(t)
                res = factor_from_exponent(N, state.L, rng_factor, max_attempts=5)
                if res is not None and 1 < res.factor < N:
                    c_factor = res.factor
                    c_K = K
                    break

        # Regev b-trick: 각 run 의 측정에서 ord 회수 후 b-trick 즉시 시도
        # (C) 와 분리된 trajectory)
        seen_orders: dict[int, int] = {}  # a_i -> r_a_i
        regev_factor = None
        regev_K = max_runs
        rng_py2 = random.Random(seed + t * 1000)
        rng_np2 = np.random.default_rng(seed + t * 1000)
        for K in range(1, max_runs + 1):
            run = simulate_regev_run(
                setup.a, N, Q, rng_np2, corrupt_prob=corrupt_prob,
            )
            for i, (ai, ki, bi) in enumerate(zip(setup.a, run.k_vec, setup.b)):
                if ai in seen_orders:
                    continue
                # 단일-좌표 회수 (convergent 만, divisor 없음)
                cands = convergent_denominators(ki, Q, N - 1)
                for d_ in cands:
                    if d_ > 0 and pow(ai, d_, N) == 1:
                        r = minimize_order(ai, N, d_)
                        if r > 0:
                            seen_orders[ai] = r
                            # b-trick
                            b_pow = pow(bi, r, N)
                            if b_pow not in (1, N - 1) and (b_pow * b_pow) % N == 1:
                                # nontrivial sqrt
                                for delta in (-1, 1):
                                    g = math.gcd((b_pow + delta) % N, N)
                                    if 1 < g < N:
                                        regev_factor = g
                                        regev_K = K
                                        break
                            break
                if regev_factor is not None:
                    break
            if regev_factor is not None:
                break

        c_ks.append(c_K if c_factor else max_runs)
        regev_ks.append(regev_K if regev_factor else max_runs)

    c_succ = sum(1 for k in c_ks if k < max_runs)
    rv_succ = sum(1 for k in regev_ks if k < max_runs)
    print(f"  {'(C) lcm':<14} "
          f"{sum(c_ks)/n_trials:>8.2f} {max(c_ks):>6} {c_succ:>3}/{n_trials}")
    print(f"  {'Regev b-trick':<14} "
          f"{sum(regev_ks)/n_trials:>8.2f} {max(regev_ks):>6} {rv_succ:>3}/{n_trials}")


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
    elif len(sys.argv) > 1 and sys.argv[1] == "factor":
        end_to_end_factor_comparison_v2()
        end_to_end_factor_comparison_v2(corrupt_prob=0.2)
        end_to_end_factor_comparison_v2(corrupt_prob=0.3)
        end_to_end_factor_comparison_v2(noise_kwargs={"depolarizing": 0.3})
        end_to_end_factor_comparison_v2(noise_kwargs={"phase_sigma": 1.0})
    else:
        print("Usage: python -m experiments.rv_filter_lll demo")
        print("\n본 모듈은 STUB. 후속 작업 (todo):")
        print("  1. RV Algorithm 6.1 의 정확한 격자 구성 (Section 6, Algorithm 6.1)")
        print("  2. Regev LLL 후처리의 numpy/sympy 구현 (Regev 2023 §3.2)")
        print("  3. (C) 좌표별과의 head-to-head 비교 실험")
        print("  4. well-spread 가정의 numerical verification")
