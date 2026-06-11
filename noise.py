"""
노이즈가 있는 양자 측정의 시뮬레이션 (세 가지 모델).

1) **Depolarizing** — 측정 k 자체가 확률 p 로 균등 무작위.
2) **Phase decoherence** — iQFT 직전 계산 레지스터 amp 에 가우시안 위상 노이즈.
   peak 가 분산되어 (A)/(B) 의 회수율 저하.
3) **Modular exponentiation error** — f(x) = a^x mod N 의 일부 값이 잘못 계산됨.
   주기 구조가 부분적으로 파괴됨 → (C) 의 견고함을 진짜로 시험.

세 모델의 차이:
- (1), (2) 는 단일 측정의 k 만 영향.
- (3) 은 후처리에서 검증되는 a^d ≡ 1 조건의 *전제* 인 a^x mod N 자체에 영향.
  (C) 가 의존하는 구조 정보 (r_a | λ(N)) 가 깨질 가능성.
"""

from __future__ import annotations
import math
from fractions import Fraction

import numpy as np

from shor import simulate_period_finding, _counting_qubits, PeriodMeasurement


def simulate_period_finding_noisy(
    a: int,
    N: int,
    t: int | None = None,
    rng: np.random.Generator | None = None,
    depolarizing: float = 0.0,
    phase_sigma: float = 0.0,
    modexp_error: float = 0.0,
) -> PeriodMeasurement:
    """세 가지 노이즈 모델 (개별 또는 조합).

    depolarizing ∈ [0, 1]: k 가 확률 p 로 균등 무작위.
    phase_sigma ≥ 0: iQFT 전 amp 에 N(0, σ²) 가우시안 위상 (라디안). 0 = 노이즈 없음.
    modexp_error ∈ [0, 1]: f(x) 의 확률 q 로 [0, N) 균등 무작위 값으로 교체.
    """
    if math.gcd(a, N) != 1:
        raise ValueError("gcd(a,N) must be 1")
    rng = rng or np.random.default_rng()
    t = t or _counting_qubits(N)
    Q = 1 << t

    # (1) Depolarizing: k 직접 교체
    if depolarizing > 0 and rng.random() < depolarizing:
        k = int(rng.integers(0, Q))
        frac = Fraction(k, Q).limit_denominator(N - 1)
        return PeriodMeasurement(
            k=k, Q=Q, y0=-1, fraction=frac,
            period_candidate=frac.denominator,
        )

    # (3) ModExp 오류 또는 (2) 위상 노이즈가 있으면 직접 시뮬
    if modexp_error > 0 or phase_sigma > 0:
        # f(x) 직접 계산 (오류 주입 가능)
        vals = np.empty(Q, dtype=np.int64)
        cur = 1
        for x in range(Q):
            vals[x] = cur
            cur = (cur * a) % N
        if modexp_error > 0:
            mask = rng.random(Q) < modexp_error
            n_err = int(mask.sum())
            if n_err > 0:
                vals[mask] = rng.integers(0, N, size=n_err)

        # 작업 레지스터 측정
        idx = int(rng.integers(0, Q))
        y0 = int(vals[idx])
        xs = np.flatnonzero(vals == y0)

        # 부분상태 진폭
        amps = np.zeros(Q, dtype=np.complex128)
        amps[xs] = 1.0 / math.sqrt(len(xs))

        # (2) 위상 노이즈: 각 basis state 에 random phase
        if phase_sigma > 0:
            phases = rng.normal(0.0, phase_sigma, size=Q)
            amps = amps * np.exp(1j * phases)

        qft = np.fft.fft(amps) / math.sqrt(Q)
        probs = np.abs(qft) ** 2
        s = probs.sum()
        if s > 0:
            probs /= s
        else:
            probs = np.ones(Q) / Q

        k = int(rng.choice(Q, p=probs))
        frac = Fraction(k, Q).limit_denominator(N - 1)
        return PeriodMeasurement(
            k=k, Q=Q, y0=y0, fraction=frac,
            period_candidate=frac.denominator,
        )

    # 모든 노이즈 0 → 원본 동작
    return simulate_period_finding(a, N, t=t, rng=rng)
