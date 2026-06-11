"""
노이즈가 있는 양자 측정의 시뮬레이션.

가장 단순한 모델: **depolarizing 채널**
    P_noisy(k) = (1 - p) · P_true(k) + p · (1/Q)

확률 p 로 측정 결과가 균등 무작위 k 가 되고, 그렇지 않으면 표준 측정 분포.
실제 하드웨어의 게이트 오류·디코히어런스를 거시적으로 근사한 모델이다.
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
) -> PeriodMeasurement:
    """depolarizing 노이즈가 있는 양자 측정.

    depolarizing ∈ [0, 1]:
        0 → 노이즈 없음 (원본 simulate_period_finding 와 동일).
        1 → 완전한 노이즈 (k 균등 무작위).
    """
    if math.gcd(a, N) != 1:
        raise ValueError("gcd(a,N) must be 1")
    rng = rng or np.random.default_rng()
    t = t or _counting_qubits(N)
    Q = 1 << t

    if depolarizing > 0 and rng.random() < depolarizing:
        # 균등 무작위 k. 다른 필드는 후처리에서 의미 없도록 채움.
        k = int(rng.integers(0, Q))
        frac = Fraction(k, Q).limit_denominator(N - 1)
        return PeriodMeasurement(
            k=k, Q=Q, y0=-1, fraction=frac,
            period_candidate=frac.denominator,
        )

    return simulate_period_finding(a, N, t=t, rng=rng)
