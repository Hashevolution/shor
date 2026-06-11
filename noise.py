"""
노이즈가 있는 양자 측정의 시뮬레이션 (6가지 모델).

모델 분류
---------
A) 측정 결과만 영향 (k 자체에 작용)
   - depolarizing: k 가 확률 p 로 균등 무작위
   - readout_flip: k 의 각 비트를 확률 p 로 뒤집음 (XOR 노이즈)
   - bias_zero: 확률 p 로 k = 0 강제 (adversarial readout)

B) iQFT 입력 amp 에 영향 (peak shape 왜곡)
   - phase_sigma: amp 에 N(0, σ²) 가우시안 위상
   - amplitude_damp: amp[x] *= exp(-γx) (T1 모델, |0⟩ 으로 감쇠 편향)

C) modular exponentiation 자체에 영향 (구조 부분 파괴)
   - modexp_error: f(x) = a^x mod N 의 확률 q 로 무작위 값

§7.8 의 (C)-determinism 정리는 A, B 군 모두에 보편적으로 적용되어야 함.
C 군 (구조 노이즈) 도 측정 분포만 통해 L 회수에 영향 → 정리 적용 범위 안.
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
    readout_flip: float = 0.0,
    bias_zero: float = 0.0,
    phase_sigma: float = 0.0,
    amplitude_damp: float = 0.0,
    modexp_error: float = 0.0,
) -> PeriodMeasurement:
    """6가지 노이즈 모델 (개별 또는 조합).

    depolarizing ∈ [0, 1]: k 가 확률 p 로 균등 무작위.
    readout_flip ∈ [0, 1]: k 의 각 비트 위치마다 확률 p 로 flip.
    bias_zero ∈ [0, 1]: 확률 p 로 k = 0 강제.
    phase_sigma ≥ 0: iQFT 직전 amp 에 N(0, σ²) 위상 (라디안).
    amplitude_damp ∈ [0, 1]: amp[x] *= exp(-γ x), γ = amplitude_damp.
    modexp_error ∈ [0, 1]: f(x) 의 확률 q 로 [0, N) 균등 무작위 값.
    """
    if math.gcd(a, N) != 1:
        raise ValueError("gcd(a,N) must be 1")
    rng = rng or np.random.default_rng()
    t = t or _counting_qubits(N)
    Q = 1 << t

    # (A1) Depolarizing: k 직접 교체 (최우선)
    if depolarizing > 0 and rng.random() < depolarizing:
        k = int(rng.integers(0, Q))
        return _wrap(k, Q, N, y0=-1)

    # (A2) Bias zero
    if bias_zero > 0 and rng.random() < bias_zero:
        return _wrap(0, Q, N, y0=-1)

    # B, C 모델은 직접 시뮬 필요
    need_full_sim = (
        modexp_error > 0 or phase_sigma > 0 or amplitude_damp > 0
    )

    if need_full_sim:
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

        idx = int(rng.integers(0, Q))
        y0 = int(vals[idx])
        xs = np.flatnonzero(vals == y0)

        amps = np.zeros(Q, dtype=np.complex128)
        amps[xs] = 1.0 / math.sqrt(len(xs))

        # (B1) phase decoherence
        if phase_sigma > 0:
            phases = rng.normal(0.0, phase_sigma, size=Q)
            amps = amps * np.exp(1j * phases)

        # (B2) amplitude damping
        if amplitude_damp > 0:
            gamma = amplitude_damp
            decay = np.exp(-gamma * np.arange(Q))
            amps = amps * decay

        qft = np.fft.fft(amps) / math.sqrt(Q)
        probs = np.abs(qft) ** 2
        s = probs.sum()
        probs = probs / s if s > 0 else np.ones(Q) / Q
        k = int(rng.choice(Q, p=probs))
    else:
        # 노이즈 없음 → 원본 동작
        m_orig = simulate_period_finding(a, N, t=t, rng=rng)
        k, y0 = m_orig.k, m_orig.y0

    # (A3) Readout bit flip (모든 모델 위에 추가 적용 가능)
    if readout_flip > 0:
        for bit in range(t):
            if rng.random() < readout_flip:
                k ^= (1 << bit)

    return _wrap(k, Q, N, y0=y0)


def _wrap(k: int, Q: int, N: int, y0: int = -1) -> PeriodMeasurement:
    """k 로부터 PeriodMeasurement 생성."""
    frac = Fraction(k, Q).limit_denominator(N - 1)
    return PeriodMeasurement(
        k=k, Q=Q, y0=y0, fraction=frac,
        period_candidate=frac.denominator,
    )
