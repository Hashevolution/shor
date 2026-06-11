"""
쇼어 알고리즘의 양자 부분 (주기 찾기)을 numpy 상태벡터로 시뮬레이션.

설계 노트
---------
이론대로는 (계산 t큐비트) ⊗ (작업 n큐비트) 합동 상태벡터를 만들어 ModExp →
inverse QFT를 가해야 한다. 이는 차원 Q·2^n = 2^(t+n)의 벡터를 요구해
N ≳ 50 정도에서도 무거워진다.

여기서는 "**작업 레지스터를 먼저 측정**"한다. 측정 결과 y₀가 무엇이든,
계산 레지스터는 {x : a^x ≡ y₀ mod N}에 대한 균등 중첩으로 붕괴한다.
이 부분상태에 역 QFT를 적용한 측정 분포는 전체 회로의 한계 분포와 정확히 동일하다.

역 QFT는 정규화된 DFT이므로 numpy.fft.fft / sqrt(Q)로 처리한다.
"""

from __future__ import annotations
import math
import random
import sys
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from classical import shor_reduce, FactorResult


# ────────────────────────────────────────────────────────────────────────────
# 양자 주기 찾기 시뮬레이션
# ────────────────────────────────────────────────────────────────────────────

@dataclass
class PeriodMeasurement:
    k: int          # 계산 레지스터 측정값
    Q: int          # 2^t (계산 레지스터 크기)
    y0: int         # 작업 레지스터 측정값
    fraction: Fraction
    period_candidate: int


def _counting_qubits(N: int) -> int:
    """t = 2·⌈log₂ N⌉. Q ≥ N²이 되어 연분수 회수가 보장된다."""
    n = max(1, (N - 1).bit_length())
    return 2 * n


def simulate_period_finding(
    a: int,
    N: int,
    t: int | None = None,
    rng: np.random.Generator | None = None,
) -> PeriodMeasurement:
    """주기 찾기 회로를 한 번 실행한 효과를 시뮬레이션.

    Parameters
    ----------
    a : 곱셈 위수를 찾고자 하는 정수, gcd(a,N)=1.
    N : 합성수.
    t : 계산 레지스터 큐비트 수. 기본 2·⌈log₂ N⌉.
    rng : numpy Generator (시드 고정용).
    """
    if math.gcd(a, N) != 1:
        raise ValueError("gcd(a,N) must be 1")
    rng = rng or np.random.default_rng()
    t = t or _counting_qubits(N)
    Q = 1 << t

    # 1) f(x) = a^x mod N (x = 0..Q-1)
    vals = np.empty(Q, dtype=np.int64)
    cur = 1
    for x in range(Q):
        vals[x] = cur
        cur = (cur * a) % N

    # 2) 작업 레지스터 측정. 결합 상태가 균등 분포이므로,
    #    y0를 vals에서 균등 무작위 추출하면 동등한 분포.
    idx = int(rng.integers(0, Q))
    y0 = int(vals[idx])
    xs = np.flatnonzero(vals == y0)  # 등차수열 x0, x0+r, x0+2r, ...

    # 3) 계산 레지스터의 부분상태 진폭 (균등 중첩)
    amps = np.zeros(Q, dtype=np.complex128)
    amps[xs] = 1.0 / math.sqrt(len(xs))

    # 4) 역 QFT: 컨벤션상 numpy.fft.fft가 sum_x amps[x] e^(-2πixk/Q)
    qft = np.fft.fft(amps) / math.sqrt(Q)
    probs = np.abs(qft) ** 2
    probs /= probs.sum()  # 수치오차 보정

    # 5) k 표본추출
    k = int(rng.choice(Q, p=probs))

    # 6) 연분수로 k/Q ≈ j/r에서 r 후보 추출
    frac = Fraction(k, Q).limit_denominator(N - 1)
    r_candidate = frac.denominator
    return PeriodMeasurement(k=k, Q=Q, y0=y0, fraction=frac,
                             period_candidate=r_candidate)


def quantum_order(
    a: int,
    N: int,
    shots: int = 8,
    t: int | None = None,
    rng: np.random.Generator | None = None,
) -> int:
    """양자 주기 찾기를 여러 번 측정해 a의 N에 대한 위수를 회수.

    각 측정은 r의 약수만 줄 수도 있으므로 (j와 r의 공약수가 있을 때),
    후보들의 lcm을 취하고 a^r ≡ 1 mod N으로 검증한다.
    """
    rng = rng or np.random.default_rng()
    r_guess = 1
    for _ in range(shots):
        m = simulate_period_finding(a, N, t=t, rng=rng)
        c = m.period_candidate
        if c == 0:
            continue
        r_guess = math.lcm(r_guess, c)
        if r_guess > 0 and pow(a, r_guess, N) == 1:
            # 진짜 위수의 배수일 수도 있으니, 가장 작은 약수로 줄여본다.
            r = r_guess
            for p in _small_primes(r):
                while r % p == 0 and pow(a, r // p, N) == 1:
                    r //= p
            return r
    return 0  # 실패


def _small_primes(n: int) -> list[int]:
    """n의 소인수들 (작은 수에서만 사용)."""
    primes: list[int] = []
    x = n
    d = 2
    while d * d <= x:
        if x % d == 0:
            primes.append(d)
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        primes.append(x)
    return primes


# ────────────────────────────────────────────────────────────────────────────
# 전체 Shor (양자 위수 사용)
# ────────────────────────────────────────────────────────────────────────────

def shor_quantum(N: int, max_attempts: int = 30, seed: int | None = None) -> FactorResult | None:
    """고전 환원 + 양자 주기 찾기."""
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    def order_fn(a: int, N: int) -> int:
        r = quantum_order(a, N, rng=rng_np)
        if r == 0:
            raise RuntimeError("주기 회수 실패")
        return r

    return shor_reduce(N, order_fn=order_fn, max_attempts=max_attempts, rng=rng_py)


# ────────────────────────────────────────────────────────────────────────────
# 분석 보조: 측정 분포 직접 계산
# ────────────────────────────────────────────────────────────────────────────

def measurement_distribution(a: int, N: int, t: int | None = None) -> np.ndarray:
    """y0에 대해 주변화한 k의 확률 분포 P(k). 시각화/분석용."""
    t = t or _counting_qubits(N)
    Q = 1 << t
    vals = np.array([pow(a, x, N) for x in range(Q)], dtype=np.int64)

    probs = np.zeros(Q, dtype=np.float64)
    for y0 in np.unique(vals):
        xs = np.flatnonzero(vals == y0)
        amps = np.zeros(Q, dtype=np.complex128)
        amps[xs] = 1.0 / math.sqrt(len(xs))
        qft = np.fft.fft(amps) / math.sqrt(Q)
        # y0가 측정될 확률 = len(xs)/Q
        probs += (len(xs) / Q) * np.abs(qft) ** 2
    return probs / probs.sum()


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    print(f"=== 쇼어 양자 시뮬레이션: N={N} (seed={seed}) ===")
    result = shor_quantum(N, seed=seed)
    print(result)
