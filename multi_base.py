"""
다중 base를 이용한 Carmichael λ(N) 점근 회수와 위수 추정 회수율 향상.

핵심 아이디어
-------------
모든 a ∈ (Z/N)* 의 위수 r_a 는 Carmichael 함수 λ(N) 을 나눈다.
여러 base 의 위수 lcm L = lcm(r_{a_1}, r_{a_2}, ...) 는 λ(N) 에 점근적으로 수렴.

L 이 일단 (Z/N)* 의 exponent 가 되면 (random b 에 대해 b^L ≡ 1 mod N),
이후 새 base a 의 위수 r_a 는 L 의 약수이므로 *추가 측정 없이* 고전적으로 회수 가능:
    r_a = min{ d | L : a^d ≡ 1 mod N }

단일 측정에서도 후처리를 강화:
    - k/Q 의 연분수 수렴값 *전체* 를 후보로 시도 (Knill-Mosca)
    - 누적 L 의 약수도 후보에 추가
    - 각 후보 d 에 대해 a^d ≡ 1 검증

또한 L 이 확정되면 회로 한 번 더 돌릴 필요 없이 Miller-Rabin 식으로 인수 추출.
"""

from __future__ import annotations
import math
import random
import sys
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from classical import FactorResult, is_prime, find_prime_power
from shor import simulate_period_finding


# ──────────────────────────────────────────────────────────────────
# 수학 도구
# ──────────────────────────────────────────────────────────────────

def divisors(n: int) -> list[int]:
    """n 의 양의 약수 (오름차순)."""
    if n <= 0:
        return []
    small: list[int] = []
    large: list[int] = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            if i != n // i:
                large.append(n // i)
        i += 1
    return small + large[::-1]


def prime_factors(n: int) -> list[int]:
    """n 의 서로 다른 소인수."""
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


def convergent_denominators(k: int, Q: int, max_denom: int) -> list[int]:
    """k/Q 의 연분수 수렴값 분모들 (max_denom 이하).

    Fraction(k, Q).limit_denominator(M) 은 best approximation 한 개만 주지만,
    실제 r 회수에는 모든 수렴값 분모를 시도하는 게 안정적이다.
    """
    a, b = k, Q
    h_prev, h_curr = 0, 1
    k_prev, k_curr = 1, 0
    denoms: list[int] = []
    while b != 0:
        q = a // b
        a, b = b, a - q * b
        h_prev, h_curr = h_curr, q * h_curr + h_prev
        k_prev, k_curr = k_curr, q * k_curr + k_prev
        if 0 < k_curr <= max_denom:
            denoms.append(k_curr)
        elif k_curr > max_denom:
            break
    return denoms


def minimize_order(a: int, N: int, candidate: int) -> int:
    """candidate 가 a 의 mod N 에서 exponent (a^candidate ≡ 1) 일 때
    실제 위수로 축소. 작은 소인수로 나눠가며 검증."""
    if candidate <= 0 or pow(a, candidate, N) != 1:
        return 0
    r = candidate
    for p in prime_factors(r):
        while r % p == 0 and pow(a, r // p, N) == 1:
            r //= p
    return r


# ──────────────────────────────────────────────────────────────────
# 다중 base 누적 상태
# ──────────────────────────────────────────────────────────────────

@dataclass
class MultiBaseState:
    """여러 base 의 위수를 누적해 λ(N) 후보 L 을 유지."""
    L: int = 1                                  # 지금까지 본 r_a 들의 lcm. λ(N) 의 약수.
    orders: dict[int, int] = field(default_factory=dict)  # a → r_a
    measurements: int = 0                       # 실제 양자 측정 호출 횟수
    classical_recoveries: int = 0               # 측정 없이 L 의 약수로 회수한 횟수

    def update(self, a: int, r_a: int) -> None:
        if r_a > 0:
            self.orders[a] = r_a
            self.L = math.lcm(self.L, r_a)

    def is_exponent_for(self, a: int, N: int) -> bool:
        """현재 L 이 a 에 대한 exponent 인가."""
        return self.L > 1 and pow(a, self.L, N) == 1


# ──────────────────────────────────────────────────────────────────
# 측정 후처리: 수렴값 + 누적 L 활용
# ──────────────────────────────────────────────────────────────────

def order_from_measurement(
    a: int, N: int, k: int, Q: int, state: MultiBaseState,
) -> int:
    """한 측정 k 에서 a 의 위수 후보 회수.

    후보 풀: k/Q 의 모든 연분수 수렴값 분모 + 누적 L 의 약수.
    각 후보 d 에 대해 a^d ≡ 1 검증. 최소 유효 d 를 반환 후 최소화.
    """
    candidates: set[int] = set(convergent_denominators(k, Q, N - 1))
    if state.L > 1:
        candidates.update(divisors(state.L))

    valid = [d for d in candidates if d > 0 and pow(a, d, N) == 1]
    if not valid:
        return 0
    return minimize_order(a, N, min(valid))


# ──────────────────────────────────────────────────────────────────
# 다중 base 위수 회수
# ──────────────────────────────────────────────────────────────────

def quantum_order_multi(
    a: int, N: int, state: MultiBaseState,
    shots: int = 4,
    rng: np.random.Generator | None = None,
) -> int:
    """누적 L 활용한 위수 회수.

    빠른 길: L 이 이미 a 의 exponent 면 측정 없이 즉시 r_a 회수.
    느린 길: 양자 측정 후 후처리. 매 측정 결과를 누적 L 로 보강.
    """
    rng = rng or np.random.default_rng()

    if state.is_exponent_for(a, N):
        state.classical_recoveries += 1
        return minimize_order(a, N, state.L)

    r_guess = 0
    for _ in range(shots):
        m = simulate_period_finding(a, N, rng=rng)
        state.measurements += 1
        d = order_from_measurement(a, N, m.k, m.Q, state)
        if d == 0:
            continue
        r_guess = d if r_guess == 0 else math.lcm(r_guess, d)
        if pow(a, r_guess, N) == 1:
            r_guess = minimize_order(a, N, r_guess)
            return r_guess
    return r_guess


# ──────────────────────────────────────────────────────────────────
# exponent → 인수분해 (Miller-Rabin 식)
# ──────────────────────────────────────────────────────────────────

def factor_from_exponent(
    N: int, L: int, rng: random.Random, max_attempts: int = 20,
) -> Optional[FactorResult]:
    """exponent L (모든 a 에 대해 a^L ≡ 1 mod N) 에서 N 의 인수 추출.

    L = 2^t · m (m 홀수) 로 분해.
    임의 a 에 대해 시퀀스 a^m, a^(2m), ..., a^(2^t m) = 1 에서
    처음 1 이 되는 직전 값이 1 의 자명하지 않은 제곱근 (±1 아님) 이면
    그 값으로 gcd(x±1, N) 가 인수.
    """
    if L <= 0:
        return None
    t, m = 0, L
    while m % 2 == 0:
        m //= 2
        t += 1

    for attempt in range(1, max_attempts + 1):
        a = rng.randrange(2, N)
        g = math.gcd(a, N)
        if g > 1:
            return FactorResult(g, N // g, "exponent_gcd", a=a, attempts=attempt)
        if pow(a, L, N) != 1:
            continue  # L 이 a 에 대해 exponent 가 아님

        x = pow(a, m, N)
        if x == 1:
            continue
        for _ in range(t):
            if x == N - 1:
                break
            y = pow(x, 2, N)
            if y == 1:
                for cand in (math.gcd(x - 1, N), math.gcd(x + 1, N)):
                    if 1 < cand < N:
                        return FactorResult(cand, N // cand,
                                            "exponent_factor", a=a, r=L,
                                            attempts=attempt)
                break
            x = y
    return None


# ──────────────────────────────────────────────────────────────────
# 전체 파이프라인
# ──────────────────────────────────────────────────────────────────

def shor_quantum_multi(
    N: int, max_bases: int = 10, shots_per_base: int = 3,
    seed: Optional[int] = None,
) -> tuple[Optional[FactorResult], MultiBaseState]:
    """다중 base 누적 + exponent 활용한 쇼어 인수분해.

    반환: (result, state) — state.measurements 로 양자 측정 횟수 확인,
    state.classical_recoveries 로 L 활용 회수 횟수 확인.
    """
    state = MultiBaseState()

    if N < 2:
        return None, state
    if N % 2 == 0:
        return FactorResult(2, N // 2, "even"), state
    if is_prime(N):
        return None, state
    mp = find_prime_power(N)
    if mp is not None:
        return FactorResult(mp, N // mp, "prime_power"), state

    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    for attempt in range(1, max_bases + 1):
        a = rng_py.randrange(2, N)
        g = math.gcd(a, N)
        if g > 1:
            return FactorResult(g, N // g, "gcd_shortcut", a=a, attempts=attempt), state

        r = quantum_order_multi(a, N, state, shots=shots_per_base, rng=rng_np)
        if r == 0:
            continue
        state.update(a, r)

        # 1) 표준 쇼어 환원
        if r % 2 == 0:
            x = pow(a, r // 2, N)
            if x != N - 1:
                for cand in (math.gcd(x - 1, N), math.gcd(x + 1, N)):
                    if 1 < cand < N:
                        return FactorResult(cand, N // cand,
                                            "period_multi", a=a, r=r,
                                            attempts=attempt), state

        # 2) 누적 L 활용한 Miller-Rabin 식 추가 시도
        result = factor_from_exponent(N, state.L, rng_py, max_attempts=5)
        if result is not None:
            result.attempts = attempt
            return result, state

    return None, state


# ──────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 33
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    print(f"=== 다중 base 쇼어 시뮬레이션: N={N} (seed={seed}) ===")
    result, state = shor_quantum_multi(N, seed=seed)
    print(f"결과: {result}")
    print(f"양자 측정: {state.measurements}, 고전 회수: {state.classical_recoveries}")
    print(f"누적 L: {state.L}, 본 base 수: {len(state.orders)}")
    print(f"본 위수들: {state.orders}")
