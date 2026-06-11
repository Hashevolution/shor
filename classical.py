"""
쇼어 알고리즘의 고전 부분 + 고전적 위수 계산 baseline.

양자 시뮬레이션의 정답을 검증하고, 환원 단계(짝수 처리, prime power 검출,
gcd 추출)를 재사용 가능한 형태로 모은다.
"""

from __future__ import annotations
import math
import random
import sys
from dataclasses import dataclass


def is_prime(n: int) -> bool:
    """결정적 Miller-Rabin (n < 3.3 * 10^24까지 정확)."""
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def integer_kth_root(n: int, k: int) -> int:
    """floor(n^(1/k))를 정확히 계산 (정수 이분탐색)."""
    if k == 1:
        return n
    lo, hi = 1, 1 << ((n.bit_length() + k - 1) // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def find_prime_power(n: int) -> int | None:
    """n = m^k (k ≥ 2, m ≥ 2)이면 m 반환, 아니면 None."""
    for k in range(2, n.bit_length() + 1):
        m = integer_kth_root(n, k)
        if m >= 2 and m ** k == n:
            return m
    return None


def classical_order(a: int, N: int, limit: int | None = None) -> int:
    """a의 mod N에서의 곱셈 위수 r (가장 작은 r > 0, a^r ≡ 1)."""
    if math.gcd(a, N) != 1:
        raise ValueError(f"gcd({a},{N})≠1이라 위수가 정의되지 않음")
    limit = limit or N
    x, r = a % N, 1
    while x != 1:
        x = (x * a) % N
        r += 1
        if r > limit:
            raise RuntimeError(f"limit={limit}까지 주기를 못 찾음")
    return r


@dataclass
class FactorResult:
    factor: int
    cofactor: int
    method: str           # "even" | "prime_power" | "gcd_shortcut" | "period"
    a: int | None = None  # 사용한 base
    r: int | None = None  # 발견한 위수
    attempts: int = 0

    def __str__(self) -> str:
        s = f"{self.factor} × {self.cofactor} (method={self.method}"
        if self.a is not None:
            s += f", a={self.a}"
        if self.r is not None:
            s += f", r={self.r}"
        return s + f", attempts={self.attempts})"


def shor_reduce(N: int, order_fn=classical_order, max_attempts: int = 30,
                rng: random.Random | None = None) -> FactorResult | None:
    """쇼어의 고전 환원 셸. 위수 계산은 order_fn에 위임 (양자/고전 교체 가능).

    order_fn(a, N) → r (a^r ≡ 1 mod N)
    """
    if N < 2:
        return None
    if N % 2 == 0:
        return FactorResult(2, N // 2, "even")
    if is_prime(N):
        return None
    m = find_prime_power(N)
    if m is not None:
        return FactorResult(m, N // m, "prime_power")

    rng = rng or random.Random()
    for attempt in range(1, max_attempts + 1):
        a = rng.randrange(2, N)
        g = math.gcd(a, N)
        if g > 1:
            return FactorResult(g, N // g, "gcd_shortcut", a=a, attempts=attempt)

        try:
            r = order_fn(a, N)
        except Exception:
            continue
        if r == 0 or pow(a, r, N) != 1:
            continue
        if r % 2 != 0:
            continue

        x = pow(a, r // 2, N)
        if x == N - 1:
            continue

        p = math.gcd(x - 1, N)
        q = math.gcd(x + 1, N)
        for cand in (p, q):
            if 1 < cand < N:
                return FactorResult(cand, N // cand, "period", a=a, r=r, attempts=attempt)
    return None


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        targets = [15, 21, 35, 91, 143]
    else:
        targets = [int(x) for x in argv[1:]]
    rng = random.Random(0)
    for N in targets:
        result = shor_reduce(N, rng=rng)
        print(f"N={N}: {result}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
