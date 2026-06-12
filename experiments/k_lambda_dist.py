"""
경험적 K_λ 분포 측정 (Phase 1 검증용).

K_λ := L = lcm(r_{a_1}, …, r_{a_K}) 가 처음 λ(N) 이 될 때의 K.

다양한 반소수 N=pq 에서 trials 회 시뮬레이션해 K_λ 의 평균/꼬리/분포를 수집.
Phase 1 의 정량적 정리 (정리 2) 의 예측과 비교한다.

알고리즘 노이즈 (양자 측정) 와 무관하게, **순수 군이론적** K_λ 만 측정 — 매 base
의 실제 위수 (classical_order) 를 사용. 노이즈 layer 는 paper §4 의 (C)-determinism
정리가 이미 cover.

실행:
    python -m experiments.k_lambda_dist
    python -m experiments.k_lambda_dist 15 21 35 77 143 209 437 1147 2491 4087
"""

from __future__ import annotations
import math
import random
import sys
from dataclasses import dataclass

from classical import classical_order


@dataclass
class KLambdaStats:
    N: int
    p: int
    q: int
    lam: int                  # λ(N)
    omega_lam: int            # 서로 다른 prime divisor 수 of λ(N)
    trials: int
    mean: float
    median: float
    p90: int                  # 90 분위
    p99: int                  # 99 분위
    max: int
    bound_c: float            # 정리 2 (c) 의 E[K_λ] 상한: 1 + Σ 1/(ℓ^s - 1)
    bound_b: float            # 정리 2 (b) 식 상한 (단순): log₂(ω) + 2


def carmichael_semiprime(p: int, q: int) -> int:
    """N=pq 반소수의 Carmichael λ(N) = lcm(p-1, q-1)."""
    return math.lcm(p - 1, q - 1)


def omega(n: int) -> int:
    """n 의 서로 다른 소인수 수."""
    if n < 2:
        return 0
    count = 0
    x = n
    d = 2
    while d * d <= x:
        if x % d == 0:
            count += 1
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        count += 1
    return count


def factor_semiprime(N: int) -> tuple[int, int]:
    """N 의 소인수 (p, q). N 이 반소수가 아니면 raise."""
    for p in range(2, int(N**0.5) + 1):
        if N % p == 0:
            q = N // p
            return p, q
    raise ValueError(f"{N} is prime, not a semiprime")


def v_ell(n: int, ell: int) -> int:
    """ℓ-adic valuation v_ℓ(n)."""
    v = 0
    while n % ell == 0:
        n //= ell
        v += 1
    return v


def primes_dividing(n: int) -> list[int]:
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


def theorem2_bounds(p: int, q: int) -> tuple[float, float]:
    """정리 2 의 상한값 계산.

    Returns:
        bound_c: E[K_λ] 의 sharp 상한 = 1 + Σ_ℓ 1/(ℓ^s_ℓ - 1).
        bound_b: simple 상한 = log₂(ω(λ)) + 2.
    """
    lam = math.lcm(p - 1, q - 1)
    primes = primes_dividing(lam)
    bound_c = 1.0
    for ell in primes:
        v_max = v_ell(lam, ell)
        v_p = v_ell(p - 1, ell)
        v_q = v_ell(q - 1, ell)
        s_ell = (1 if v_p == v_max else 0) + (1 if v_q == v_max else 0)
        bound_c += 1.0 / (ell ** s_ell - 1)
    bound_b = math.log2(len(primes)) + 2 if len(primes) >= 1 else 2.0
    return bound_c, bound_b


def measure_k_lambda(
    N: int, trials: int = 200, seed: int = 0, max_K: int = 200,
) -> KLambdaStats:
    """N=pq 에서 K_λ 의 경험 분포 수집."""
    p, q = factor_semiprime(N)
    lam = carmichael_semiprime(p, q)
    omega_lam = omega(lam)

    samples: list[int] = []
    rng = random.Random(seed)

    for trial in range(trials):
        L = 1
        for K in range(1, max_K + 1):
            for _retry in range(50):
                a = rng.randrange(2, N)
                if math.gcd(a, N) == 1:
                    break
            else:
                continue
            r_a = classical_order(a, N)
            L = math.lcm(L, r_a)
            if L == lam:
                samples.append(K)
                break
        else:
            samples.append(max_K)  # 도달 실패 (희박)

    samples.sort()
    n = len(samples)
    bound_c, bound_b = theorem2_bounds(p, q)
    return KLambdaStats(
        N=N, p=p, q=q, lam=lam, omega_lam=omega_lam, trials=n,
        mean=sum(samples) / n,
        median=samples[n // 2],
        p90=samples[int(n * 0.9)],
        p99=samples[min(int(n * 0.99), n - 1)],
        max=samples[-1],
        bound_c=bound_c,
        bound_b=bound_b,
    )


def print_table(stats: list[KLambdaStats]) -> None:
    print(f"{'N':>7}  {'p':>5} {'q':>5}  {'λ(N)':>9}  {'ω(λ)':>4}  "
          f"{'mean':>5}  {'med':>4}  {'p90':>4}  {'p99':>4}  {'max':>4}  "
          f"{'thm2c':>6}  {'thm2b':>6}  {'ok':>3}")
    for s in stats:
        ok = "✓" if s.mean <= s.bound_c + 1e-9 else "✗"
        print(f"{s.N:>7}  {s.p:>5} {s.q:>5}  {s.lam:>9}  {s.omega_lam:>4}  "
              f"{s.mean:>5.2f}  {s.median:>4}  {s.p90:>4}  {s.p99:>4}  {s.max:>4}  "
              f"{s.bound_c:>6.3f}  {s.bound_b:>6.3f}  {ok:>3}")


DEFAULT_NS = [
    15, 21, 33, 35, 77, 91, 143, 187, 209, 221, 247, 323, 391, 437,
    1147, 2491, 4087,
]


def main(argv):
    Ns = [int(x) for x in argv[1:]] if len(argv) > 1 else DEFAULT_NS
    trials = 1000
    print(f"# K_λ 경험 분포 — {trials} trials per N")
    stats = []
    for N in Ns:
        try:
            s = measure_k_lambda(N, trials=trials)
        except ValueError as e:
            print(f"# skip {N}: {e}")
            continue
        stats.append(s)
    print_table(stats)


if __name__ == "__main__":
    main(sys.argv)
