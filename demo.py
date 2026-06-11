"""
쇼어 알고리즘 데모.

실행 예:
    python demo.py                # 기본 데모: N=15, 21, 35
    python demo.py 33 35          # 특정 N들
    python demo.py --dist 15 7    # N=15, a=7에 대한 측정 분포 출력
"""

from __future__ import annotations
import sys
import math

import numpy as np

from classical import shor_reduce, classical_order
from shor import (
    shor_quantum,
    simulate_period_finding,
    measurement_distribution,
    _counting_qubits,
)


def run_factor(N: int, seed: int = 0) -> None:
    print(f"\n── N = {N} ──")
    classical = shor_reduce(N, rng=__import__("random").Random(seed))
    print(f"  [고전 baseline] {classical}")

    quantum = shor_quantum(N, seed=seed)
    print(f"  [양자 시뮬레이션] {quantum}")

    if quantum and quantum.a is not None and quantum.r is not None:
        true_r = classical_order(quantum.a, N)
        ok = "✓" if quantum.r == true_r else f"✗ (true r={true_r})"
        print(f"  검증: a={quantum.a}, 측정 r={quantum.r} {ok}")


def show_distribution(N: int, a: int) -> None:
    t = _counting_qubits(N)
    Q = 1 << t
    print(f"\n── 측정 분포: N={N}, a={a}, t={t} qubits, Q={Q} ──")
    if math.gcd(a, N) != 1:
        print(f"  gcd({a},{N}) = {math.gcd(a, N)} — 위수 정의되지 않음.")
        return
    true_r = classical_order(a, N)
    print(f"  실제 위수 r = {true_r}")
    probs = measurement_distribution(a, N)
    top = np.argsort(probs)[::-1][:8]
    print(f"  상위 측정값 (k, k/Q, 확률):")
    for k in top:
        if probs[k] < 1e-6:
            break
        ratio = k / Q
        # j/r 형태로 추정
        from fractions import Fraction
        f = Fraction(int(k), Q).limit_denominator(N - 1)
        print(f"    k={int(k):5d}  k/Q={ratio:.6f}  P={probs[k]:.4f}  ≈ {f}")
    # 이상적으로는 j·Q/r 부근에 피크가 있어야 함
    print(f"  예상 피크 위치 (j·Q/r): {[round(j*Q/true_r) for j in range(true_r)]}")


def main(argv: list[str]) -> int:
    if "--dist" in argv:
        i = argv.index("--dist")
        N = int(argv[i + 1])
        a = int(argv[i + 2])
        show_distribution(N, a)
        return 0

    targets = [int(x) for x in argv[1:]] if len(argv) > 1 else [15, 21, 35]
    for N in targets:
        run_factor(N)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
