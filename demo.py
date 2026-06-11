"""
쇼어 알고리즘 데모.

실행 예:
    python demo.py                       # 기본 데모: N=15, 21, 35
    python demo.py 33 35                 # 특정 N 들
    python demo.py --dist 15 7           # N=15, a=7 측정 분포 출력
    python demo.py --multi 91            # 다중 base 모드
    python demo.py --compare 33 35 77    # 단일 vs 다중 회수율 비교
"""

from __future__ import annotations
import sys
import math
import random
import time

import numpy as np

from classical import shor_reduce, classical_order
from shor import (
    shor_quantum,
    simulate_period_finding,
    measurement_distribution,
    quantum_order,
    _counting_qubits,
)
from multi_base import shor_quantum_multi


def run_factor(N: int, seed: int = 0) -> None:
    print(f"\n── N = {N} ──")
    classical = shor_reduce(N, rng=random.Random(seed))
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
        from fractions import Fraction
        f = Fraction(int(k), Q).limit_denominator(N - 1)
        print(f"    k={int(k):5d}  k/Q={ratio:.6f}  P={probs[k]:.4f}  ≈ {f}")
    print(f"  예상 피크 위치 (j·Q/r): {[round(j*Q/true_r) for j in range(true_r)]}")


def run_multi(N: int, seed: int = 0) -> None:
    print(f"\n── 다중 base: N = {N} ──")
    result, state = shor_quantum_multi(N, seed=seed)
    print(f"  결과: {result}")
    print(f"  양자 측정: {state.measurements}, 고전 회수: {state.classical_recoveries}")
    print(f"  누적 L: {state.L}, 본 위수: {state.orders}")


def _count_single_base(N: int, seed: int):
    """단일 base shor_quantum을 측정 카운터와 함께 실행."""
    import shor as shor_mod
    original = shor_mod.simulate_period_finding
    count = [0]

    def counting_wrapper(*args, **kwargs):
        count[0] += 1
        return original(*args, **kwargs)

    shor_mod.simulate_period_finding = counting_wrapper
    try:
        result = shor_mod.shor_quantum(N, seed=seed)
        return result, count[0]
    finally:
        shor_mod.simulate_period_finding = original


def compare(targets: list[int], trials: int = 50) -> None:
    """단일-base vs 다중-base 인수분해 비교.

    같은 시드 셋으로 두 방법을 돌려 성공률·평균 측정 횟수·평균 시간을 표로 출력.
    측정 횟수는 simulate_period_finding 호출 횟수로 정확히 셈.
    """
    print(f"\n── 비교 ({trials} trials, seeds 0..{trials-1}) ──")
    print(f"  {'N':>4}  {'method':<12}  {'success':>7}  "
          f"{'meas/trial':>10}  {'ms/trial':>9}")

    for N in targets:
        s_succ = 0
        s_meas_total = 0
        s_time = 0.0
        for seed in range(trials):
            t0 = time.perf_counter()
            res, n_meas = _count_single_base(N, seed)
            s_time += time.perf_counter() - t0
            if res is not None:
                s_succ += 1
            s_meas_total += n_meas

        m_succ = 0
        m_meas_total = 0
        m_time = 0.0
        for seed in range(trials):
            t0 = time.perf_counter()
            res, state = shor_quantum_multi(N, seed=seed)
            m_time += time.perf_counter() - t0
            if res is not None:
                m_succ += 1
                m_meas_total += state.measurements

        print(f"  {N:>4}  {'single-base':<12}  "
              f"{s_succ:>3}/{trials:<3}  "
              f"{s_meas_total/trials:>10.2f}  "
              f"{s_time/trials*1000:>9.2f}")
        print(f"  {N:>4}  {'multi-base':<12}  "
              f"{m_succ:>3}/{trials:<3}  "
              f"{m_meas_total/trials:>10.2f}  "
              f"{m_time/trials*1000:>9.2f}")


def compare_period_finding(N: int, trials: int = 200, seed: int = 0) -> None:
    """순수 위수 회수 확률 비교 (인수분해와 분리).

    각 trial: 랜덤 base a 선택, 측정 1회, 세 가지 후처리로 r 회수 시도.
      (A) limit_denominator(N-1)  — 기존 simulate_period_finding 방식
      (B) 모든 연분수 수렴값 분모 시도 (Knill-Mosca)
      (C) (B) + 이전 trial 들의 누적 L 약수 시도 (다중 base)

    회수 성공: 얻은 d 가 a 의 실제 위수 r_a 와 일치.
    """
    from fractions import Fraction
    from multi_base import (
        convergent_denominators,
        divisors,
        minimize_order,
        MultiBaseState,
    )
    from classical import classical_order

    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)
    state = MultiBaseState()  # (C) 용 누적 상태

    succ_A = succ_B = succ_C = 0
    bases_tried = 0

    for _ in range(trials):
        for _retry in range(50):
            a = rng_py.randrange(2, N)
            if math.gcd(a, N) == 1:
                break
        else:
            continue
        bases_tried += 1
        true_r = classical_order(a, N)

        m = simulate_period_finding(a, N, rng=rng_np)

        # (A) limit_denominator 한 개
        d_A = Fraction(m.k, m.Q).limit_denominator(N - 1).denominator
        ok_A = d_A > 0 and pow(a, d_A, N) == 1 and minimize_order(a, N, d_A) == true_r
        if ok_A:
            succ_A += 1

        # (B) 모든 수렴값 시도
        cands_B = [d for d in convergent_denominators(m.k, m.Q, N - 1)
                   if d > 0 and pow(a, d, N) == 1]
        ok_B = bool(cands_B) and minimize_order(a, N, min(cands_B)) == true_r
        if ok_B:
            succ_B += 1

        # (C) (B) + 누적 L 약수
        cands_C = set(cands_B)
        if state.L > 1:
            cands_C.update(d for d in divisors(state.L) if pow(a, d, N) == 1)
        ok_C = bool(cands_C) and minimize_order(a, N, min(cands_C)) == true_r
        if ok_C:
            succ_C += 1
            state.update(a, true_r)
        elif ok_B:
            state.update(a, true_r)
        elif ok_A:
            state.update(a, true_r)

    print(f"\n── 1회 측정 위수 회수율: N={N}, {bases_tried} trials ──")
    print(f"  (A) limit_denominator      : {succ_A:>4}/{bases_tried} "
          f"= {succ_A/bases_tried:6.1%}")
    print(f"  (B) 모든 수렴값             : {succ_B:>4}/{bases_tried} "
          f"= {succ_B/bases_tried:6.1%}")
    print(f"  (C) 수렴값 + 누적 L         : {succ_C:>4}/{bases_tried} "
          f"= {succ_C/bases_tried:6.1%}")
    print(f"  최종 누적 L = {state.L}")


def main(argv: list[str]) -> int:
    if "--dist" in argv:
        i = argv.index("--dist")
        N = int(argv[i + 1])
        a = int(argv[i + 2])
        show_distribution(N, a)
        return 0

    if "--multi" in argv:
        i = argv.index("--multi")
        Ns = [int(x) for x in argv[i + 1:]] if i + 1 < len(argv) else [33]
        for N in Ns:
            run_multi(N)
        return 0

    if "--compare" in argv:
        i = argv.index("--compare")
        Ns = [int(x) for x in argv[i + 1:]] if i + 1 < len(argv) else [15, 21, 33, 35]
        compare(Ns)
        return 0

    if "--period" in argv:
        i = argv.index("--period")
        Ns = [int(x) for x in argv[i + 1:]] if i + 1 < len(argv) else [33, 77, 143, 209]
        for N in Ns:
            compare_period_finding(N)
        return 0

    targets = [int(x) for x in argv[1:]] if len(argv) > 1 else [15, 21, 35]
    for N in targets:
        run_factor(N)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
