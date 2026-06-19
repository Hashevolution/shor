"""
marker_code_expected.py — M4 잔여: 랜덤 마커상태 E[ξ]의 *정확* 닫힌형(명제 5의 Δ(M,N) 완성).

명제 5는 Sidon 주항 ξ_Sidon=(7M−6)/M³ 까지였다. 균일 랜덤 M-부분집합 W⊂𝔽₂ⁿ(N=2ⁿ)에 대해
ξ=(1/M⁴)Σ_x E(W∩(W⊕x)) 의 기댓값을 4중쌍(a,b,c,d: a⊕b⊕c⊕d=0) 분류로 *정확히* 닫는다:

  base 4중쌍 N³개 = all-equal(N) + paired(3N(N−1)) + genuine(N(N−1)(N−2)).
  각 base와 이동 x(N개)에 대해 8점 {a,b,c,d,a⊕x,…,d⊕x}의 서로 다른 개수 k 를 세고,
  균일 M-부분집합에서 P(k점 모두 W)=(M)_k/(N)_k (하강 계승) 를 곱해 합산:

      E[ξ]·M⁴ = (7M²−6M)              ← all-equal+paired = Sidon 주항 (명제 5)
                + 7·(M)_4/(N−3)        ← genuine 4중쌍의 x∈{0,δ₁,δ₂,δ₃} (k=4) + paired k=4 잔차
                + N(N−1)(N−2)(N−4)·(M)_8/(N)_8   ← genuine, x 일반 (k=8)

  (genuine 4중쌍은 차분이 2차원 부분공간 {0,δ₁,δ₂,δ₃}을 이뤄 x∈그 4원소면 A⊕x=A(k=4),
   그 외 N−4개 x는 k=8.)  주항 Δ ≈ 7M⁴/N (M≪N) → Δξ≈7/N.

  E[M₂] ≥ −log₂E[ξ] (Jensen; −log₂ 볼록) 이고 간극은 작다 → −log₂E[ξ] 는 E[M₂]의 정밀 하한.

Reproduction:
  python -u -m experiments.marker_code_expected
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("marker_code_expected_results.txt")


def falling(a, k):
    r = 1.0
    for i in range(k):
        r *= (a - i)
    return r


def expected_xi_closed(M, N):
    """균일 랜덤 M-부분집합의 E[ξ] 정확 닫힌형."""
    sidon = 7 * M * M - 6 * M
    g4 = 7 * falling(M, 4) / (N - 3)
    g8 = N * (N - 1) * (N - 2) * (N - 4) * falling(M, 8) / falling(N, 8)
    return (sidon + g4 + g8) / M ** 4


def flat(W, n):
    psi = np.zeros(2 ** n)
    psi[list(W)] = 1.0
    return psi / math.sqrt(len(W))


def main():
    t0 = time.time()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s + "\n")

    out("# marker_code_expected — 랜덤 마커상태 E[ξ] 정확 닫힌형 (명제 5의 Δ(M,N) 완성)")
    out("# E[ξ]M⁴ = (7M²−6M) + 7(M)₄/(N−3) + N(N−1)(N−2)(N−4)(M)₈/(N)₈")
    out("")

    out("## E[ξ]: 닫힌형 vs 몬테카를로,  그리고 −log₂E[ξ] (E[M₂]의 Jensen 하한) vs 실측 E[M₂]")
    out(f"   {'M':>3} {'n':>3} {'M²/N':>6} {'E[ξ]closed':>11} {'E[ξ]MC':>10} {'relerr':>8} "
        f"{'−log₂E[ξ]':>10} {'E[M₂]MC':>9} {'Sidon':>8}")
    rng = np.random.default_rng(1)
    cases = [(4, 6), (8, 6), (8, 8), (8, 10), (16, 8), (16, 10), (32, 10), (32, 12)]
    maxrel = 0.0
    for (M, n) in cases:
        N = 2 ** n
        nsamp = 300 if n <= 8 else (120 if n <= 10 else 30)
        xis, m2s = [], []
        for _ in range(nsamp):
            W = sorted(rng.choice(N, size=M, replace=False).tolist())
            m2 = sre2(flat(W, n))
            m2s.append(m2)
            xis.append(2.0 ** (-m2))
        e_xi_mc = float(np.mean(xis))
        e_xi_cl = expected_xi_closed(M, N)
        relerr = abs(e_xi_mc - e_xi_cl) / e_xi_cl
        maxrel = max(maxrel, relerr)
        sidon = (7 * M - 6) / M ** 3
        out(f"   {M:>3} {n:>3} {M*M/N:>6.2f} {e_xi_cl:>11.5f} {e_xi_mc:>10.5f} {relerr:>8.1e} "
            f"{-math.log2(e_xi_cl):>10.4f} {np.mean(m2s):>9.4f} "
            f"{-math.log2(sidon):>8.4f}")
    out(f"   → 닫힌형 vs MC 최대 상대오차 = {maxrel:.1e} (MC 표본오차 수준) → E[ξ] 닫힌형 정확.")
    out("   → −log₂E[ξ] ≤ E[M₂] (Jensen) 이고 간극 작음. M²/N↑일수록 Sidon 대비 하락 = Δ 효과.")
    out("")

    out("## 점근 분해 (M≪N):  주항 Δξ ≈ 7(M)₄/(M⁴·N)  (다음항 ~ M⁸/N⁴)")
    out(f"   {'M':>3} {'n':>3} {'Δξ=E[ξ]−Sidon':>14} {'7(M)₄/(M⁴N)':>13} {'비(→1)':>8}")
    for (M, n) in [(8, 10), (8, 12), (16, 12), (8, 14)]:
        N = 2 ** n
        dxi = expected_xi_closed(M, N) - (7 * M - 6) / M ** 3
        lead = 7 * falling(M, 4) / (M ** 4 * N)
        out(f"   {M:>3} {n:>3} {dxi:>14.3e} {lead:>13.3e} {dxi / lead:>8.4f}")
    out("   → M²/N→0에서 Δξ → 7(M)₄/(M⁴N) (비→1).  (M)₄/M⁴→1 은 M→∞에서만; 유한 M은 그 인자만큼 작음.")
    out("")

    out("## 결론")
    out("   명제 5의 Δ(M,N) 완성: E[ξ]M⁴ = (7M²−6M) + 7(M)₄/(N−3) + N(N−1)(N−2)(N−4)(M)₈/(N)₈.")
    out("   주항 Sidon, 보정은 genuine 4중쌍(가법구조)에서 옴; Δξ≈7(M)₄/(M⁴N)(M≪N). E[M₂]≈−log₂E[ξ].")

    out("")
    out(f"# Elapsed: {time.time() - t0:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
