"""
oracle_tcount_magic.py — M4 과제 C+D: "오라클이 숨긴 magic = 게이트분해의 비클리포드(T) 비용".

동기 (magic-and-quantum-speedup.md T3(b); 따름정리 1):
  그래프상태 |ψ_f⟩=(1/√N)Σ_x|x⟩|f(x)⟩ 는 M₂=0 ⟺ f 𝔽₂-아핀(따름정리 1). 한편 오라클
  U_f:|x⟩|y⟩→|x⟩|y⊕f(x)⟩ 를 게이트로 짜면, f의 ANF(대수정규형) 단항식별로:
     • 1차항 x_i        → CNOT      (클리포드, T=0)
     • 2차항 x_i x_j    → Toffoli   (비클리포드, T 필요)
     • ≥3차항          → 다중제어  (T 더 필요)
  ⟹ **오라클의 비클리포드(T) 비용은 f의 *비선형(차수≥2) ANF*가 통제**한다. 따라서

      M₂(|ψ_f⟩) > 0  ⟺  f 비선형  ⟺  U_f 가 비클리포드(T) 게이트를 요구.

  즉 "상태가 품은 magic"과 "오라클 게이트분해가 치르는 T 비용"은 *같은 비선형성*에서 나오며,
  둘 다 f 아핀에서만 0이다. 비선형 ANF 단항식 수 T_proxy(=각 출력비트 ANF의 차수≥2 항 합)를
  Toffoli/T 비용의 대용(하한 신호)으로 쓴다.

차별화 (2507.16543 "Quantum Dark Magic"의 Clifford-가림과 구분):
  그쪽은 *Clifford 켤레*가 magic을 가린다(permutation-agnostic 거리로 드러냄). 여기서는
  *오라클/FFT 추상화*가 비선형성의 magic을 블랙박스로 가린다 — 메커니즘이 다르다.

Reproduction:
  python -u -m experiments.oracle_tcount_magic
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("oracle_tcount_magic_results.txt")


# ── 도구 ─────────────────────────────────────────────────────────────────────
def graph_state(values, n_in, n_out):
    psi = np.zeros(2 ** (n_in + n_out))
    for x in range(2 ** n_in):
        psi[(x << n_out) | values[x]] = 1.0
    return psi / np.sqrt(2 ** n_in)


def anf_coeffs(truth: np.ndarray, n: int) -> np.ndarray:
    """단일 출력비트 진리표 → ANF 계수(이진 Möbius 변환).  coeff[S]=⊕_{x⊆S} g(x)."""
    a = truth.astype(np.uint8).copy()
    for i in range(n):
        bit = 1 << i
        idx = np.arange(2 ** n)
        sel = (idx & bit) != 0
        a[sel] ^= a[idx[sel] ^ bit]
    return a


def nonlinear_monomials(values, n_in, n_out) -> int:
    """f의 모든 출력비트 ANF에서 차수≥2 단항식 개수 합 (= Toffoli/T 비용 대용)."""
    vals = np.asarray(values)
    deg = np.array([bin(s).count("1") for s in range(2 ** n_in)])
    total = 0
    for k in range(n_out):
        truth = ((vals >> k) & 1).astype(np.uint8)
        coeff = anf_coeffs(truth, n_in)
        total += int(np.sum((coeff != 0) & (deg >= 2)))
    return total


def simon_linear(n, s):
    p = (s & -s).bit_length() - 1
    return [x ^ (s if (x >> p) & 1 else 0) for x in range(2 ** n)]


def add_and_terms(values, n, terms):
    out = list(values)
    for x in range(2 ** n):
        for (i, j, k) in terms:
            if (x >> i) & 1 and (x >> j) & 1:
                out[x] ^= (1 << k)
    return out


# ── 실험 ─────────────────────────────────────────────────────────────────────
def main():
    t0 = time.time()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s + "\n")

    out("# oracle_tcount_magic — 상태 magic(M₂) ↔ 오라클 비클리포드(T) 비용 (비선형 ANF 단항식 수)")
    out("# 명제: M₂(|ψ_f⟩)>0 ⟺ f 비선형 ⟺ U_f 가 T 게이트 요구.  T_proxy = 차수≥2 ANF 단항식 수.")
    out("")

    # 0) ANF 도구 검증: 선형 f → 비선형 단항식 0
    out("## 0. 도구 검증 — 선형/아핀 f 는 비선형 ANF 단항식 0 (⟹ 클리포드 오라클)")
    rng = np.random.default_rng(0)
    ok = True
    for n in [3, 4, 5]:
        A = rng.integers(0, 2, size=(n, n))
        b = rng.integers(0, 2, size=n)
        vals = []
        for x in range(2 ** n):
            xb = np.array([(x >> i) & 1 for i in range(n)])
            y = (A @ xb + b) % 2
            vals.append(int(sum(int(y[i]) << i for i in range(n))))
        nl = nonlinear_monomials(vals, n, n)
        ok = ok and nl == 0
        out(f"   랜덤 아핀 f (n={n}): 비선형 단항식 = {nl},  M₂ = {sre2(graph_state(vals, n, n)):.4f}")
    out(f"   → 아핀 f 는 T_proxy=0 이고 M₂=0 (둘 다 0).  도구 정상: {ok}")
    out("")

    # 1) T_proxy=0 ⟺ M₂=0,  그리고 동반 증가
    out("## 1. 비선형 항 추가 → T_proxy 와 M₂ 가 함께 0에서 켜지고 함께 증가 (n=4, s=3)")
    n, s = 4, 3
    base = simon_linear(n, s)
    out(f"   {'추가 AND항':>9} {'T_proxy':>8} {'M2':>9}")
    for terms in [[], [(0, 1, 2)], [(0, 1, 2), (0, 2, 3)],
                  [(0, 1, 2), (0, 2, 3), (1, 3, 0)],
                  [(0, 1, 2), (0, 2, 3), (1, 3, 0), (2, 3, 1)]]:
        vals = add_and_terms(base, n, terms)
        tp = nonlinear_monomials(vals, n, n)
        out(f"   {len(terms):>9} {tp:>8} {sre2(graph_state(vals, n, n)):>9.4f}")
    out("   → T_proxy=0 ⟺ M₂=0; 비선형성↑ 이면 둘 다 ↑ (상태 magic ↔ 오라클 T 비용 동반).")
    out("")

    # 2) 끝점 대조: Simon(선형) vs Shor(modexp 비선형)
    out("## 2. 끝점 — Simon(선형: T_proxy 0, M₂ 0)  vs  Shor modexp(비선형: T_proxy↑, M₂↑)")
    out(f"   {'알고리즘':<18} {'in':>3} {'out':>4} {'T_proxy':>8} {'M2':>9}")
    for n in [3, 4]:
        vals = simon_linear(n, 3)
        out(f"   {'Simon (선형오라클)':<18} {n:>3} {n:>4} "
            f"{nonlinear_monomials(vals, n, n):>8} {sre2(graph_state(vals, n, n)):>9.4f}")
    for (N, a, t) in [(15, 7, 5), (15, 2, 5), (21, 2, 6)]:
        nout = N.bit_length()
        vals = [pow(a, x, N) for x in range(2 ** t)]
        tp = nonlinear_monomials(vals, t, nout)
        out(f"   {'Shor a=%d mod %d' % (a, N):<18} {t:>3} {nout:>4} "
            f"{tp:>8} {sre2(graph_state(vals, t, nout)):>9.4f}")
    out("   → Simon: 선형이라 T_proxy=0·M₂=0 인데도 지수 *쿼리* 우위(magic 불필요).")
    out("     Shor: modexp 비선형이라 오라클이 T 비용·magic 둘 다 강제. 같은 비선형성이 원천.")
    out("")

    out("## 결론 — 과제 C+D 통합")
    out("   • 상태 magic M₂(|ψ_f⟩) 과 오라클 게이트분해의 비클리포드(T) 비용은 *같은 비선형성*에서")
    out("     나오며 둘 다 f 아핀에서만 0 (따름정리 1 + ANF 차수≥2 ⟺ Toffoli/T).")
    out("   • '쿼리 오라클이 magic을 숨긴다'(T3b) = 'T 비용을 블랙박스로 숨긴다'와 동치 — Simon(0)")
    out("     vs Shor(>0)로 정량화. (2507.16543의 Clifford-가림과는 다른 메커니즘: 오라클 추상화.)")

    out("")
    out(f"# Elapsed: {time.time() - t0:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
