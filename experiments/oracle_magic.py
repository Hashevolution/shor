"""
oracle_magic.py — T3(b): "magic은 문제의 비선형성에 살고, 오라클 블랙박스가 그걸 숨긴다".

동기 (magic-and-quantum-speedup.md, T3 사다리의 클리포드 끝점 + 오라클-가림 가설):
  쿼리 모델의 속도우위는 *오라클 상대적*이라 오라클 내부의 magic 비용을 숨긴다. 이를
  함수 그래프 상태로 정량화한다:
      |ψ_f⟩ = (1/√N) Σ_x |x⟩|f(x)⟩       (= H^n 후 오라클을 친 상태; Simon/Shor 공통 중간상태)

핵심 결과 (magic.sre2 로 직접 계산):
  1. **M₂(|ψ_f⟩) = 0  ⟺  f 가 𝔽₂-아핀(선형)**. f가 비선형이면 M₂>0, 비선형항(AND) 수에 따라 증가.
  2. **Simon**: 숨은문자열 s의 2-대-1 함수를 *선형*으로 잡을 수 있다 (f(x)=x⊕(xₚ·s)). 그러면
     오라클이 클리포드 → **회로 전체 M₂≡0인데도 지수 쿼리 속도우위**. ⟹ magic은 (쿼리모델)
     지수 속도우위에 *필요하지 않다*.
  3. **Shor**: f(x)=aˣ mod N 은 본질적으로 비선형(곱셈적) → 오라클이 magic을 *강제*하고,
     M₂>0이 문제크기와 함께 증가 (= 측정전 Shor 중간상태의 magic).
  ⟹ "보이는 회로의 magic"이 아니라 "문제의 비선형성"이 자원의 원천. 쿼리 오라클(과 shor.py의
     np.fft 지름길)은 그 magic을 블랙박스로 *숨긴다*. (2507.16543의 Clifford-가림과는 구분:
     이쪽은 *오라클/추상화*가 가리는 것.)

Reproduction:
  python -u -m experiments.oracle_magic
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("oracle_magic_results.txt")


def graph_state(values: list[int], n_in: int, n_out: int) -> np.ndarray:
    """|ψ_f⟩ = (1/√2ⁿ) Σ_x |x⟩|f(x)⟩,  values[x]=f(x) ∈ [0,2^{n_out})."""
    psi = np.zeros(2 ** (n_in + n_out))
    for x in range(2 ** n_in):
        psi[(x << n_out) | values[x]] = 1.0
    return psi / np.sqrt(2 ** n_in)


def simon_linear_oracle(n: int, s: int) -> list[int]:
    """선형 Simon 오라클 f(x)=x⊕(xₚ·s), p=s의 최하위 set bit. 2-대-1, kernel {0,s}."""
    p = (s & -s).bit_length() - 1
    return [x ^ (s if (x >> p) & 1 else 0) for x in range(2 ** n)]


def add_and_terms(values: list[int], n: int, terms: list[tuple[int, int, int]]) -> list[int]:
    """선형 f에 비선형 AND 항 ⊕(x_i·x_j)·e_k 들을 추가."""
    out = list(values)
    for x in range(2 ** n):
        for (i, j, k) in terms:
            if (x >> i) & 1 and (x >> j) & 1:
                out[x] ^= (1 << k)
    return out


def main() -> None:
    t0 = time.time()
    lines: list[str] = []

    def out(s: str = "") -> None:
        print(s)
        lines.append(s + "\n")

    out("# oracle_magic — 함수 그래프상태 |ψ_f⟩=Σ_x|x⟩|f(x)⟩ 의 magic (M₂, bits)")
    out("# 명제: M₂=0 ⟺ f 가 𝔽₂-아핀.  magic = 문제의 비선형성, 오라클이 숨김.")
    out("")

    # 1) Simon 선형 오라클 → M₂ = 0 (여러 n, s)
    out("## 1. Simon 선형 오라클 f(x)=x⊕(xₚ·s)  → M₂=0 기대 (클리포드 끝점)")
    out(f"   {'n':>3} {'s':>4} {'qubits':>7} {'M2':>9}")
    for n in [2, 3, 4]:
        for s in sorted({1, 3, (1 << (n - 1)) | 1}):
            if s >= 2 ** n:
                continue
            vals = simon_linear_oracle(n, s)
            out(f"   {n:>3} {s:>4} {2*n:>7} {sre2(graph_state(vals, n, n)):>9.4f}")
    out("   → 모두 0: Simon은 선형(클리포드) 오라클로 풀린다 ⟹ 지수 쿼리 속도우위 + magic 0.")
    out("")

    # 2) 비선형성 도입 → magic 발생, 항 수에 따라 증가
    out("## 2. 비선형 AND 항 추가 → M₂>0,  비선형항 수에 따라 증가  (n=4, s=3)")
    n, s = 4, 3
    base = simon_linear_oracle(n, s)
    out(f"   {'AND항수':>7} {'M2':>9}")
    for terms in [[], [(0, 1, 2)], [(0, 1, 2), (0, 2, 3)], [(0, 1, 2), (0, 2, 3), (1, 3, 0)]]:
        vals = add_and_terms(base, n, terms)
        out(f"   {len(terms):>7} {sre2(graph_state(vals, n, n)):>9.4f}")
    out("   → magic은 f의 비선형성(2차항)의 양과 함께 자란다.")
    out("")

    # 3) Shor modexp 오라클 = 본질적 비선형 → M₂>0, 문제크기와 함께 증가
    out("## 3. Shor 오라클 f(x)=aˣ mod N  (= 측정전 중간상태) → M₂>0, 크기 따라 증가")
    out(f"   {'N':>4} {'a':>3} {'r':>3} {'qubits':>7} {'M2':>9}")
    for (N, a, t) in [(15, 7, 6), (15, 2, 6), (21, 2, 7), (33, 5, 7), (35, 3, 8)]:
        nout = N.bit_length()
        vals = [pow(a, x, N) for x in range(2 ** t)]
        r = 1
        while pow(a, r, N) != 1:
            r += 1
        out(f"   {N:>4} {a:>3} {r:>3} {t+nout:>7} {sre2(graph_state(vals, t, nout)):>9.4f}")
    out("   → modexp는 곱셈적=비선형 ⟹ 오라클이 magic을 강제. Simon(선형, 0)과 정반대.")
    out("")

    # 결론
    out("## 결론 — 오라클-가림 가설 (T3b)")
    out("   • M₂(|ψ_f⟩)=0 ⟺ f 아핀.  magic의 원천은 회로 표면이 아니라 '문제의 비선형성'.")
    out("   • 쿼리 모델 속도우위는 오라클 내부의 magic 비용을 숨긴다(Simon: 선형 오라클로 0;")
    out("     Shor: modexp 비선형으로 강제). shor.py의 np.fft 지름길이 R_k magic을 숨기는 것과")
    out("     같은 현상 — 블랙박스(오라클/FFT)가 magic을 가린다.")
    out("   • 속도우위 사다리: Simon/BV(아핀,0) → Grover(유한,밀도→0) → Shor(비선형,크기비례→최대).")

    elapsed = time.time() - t0
    out("")
    out(f"# Elapsed: {elapsed:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
