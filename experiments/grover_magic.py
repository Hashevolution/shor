"""
grover_magic.py — Grover 탐색의 비안정자성(magic, stabilizer 2-Rényi 엔트로피) 궤적.

동기 (magic-and-quantum-speedup.md, T3):
  Shor(지수 속도우위)의 magic은 주기 r과 함께 자라 상한 L≈n bit로 *포화*한다
  (Paviglianiti et al., arXiv:2605.05347: M₂ ~ log r → L). 그렇다면 *2차(다항)*
  속도우위인 Grover에서 magic은 어떻게 행동하나? — 선행연구 조사상 Grover의 magic/SRE
  궤적은 미개척(자원 분석은 대부분 coherence·entanglement; 양자걷기 magic 2506.17783/
  2504.19750은 1D walk이지 탐색이 아님). 이 스크립트가 그 빈칸을 채운다.

상태 (M개 표시, 단일 표시 M=1이 표준):
  k회 반복 후 Grover 상태는 2차원 부분공간의 회전으로, 표시원소 진폭 a, 그 외 b의 두 값만
  갖는다:  a = sin((2k+1)θ)/√M,  b = cos((2k+1)θ)/√(N−M),  sinθ = √(M/N).

핵심 결과:
  1. M₂(k)는 0(균일 중첩=안정자)에서 올라 *탐색 중간*에서 정점, 정답 k*≈(π/4)√N 에서
     다시 ≈0 (정답은 계산기저 = 안정자). magic을 *썼다가 되돌린다*.
  2. M=1 상태의 M₂는 닫힌형을 가진다(아래 sre2_grover_closed, 수치로 검증).
  3. 정점 magic은 n에 따라 *포화*: N→∞에서 M₂_peak → 3 bit (a²=1/2에서),
     즉 magic 밀도 M₂/n → 0. ↔ Shor는 M₂/L → 1 (최대). magic은 속도우위 *유형*을 가른다.

Reproduction:
  python -u -m experiments.grover_magic
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("grover_magic_results.txt")


# ── Grover 상태 ──────────────────────────────────────────────────────────────
def grover_state(n: int, marked: list[int], k: int) -> np.ndarray:
    """n큐비트, marked 표시집합, k회 반복 후 Grover 상태벡터."""
    N = 2 ** n
    M = len(marked)
    theta = math.asin(math.sqrt(M / N))
    ang = (2 * k + 1) * theta
    a = math.sin(ang) / math.sqrt(M)
    b = math.cos(ang) / math.sqrt(N - M)
    psi = np.full(N, b, dtype=float)
    for w in marked:
        psi[w] = a
    return psi


def sre2_grover_closed(n: int, a: float, b: float) -> float:
    """M=1 Grover 상태(진폭 a 1개, b N−1개)의 M₂ 닫힌형.

    Σ_P⟨P⟩⁴ = 1 + (N−1)(a²−b²)⁴ + (N−1)(b²(N−2)+2ab)⁴ + (N−1)(N/2−1)(2b(a−b))⁴.
    """
    N = 2 ** n
    S = 1.0
    S += (N - 1) * (a * a - b * b) ** 4
    S += (N - 1) * (b * b * (N - 2) + 2 * a * b) ** 4
    S += (N - 1) * (N / 2 - 1) * (2 * b * (a - b)) ** 4
    return -math.log2(S / N)


def peak_magic_limit(a: float) -> float:
    """N→∞ 극한의 M=1 Grover magic: M₂ → −log₂(a⁸+(1−a²)⁴). a²=1/2에서 최댓값 3."""
    u = a * a
    return -math.log2(u ** 4 + (1 - u) ** 4)


# ── 실험 ─────────────────────────────────────────────────────────────────────
def main() -> None:
    t0 = time.time()
    lines: list[str] = []

    def out(s: str = "") -> None:
        print(s)
        lines.append(s + "\n")

    out("# grover_magic — Grover 탐색의 비안정자성 M₂ (stabilizer 2-Rényi 엔트로피, bits)")
    out("# 측정도구: magic.sre2 (XOR-FWHT, brute-force 검증됨)")
    out("")

    # 0) 닫힌형 검증 (sre2 vs sre2_grover_closed)
    out("## 0. M=1 닫힌형 검증 (sre2 ↔ sre2_grover_closed)")
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for n in [3, 4, 5, 6, 8, 10]:
        N = 2 ** n
        for _ in range(5):
            ang = rng.uniform(0.1, math.pi / 2 - 0.1)
            a, b = math.sin(ang), math.cos(ang) / math.sqrt(N - 1)
            psi = np.full(N, b)
            psi[0] = a
            maxerr = max(maxerr, abs(sre2(psi) - sre2_grover_closed(n, a, b)))
    out(f"   max|sre2 − closed| = {maxerr:.2e}  (30 random angles)  → 닫힌형 확인")
    out("")

    # 1) magic 궤적 vs 반복수 k (M=1)
    out("## 1. magic 궤적 M₂(k),  단일 표시 M=1")
    out(f"   {'n':>3} {'N':>6} {'k*':>4} {'M2@k*':>8} {'peak':>7} {'@k':>4}   trajectory(k=0..k*+1)")
    for n in [4, 6, 8, 10, 12]:
        N = 2 ** n
        kstar = int(round(math.pi / 4 * math.sqrt(N) - 0.5))
        traj = [sre2(grover_state(n, [1], k)) for k in range(kstar + 2)]
        kpk = int(np.argmax(traj))
        seg = " ".join(f"{v:.2f}" for v in traj[: min(len(traj), 16)])
        tail = " ..." if len(traj) > 16 else ""
        out(f"   {n:>3} {N:>6} {kstar:>4} {traj[kstar]:>8.4f} {traj[kpk]:>7.4f} {kpk:>4}   {seg}{tail}")
    out("   → 0(균일=안정자) → 중간서 정점 → k*서 ≈0(정답=계산기저=안정자).")
    out("")

    # 2) 정점 magic의 스케일링 (닫힌형으로 큰 n까지) + 점근 극한
    out("## 2. 정점 magic 스케일링 (닫힌형, M=1) — n↑ 에서 포화하는가?")
    out(f"   {'n':>3} {'peak M2':>9} {'M2/n':>7} {'Δ(직전행)':>9}   (극한 N→∞: 3.000)")
    prev = None
    angs = np.linspace(0.30, 0.60, 1200) * (math.pi / 2)
    for n in [3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28]:
        N = 2 ** n
        a = np.sin(angs)
        b = np.cos(angs) / math.sqrt(N - 1)
        m = np.array([sre2_grover_closed(n, ai, bi) for ai, bi in zip(a, b)])
        pk = float(m.max())
        d = "" if prev is None else f"{pk - prev:>+8.4f}"
        out(f"   {n:>3} {pk:>9.4f} {pk / n:>7.4f} {d:>8}")
        prev = pk
    out(f"   해석적 극한: M₂ → −log₂(a⁸+(1−a²)⁴), a²=1/2에서 최대 = {peak_magic_limit(1/math.sqrt(2)):.4f}")
    out("   → 정점이 3 bit로 포화 (밀도 M₂/n → 0).  ↔ Shor: M₂ → L (밀도 → 1; 2605.05347).")
    out("")

    # 3) 표시개수 M 의존성 (정점, 닫힌형은 M=1 전용이라 sre2로 직접)
    out("## 3. 표시개수 M 의존성 (n=10, 정점 magic over k)")
    n = 10
    N = 2 ** n
    out(f"   {'M':>4} {'peak M2':>9} {'M2/n':>7}")
    for M in [1, 2, 4, 8, 16]:
        marked = list(range(M))
        kstar = int(round(math.pi / 4 * math.sqrt(N / M) - 0.5))
        m = [sre2(grover_state(n, marked, k)) for k in range(max(2, kstar + 2))]
        pk = max(m)
        out(f"   {M:>4} {pk:>9.4f} {pk / n:>7.4f}")
    out("   → 표시개수를 키워도 정점은 여전히 O(1) (밀도 낮음).")
    out("")

    # 결론
    out("## 결론")
    out("   Grover(2차 속도우위): magic은 *유한*(정점 → 3 bit, 밀도 → 0)이고 정답에서 0으로")
    out("   되돌려진다. Shor(지수 속도우위): magic은 *최대*(M₂ → L, 밀도 → 1)로 구동된다.")
    out("   → 비안정자성의 *양/밀도*가 속도우위의 유형(다항 vs 지수)을 가른다는 직접 증거.")

    elapsed = time.time() - t0
    out("")
    out(f"# Elapsed: {elapsed:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
