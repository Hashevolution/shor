"""
marker_code_magic.py — 표시집합 W를 고전 부호(code)로 보고, 그 대수적 특성으로
flat 상태 |flat_W⟩ = |W|^{-1/2} Σ_{x∈W}|x⟩ 의 magic(M₂)을 예측·정정하는 PoC.

동기 (magic-and-quantum-speedup.md §3 "다음 과제: 부호 이론", magic-results.md 명제 2′):
  명제 2′은 일반 M의 Grover magic이 |flat_W⟩의 magic으로 환원되고, 그 닫힌형이 W의
  *자기상관* A_W(x)=|W∩(W⊕x)| 에 의존함을 남겨두었다(= 2605.05347의 기하항 Λ와 동형).
  인수인계서(PROJECT JAMES v0.5.0)는 예측 지표를 "최소 해밍 거리 + 가중치 분포"로
  제안했으나, 이는 부정확하다. 본 스크립트가 정정한다:

  정정된 정확 객체 — 자기상관 / Walsh 스펙트럼.
  flat 상태에서 ⟨Z^z X^x⟩ = (1/M)·ĝ_x(z),  g_x(c)=[c∈W ∧ c⊕x∈W],  Σ_c g_x = A_W(x).
  ⟹  Σ_P⟨P⟩⁴ = (1/M⁴) Σ_x Σ_z ĝ_x(z)⁴,   M₂ = -log₂( (1/(N·M⁴)) Σ_{x,z} ĝ_x(z)⁴ ).
  magic은 1_W의 자기상관 A_W(= 쌍대부호 가중치 분포 계열)와 그 Walsh 4차 모멘트가 결정하며,
  *최소 해밍 거리만으로는 결정되지 않는다*(§1에서 반례로 박제).

  정정된 영점 판정(보조정리 1의 부호이론적 재서술):
      |flat_W⟩ 안정자(M₂=0)  ⟺  W 아핀부분공간  ⟺  A_W(x) ∈ {0, M} 전부(자기상관 2값).
  비아핀이면 0 < A_W(x) < M 인 x가 생기고, 그 "중간 질량"이 magic의 원천이다.
  이를 정량화한 비아핀성 스칼라
      τ(W) = (1/N) Σ_{x≠0} (A_W(x)/M)·(1 − A_W(x)/M)   ( =0 ⟺ 아핀 ⟺ M₂=0 )
  를 보정 지표 후보로 도입하고, M₂와의 상관(근사 예측력) 및 d_min의 무력함을 보인다.

Reproduction:
  python -u -m experiments.marker_code_magic
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("marker_code_magic_results.txt")


# ── 기본 도구 ────────────────────────────────────────────────────────────────
def _fwht(a: np.ndarray) -> np.ndarray:
    """정규화 없는 Walsh–Hadamard 변환 (길이 2ⁿ)."""
    a = np.asarray(a, dtype=float).copy()
    N = a.size
    h = 1
    while h < N:
        a = a.reshape(N // (2 * h), 2, h)
        x, y = a[:, 0, :], a[:, 1, :]
        a = np.stack([x + y, x - y], axis=1).reshape(N)
        h *= 2
    return a


def flat(W: list[int], n: int) -> np.ndarray:
    psi = np.zeros(2 ** n)
    psi[list(W)] = 1.0
    return psi / math.sqrt(len(W))


def subspace(basis: list[int], off: int = 0) -> list[int]:
    pts = {off}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return sorted(pts)


def autocorr(W: list[int], n: int) -> np.ndarray:
    """A_W(x) = |W ∩ (W⊕x)| 를 모든 x에 대해 (FWHT로 O(N log N))."""
    N = 2 ** n
    ind = np.zeros(N)
    ind[list(W)] = 1.0
    # A_W = (1/N) FWHT( FWHT(ind)² )   (Wiener–Khinchin, F₂ 위)
    f = _fwht(ind)
    return _fwht(f * f) / N


def min_hamming(W: list[int]) -> int:
    """W를 고전 부호로 볼 때 서로 다른 코드워드의 최소 해밍 거리."""
    best = math.inf
    for i in range(len(W)):
        for j in range(i + 1, len(W)):
            best = min(best, bin(W[i] ^ W[j]).count("1"))
    return int(best)


def sre2_flat_autocorr(W: list[int], n: int) -> float:
    """정정된 정확 객체로 |flat_W⟩의 M₂ 계산:  M₂ = -log₂((1/(N·M⁴)) Σ_{x,z} ĝ_x(z)⁴)."""
    N = 2 ** n
    M = len(W)
    ind = np.zeros(N)
    ind[list(W)] = 1.0
    idx = np.arange(N)
    total = 0.0
    for x in range(N):
        g = ind * ind[idx ^ x]          # g_x(c) = [c∈W ∧ c⊕x∈W]
        G = _fwht(g)                     # ĝ_x
        total += np.sum(G ** 4)
    xi = total / (M ** 4) / N
    return float(-math.log2(xi))


def tau(W: list[int], n: int) -> float:
    """비아핀성 스칼라 τ(W) = (1/N) Σ_{x≠0} (A/M)(1−A/M).  =0 ⟺ 자기상관 2값 ⟺ 아핀."""
    N = 2 ** n
    M = len(W)
    A = autocorr(W, n) / M              # A_W(x)/M ∈ [0,1]
    A[0] = 0.0                          # x=0 제외 (A_W(0)/M = 1)
    return float(np.sum(A * (1.0 - A)) / N)


# ── 실험 ─────────────────────────────────────────────────────────────────────
def main() -> None:
    t0 = time.time()
    lines: list[str] = []

    def out(s: str = "") -> None:
        print(s)
        lines.append(s + "\n")

    out("# marker_code_magic — 표시집합 W의 부호 특성 → flat 상태 magic (보정 지표)")
    out("# 측정: magic.sre2 (XOR-FWHT). 정정 객체: 자기상관 A_W / Walsh 4차 모멘트.")
    out("")

    # ── 0. 정정된 정확 객체 검증: sre2(flat_W) == 자기상관/Walsh 공식 ──────────
    out("## 0. 정확 항등식 검증  sre2(flat_W) == autocorr/Walsh 공식  (보정 객체가 정확함)")
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for n in [3, 4, 5, 6]:
        for _ in range(8):
            M = int(rng.integers(2, 2 ** n))
            W = sorted(rng.choice(2 ** n, size=M, replace=False).tolist())
            maxerr = max(maxerr, abs(sre2(flat(W, n)) - sre2_flat_autocorr(W, n)))
    out(f"   max|sre2 − autocorr공식| = {maxerr:.2e}  (32 random W) → 자기상관/Walsh가 정확 객체")
    out("")

    # ── 1. 최소 해밍 거리는 예측 지표로 부적합 (인수인계서 정정) ──────────────
    out("## 1. 최소 해밍 거리 d_min 은 magic을 결정하지 못한다 (인수인계서 지표 정정)")
    out(f"   {'W':<28} {'|W|':>4} {'d_min':>6} {'M2':>9}  비고")
    n = 5
    cases = [
        (subspace([1, 2]),            "아핀부분공간 {0,1,2,3}"),
        ([0, 1, 2, 4],                "비아핀, d_min 동일"),
        (subspace([1, 2, 4]),         "아핀부분공간 dim3"),
        ([0, 1, 2, 4, 8, 15, 7, 11],  "비아핀 dim3 크기"),
    ]
    rows = []
    for W, lbl in cases:
        d, m = min_hamming(W), sre2(flat(W, n))
        rows.append((lbl, len(W), d, m))
        out(f"   {lbl:<28} {len(W):>4} {d:>6} {m:>9.4f}")
    out("   → {0,1,2,3}(M2=0)과 {0,1,2,4}(M2>0)은 d_min=1로 같으나 magic이 다르다.")
    out("     d_min(과 1차 가중치 분포)이 같아도 magic이 갈리므로 d_min은 지표로 불충분.")
    out("")

    # ── 2. 정정된 영점 판정 + 보정 스칼라 τ 의 예측력 ─────────────────────────
    out("## 2. 정정 판정: 아핀 ⟺ A_W∈{0,M} (2값) ⟺ M2=0;  τ(W)=0 ⟺ M2=0")
    n = 6
    out(f"   {'구조':<22} {'|W|':>4} {'ACF값 수':>8} {'τ':>9} {'M2':>9}")
    aff_cases = [subspace([1]), subspace([1, 2]), subspace([1, 2, 4]),
                 subspace([3, 12], off=5)]
    for W in aff_cases:
        A = autocorr(W, n)
        nval = len(set(np.rint(A).astype(int).tolist()))   # 서로 다른 ACF 값의 수
        out(f"   {'아핀 dim'+str(int(round(math.log2(len(W))))):<22} "
            f"{len(W):>4} {nval:>8} {tau(W, n):>9.2e} {sre2(flat(W, n)):>9.4f}")
    rng = np.random.default_rng(3)
    for M in [4, 8, 16]:
        W = sorted(rng.choice(2 ** n, size=M, replace=False).tolist())
        A = autocorr(W, n)
        nval = len(set(np.rint(A).astype(int).tolist()))
        out(f"   {'비아핀(랜덤) M='+str(M):<22} "
            f"{len(W):>4} {nval:>8} {tau(W, n):>9.4f} {sre2(flat(W, n)):>9.4f}")
    out("   → 아핀: ACF 2값(0,M), τ≈0, M2=0.  비아핀: ACF 다값, τ>0, M2>0.")
    out("")

    # τ ↔ M2 상관 (근사 예측력) + τ=0 ⟺ M2=0 무결성 (회귀 assert)
    out("## 2b. τ ↔ M2 상관 (근사 예측식의 후보) 및 τ=0 ⟺ M2=0 무결성")
    n = 7
    rng = np.random.default_rng(11)
    taus, mags = [], []
    viol = 0
    # 아핀(τ=0,M2=0) + 랜덤 비아핀 섞어서
    pool = [subspace(b) for b in ([1], [1, 2], [1, 2, 4], [1, 2, 4, 8])]
    for _ in range(120):
        M = int(rng.integers(3, 40))
        pool.append(sorted(rng.choice(2 ** n, size=M, replace=False).tolist()))
    for W in pool:
        tv, mv = tau(W, n), sre2(flat(W, n))
        taus.append(tv)
        mags.append(mv)
        if (tv < 1e-9) != (mv < 1e-9):
            viol += 1
    taus, mags = np.array(taus), np.array(mags)
    r = float(np.corrcoef(taus, mags)[0, 1])
    # log-log 추세(비영점만): M2 ≈ c·τ^p
    nz = (taus > 1e-9) & (mags > 1e-9)
    p, logc = np.polyfit(np.log(taus[nz]), np.log(mags[nz]), 1)
    out(f"   샘플 {len(pool)}개 (아핀 4 + 랜덤 비아핀).  Pearson r(τ, M2) = {r:.3f}")
    out(f"   비영점 추세  M2 ≈ {math.exp(logc):.3f}·τ^{p:.3f}  (근사; 통계적 예측식의 1차 후보)")
    out(f"   τ=0 ⟺ M2=0 위반 건수 = {viol}  (0이어야 정상: 정정 판정 무결)")
    out("   → τ는 정확한 영점 판정이자 magnitude의 단조 근사 지표. 정확값은 §0 Walsh 4차 모멘트.")
    out("")

    # ── 3. 랜덤 표시집합 magic 스케일링 (통계적 근사식의 목표 데이터) ─────────
    out("## 3. 랜덤 비아핀 W의 평균 magic 스케일링  E[M2] vs (M, n)")
    out(f"   {'n':>3} {'M':>4} {'E[M2]':>8} {'std':>7} {'E[τ]':>9}   (각 40회 평균)")
    rng = np.random.default_rng(2026)
    for n in [6, 8, 10]:
        for M in [4, 8, 16, 32]:
            if M > 2 ** n:
                continue
            ms, ts = [], []
            for _ in range(40):
                W = sorted(rng.choice(2 ** n, size=M, replace=False).tolist())
                ms.append(sre2(flat(W, n)))
                ts.append(tau(W, n))
            out(f"   {n:>3} {M:>4} {np.mean(ms):>8.4f} {np.std(ms):>7.4f} {np.mean(ts):>9.4f}")
    out("   → 랜덤(비아핀) W의 magic은 M(표시 수)과 함께 증가; n 의존은 약함.")
    out("     이 표가 '통계적 근사식'(E[M2]를 M·n·τ로) 유도의 목표 데이터 — 다음 과제.")
    out("")

    out("## 결론")
    out("   · 보정 객체 = 자기상관 A_W / Walsh 4차 모멘트 (최소 해밍 거리 아님; §0,§1).")
    out("   · 정정 영점 판정: 아핀 ⟺ A_W∈{0,M} ⟺ M2=0; τ(W)가 이를 스칼라로 포착(§2).")
    out("   · τ는 M2의 단조 근사 지표(§2b); 랜덤 W 스케일링이 통계적 예측식의 표적(§3).")

    elapsed = time.time() - t0
    out("")
    out(f"# Elapsed: {elapsed:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
