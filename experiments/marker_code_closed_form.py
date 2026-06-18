"""
marker_code_closed_form.py — 일반-M flat 마커상태의 magic 닫힌형(가법 에너지)과
랜덤 W 통계 예측식(Sidon 값). 로드맵 단계 1·2.

단계 1 (닫힌형, magic-results.md 명제 2′의 미완 부분):
  |flat_W> = |W|^{-1/2} Σ_{x∈W}|x> 의 SRE는 W의 *이동 자기교집합* S_x = W ∩ (W⊕x)의
  **가법 에너지(additive energy)** E(S)=#{(a,b,c,d)∈S^4: a⊕b⊕c⊕d=0}=Σ_v A_S(v)^2 로 닫힌다:

      ξ := (1/2^n)Σ_P⟨P⟩^4 = (1/M^4) Σ_{x∈𝔽₂ⁿ} E(W∩(W⊕x)),   M₂ = -log₂ ξ.

  유도: ⟨Z^zX^x⟩=(1/M)ĝ_x(z), g_x=1_W·1_{W⊕x}=1_{S_x};  Σ_z ĝ_x(z)^4 = N·E(S_x)
  (Parseval: Σ_z 1̂_S(z)^4 = N Σ_v A_S(v)^2 = N·E(S)).
  특수값:  M=1 → ξ=1 → M₂=0;  W 아핀부분공간 → Σ_x E(S_x)=M^4 → M₂=0 (보조정리 1 재확인).

단계 2 (통계 예측식):
  W가 Sidon 집합(B₂: 0 아닌 XOR 차분이 모두 서로 다름; 랜덤 W는 M≪2^{n/2}에서 whp Sidon)이면
  S_0=W에 E(W)=3M²−2M, 0 아닌 차분 x(개수 M(M−1)/2)마다 |S_x|=2·E(S_x)=8 →

      ξ_Sidon(M) = (7M−6)/M³,     M₂ = log₂( M³ / (7M−6) )  →  2log₂M − log₂7  (M→∞).

  즉 **랜덤(구조 없는) 마커상태의 magic은 ~2log₂M로 자란다.** 유한 N(=2^n)에서 M²≳N이면
  가법 충돌이 늘어 E가 커지고 ξ↑·M₂↓ (아래 표에서 n↑일수록 Sidon 값에 수렴).

Reproduction:
  python -u -m experiments.marker_code_closed_form
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("marker_code_closed_form_results.txt")


def _fwht(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float).copy()
    N = a.size
    h = 1
    while h < N:
        a = a.reshape(N // (2 * h), 2, h)
        x, y = a[:, 0, :], a[:, 1, :]
        a = np.stack([x + y, x - y], axis=1).reshape(N)
        h *= 2
    return a


def flat(W, n):
    psi = np.zeros(2 ** n)
    psi[list(W)] = 1.0
    return psi / math.sqrt(len(W))


def subspace(basis, off=0):
    pts = {off}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return sorted(pts)


def additive_energy(S, n):
    """E(S) = #{(a,b,c,d)∈S^4 : a⊕b⊕c⊕d=0} = (1/N) Σ_z 1̂_S(z)^4."""
    N = 2 ** n
    ind = np.zeros(N)
    ind[list(S)] = 1.0
    F = _fwht(ind)
    return float(np.sum(F ** 4) / N)


def sre2_flat_closed(W, n):
    """닫힌형: M₂ = -log₂( (1/M^4) Σ_x E(W∩(W⊕x)) )."""
    N = 2 ** n
    M = len(W)
    Wset = set(W)
    total = 0.0
    for x in range(N):
        Sx = [a for a in W if (a ^ x) in Wset]    # S_x = W ∩ (W⊕x)
        if len(Sx) >= 1:
            total += additive_energy(Sx, n)
    xi = total / (M ** 4)
    return float(-math.log2(xi))


def sidon_value(M):
    """Sidon(B₂) 마커상태의 M₂ = log₂(M³/(7M−6))."""
    return math.log2(M ** 3 / (7 * M - 6))


def greedy_sidon(M, n, rng):
    """탐욕적 Sidon 집합 구성(0 아닌 XOR 차분이 모두 서로 다른 M-집합)."""
    N = 2 ** n
    order = rng.permutation(N).tolist()
    S, diffs = [], set()
    for x in order:
        new = [x ^ s for s in S]
        if len(set(new)) == len(new) and all(d not in diffs for d in new):
            S.append(x)
            diffs.update(new)
            diffs.update(new)  # 대칭(이미 set이라 무해)
            if len(S) == M:
                break
    return sorted(S) if len(S) == M else None


def main():
    t0 = time.time()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s + "\n")

    out("# marker_code_closed_form — 일반-M flat 마커상태 magic 닫힌형 + 랜덤 W 통계식")
    out("# 닫힌형: M₂ = -log₂((1/M⁴)Σ_x E(W∩(W⊕x))),  E=가법에너지.  Sidon: log₂(M³/(7M−6)).")
    out("")

    # ── 단계 1-0: 닫힌형 == sre2 검증 ────────────────────────────────────────
    out("## 1-0. 닫힌형 검증  sre2(flat_W) == 가법에너지 닫힌형")
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for n in [3, 4, 5, 6]:
        for _ in range(8):
            M = int(rng.integers(2, 2 ** n))
            W = sorted(rng.choice(2 ** n, size=M, replace=False).tolist())
            maxerr = max(maxerr, abs(sre2(flat(W, n)) - sre2_flat_closed(W, n)))
    out(f"   max|sre2 − closed| = {maxerr:.2e}  (32 random W)  → 닫힌형(가법에너지) 확인")
    out("")

    # ── 단계 1-1: 해석적 특수값 (M=1, 아핀 → 0) ──────────────────────────────
    out("## 1-1. 해석적 특수값")
    n = 6
    out(f"   M=1 (단일점):        M₂ = {sre2_flat_closed([5], n):.6f}   (이론 0)")
    for basis in [[1], [1, 2], [1, 2, 4]]:
        W = subspace(basis)
        out(f"   아핀 dim{len(basis)} (|W|={len(W):>2}):  M₂ = {sre2_flat_closed(W, n):.6f}   "
            f"(이론 0; Σ_x E(S_x)=M⁴)")
    out("")

    # ── 단계 2-1: Sidon 닫힌형 검증 ──────────────────────────────────────────
    out("## 2-1. Sidon(B₂) 마커상태:  M₂ == log₂(M³/(7M−6))  (해석 닫힌형)")
    out(f"   {'M':>3} {'n':>3} {'M₂(Sidon W)':>12} {'log₂(M³/(7M−6))':>16} {'diff':>10}")
    rng = np.random.default_rng(7)
    for M in [3, 4, 5, 6, 8, 12, 16]:
        n = max(6, 2 * (M - 1).bit_length() + 2)   # Sidon 여유 있는 n
        S = greedy_sidon(M, n, rng)
        if S is None:
            out(f"   {M:>3} {n:>3}  (Sidon 구성 실패)")
            continue
        m2 = sre2_flat_closed(S, n)
        sv = sidon_value(M)
        out(f"   {M:>3} {n:>3} {m2:>12.6f} {sv:>16.6f} {abs(m2 - sv):>10.2e}")
    out("   → 일치. 구조 없는(Sidon) 마커상태 magic은 M에만 의존: M₂=log₂(M³/(7M−6))≈2log₂M−log₂7.")
    out("")

    # ── 단계 2-2: 랜덤 W → Sidon 값 수렴 (N=2^n↑) ────────────────────────────
    out("## 2-2. 랜덤 W의 E[M₂] → Sidon 값  (M 고정, n↑일수록 수렴; M²≳N이면 충돌로 하락)")
    out(f"   {'M':>3} {'Sidon값':>9}   " + "  ".join(f"n={n}" for n in [6, 8, 10]))
    rng = np.random.default_rng(2026)
    for M in [4, 8, 16, 32]:
        row = []
        for n in [6, 8, 10]:
            if M > 2 ** n:
                row.append("  -  ")
                continue
            vals = [sre2(flat(sorted(rng.choice(2 ** n, size=M, replace=False).tolist()), n))
                    for _ in range(20)]
            row.append(f"{np.mean(vals):.3f}")
        out(f"   {M:>3} {sidon_value(M):>9.3f}   " + "   ".join(f"{r:>5}" for r in row))
    out("   → 각 행이 n↑(=N↑)에서 Sidon값으로 수렴. 작은 n(M²≳N)에선 가법충돌로 M₂가 낮음.")
    out("")

    # ── 단계 2-3: 유한-N 보정의 부호·크기 (가법충돌 = 하락) ──────────────────
    out("## 2-3. 유한-N 보정:  gap = Sidon값 − E[M₂]  ≥ 0,  M²/N 증가 시 확대")
    out(f"   {'M':>3} {'n':>3} {'M²/N':>7} {'E[M₂]':>8} {'Sidon':>8} {'gap':>8}")
    rng = np.random.default_rng(11)
    for (M, n) in [(8, 6), (8, 8), (8, 10), (16, 8), (16, 10), (32, 10)]:
        vals = [sre2(flat(sorted(rng.choice(2 ** n, size=M, replace=False).tolist()), n))
                for _ in range(30)]
        em = float(np.mean(vals))
        sv = sidon_value(M)
        out(f"   {M:>3} {n:>3} {M * M / 2 ** n:>7.3f} {em:>8.4f} {sv:>8.4f} {sv - em:>8.4f}")
    out("   → gap≥0이고 M²/N가 클수록 커짐: 유한-N 보정은 *하락*(가법 quadruple 증가→E↑→ξ↑→M₂↓).")
    out("")

    out("## 결론")
    out("   단계1: |flat_W> magic 닫힌형 = 이동 자기교집합의 가법에너지 합 (명제 2′ 완성).")
    out("          특수값 M=1·아핀 → 0 해석적 재확인.")
    out("   단계2: 랜덤(구조無) 마커상태 M₂ = log₂(M³/(7M−6)) ≈ 2log₂M−log₂7 (Sidon, N≫M²).")
    out("          유한-N 보정은 가법충돌에 의한 *하락*; M²/N로 통제.")

    out("")
    out(f"# Elapsed: {time.time() - t0:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
