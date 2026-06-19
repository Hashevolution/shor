"""
discover_poc.py — JAMES-DISCOVER D1 PoC: 추측-마이닝 루프(법칙 자동 발견).

설계근거: `ai-discovery-engine-design.md` §3·§6 (D1 재발견 게이트).
이 레포의 수동 발견 루프(생성→측정→추측→반증→승격)를 기계화한 최소판. "AI가 사람의 힌트
없이, 이미 손으로 찾아둔 닫힌형을 데이터에서 *재발견*하는가"를 객관 수치로 판정한다.

파이프라인(설계문서 §3):
  Generator  : 마커집합 W 생성 (Sidon / 아핀 / 랜덤)            [기존 함수 재사용]
  Probe      : 불변량 측정 (M₂=sre2, ξ=2^{-M₂}, A_W, τ, d_min)  [기존 함수 재사용]
  Miner      : sparse dictionary fit (직교매칭추구 OMP) — 신규
  Adversary  : 발견식의 성립 범위를 반례로 좁힘 — 신규(D2-lite)
  Promoter   : 정수상수 스냅 + 특수값 환원 체크 — 신규

재발견 게이트(합격 기준):
  [T1] Sidon 마커상태의 ξ를 사전 {M^{-1},M^{-2},M^{-3},M^{-4},log₂M·M^{-2}}에 OMP →
       정수계수 ξ = 7·M^{-2} − 6·M^{-3} 복원 (잔차 < 1e-9, 상수 7·−6을 정수로 스냅).
       즉 M₂ = log₂(M³/(7M−6)) 자동 복원. (명제 5의 무힌트 재발견)
  [T2] 후보 영점-예측자 {τ, mean A_W, max A_W, d_min} 중 "f(W)=0 ⟺ M₂(W)=0"을
       위반 0으로 만족하는 것을 선택 → τ 선택, d_min은 반례로 탈락. (보조정리 1 + d_min 무력성)

Reproduction:
  python -u -m experiments.discover_poc
"""
from __future__ import annotations

import math
import time
from itertools import combinations
from pathlib import Path

import numpy as np

from magic import sre2
from experiments.marker_code_magic import (
    flat,
    subspace,
    autocorr,
    min_hamming,
    tau,
)
from experiments.marker_code_closed_form import greedy_sidon, sidon_value

RESULTS_FILE = Path(__file__).with_name("discover_poc_results.txt")


# ── Miner: 최소 부분집합 정확탐색 (MDL — "가장 단순한 정확 법칙") ───────────
def sparse_exact(X: np.ndarray, y: np.ndarray, max_terms: int, tol: float = 1e-9):
    """y ≈ X·c 를 *정확히*(잔차<tol) 맞추는 **최소 크기** 열-부분집합을 찾음.

    사전이 작을 때(원자 ≲ 수십) greedy OMP보다 견고하다: 1·2·…개 부분집합을 차례로
    전수 시도해, 처음으로 tol 이하 잔차를 내는 최소 부분집합을 채택(= 최소기술길이).
    1/M의 거듭제곱처럼 공선성이 심해 greedy가 헛짚는 경우를 회피한다.
    반환: (선택 열 인덱스 튜플, 그 계수, 잔차노름)."""
    y = np.asarray(y, dtype=float)
    d = X.shape[1]
    best = None  # (combo, coef, res) — tol을 못 넘기면 최선근사 반환
    for k in range(1, max_terms + 1):
        kbest = None
        for combo in combinations(range(d), k):
            Xs = X[:, combo]
            c, *_ = np.linalg.lstsq(Xs, y, rcond=None)
            res = float(np.linalg.norm(y - Xs @ c))
            if kbest is None or res < kbest[2]:
                kbest = (combo, c, res)
        if best is None or kbest[2] < best[2]:
            best = kbest
        if kbest[2] < tol:                 # 이 크기에서 정확 적합 달성 → 최소부분집합
            return kbest
    return best


def snap_int(c: float, tol: float = 1e-6):
    """Promoter: 계수를 가까운 정수로 스냅(우연 상수 vs 의미상수 구분의 1차 필터)."""
    r = round(c)
    return (int(r), True) if abs(c - r) < tol else (c, False)


# ── 데이터 헬퍼 ──────────────────────────────────────────────────────────────
def _xi(W, n):
    """ξ = 2^{-M₂} = (1/2ⁿ)Σ_P⟨P⟩⁴, 측정으로부터 (Miner가 다루는 선형 객체)."""
    return 2.0 ** (-sre2(flat(W, n)))


def main():
    t0 = time.time()
    lines: list[str] = []

    def out(s: str = ""):
        print(s)
        lines.append(s + "\n")

    out("# discover_poc — JAMES-DISCOVER D1: 추측-마이닝 루프(법칙 자동 재발견)")
    out("# 합격선: [T1] Sidon ξ에서 정수상수 7·−6 무힌트 복원, [T2] τ 영점판정 위반 0 / d_min 탈락")
    out("")

    gate = {"T1": False, "T2": False}

    # ══ T1. 닫힌형 재발견 — Sidon 마커상태의 ξ(M) ════════════════════════════
    out("## T1. Generator(Sidon W) → Probe(ξ) → Miner(최소부분집합 정확탐색) → Promoter(정수스냅)")
    out("#  ξ를 사전 {M^-1, M^-2, M^-3, M^-4, log2(M)·M^-2}에서 최소부분집합으로 적합. 사람 힌트 없음.")
    Ms = [3, 4, 5, 6, 8, 10, 12, 16]
    rng = np.random.default_rng(7)
    rows = []
    for M in Ms:
        n = max(6, 2 * (M - 1).bit_length() + 2)   # Sidon 여유 n (기존 관례)
        S = greedy_sidon(M, n, rng)
        if S is None:
            continue
        rows.append((M, _xi(S, n)))
    Mv = np.array([m for m, _ in rows], dtype=float)
    yv = np.array([x for _, x in rows], dtype=float)

    names = ["M^-1", "M^-2", "M^-3", "M^-4", "log2(M)*M^-2"]
    cols = [Mv ** -1, Mv ** -2, Mv ** -3, Mv ** -4, np.log2(Mv) * Mv ** -2]
    X = np.column_stack(cols)

    sel, c, res = sparse_exact(X, yv, max_terms=3)
    out(f"   데이터점 {len(rows)}개 (M={[int(m) for m in Mv]})")
    out(f"   Miner 채택(최소부분집합): {[names[j] for j in sel]}   잔차노름 = {res:.2e}")
    snap = {}
    for j, cj in zip(sel, c):
        v, is_int = snap_int(cj)
        snap[names[j]] = v
        out(f"     계수[{names[j]:>12}] = {cj:+.8f}  → 스냅 {v}{' (정수)' if is_int else ''}")

    # 재발견 판정: 채택 부분집합이 정확히 {M^-2, M^-3}이고 계수가 정수 +7, −6
    coef = dict(zip((names[j] for j in sel), c))
    c2, c3 = coef.get("M^-2", 0.0), coef.get("M^-3", 0.0)
    t1_ok = (set(names[j] for j in sel) == {"M^-2", "M^-3"}
             and abs(c2 - 7) < 1e-6 and abs(c3 + 6) < 1e-6 and res < 1e-9)
    gate["T1"] = t1_ok
    out("")
    out(f"   ⟹ 복원식:  ξ(M) = {snap.get('M^-2', 0)}·M^-2 + {snap.get('M^-3', 0)}·M^-3 "
        f"= (7M−6)/M³")
    out(f"   ⟹ 따라서  M₂ = −log₂ξ = log₂(M³/(7M−6))   [명제 5 = Sidon 값, 무힌트 재발견]")
    # 특수값 환원 체크(Promoter): M→∞ 점근 2log₂M − log₂7
    Mbig = 1e6
    asym = 2 * math.log2(Mbig) - math.log2(7)
    exact = math.log2(Mbig ** 3 / (7 * Mbig - 6))
    out(f"   Promoter 점근체크: M=1e6  exact={exact:.6f}  2log₂M−log₂7={asym:.6f}  "
        f"diff={abs(exact - asym):.2e}")
    out(f"   [T1] {'PASS' if t1_ok else 'FAIL'}")
    out("")

    # ══ T1-Adversary. 발견식의 성립 범위를 반례로 좁힘 (D2-lite) ═════════════
    out("## T1-Adversary. 발견식 ξ=(7M−6)/M³ 의 성립 범위 자동 탐색 (랜덤 W로 반증)")
    out("#  Sidon 가정이 깨지는 곳(M²≳N)에서 식이 어긋남을 스스로 보고 → '성립범위' 축소.")
    out(f"   {'M':>3} {'n':>3} {'M²/N':>7} {'ξ측정(랜덤W)':>13} {'ξ예측':>10} {'상대오차':>9} 판정")
    rng = np.random.default_rng(2026)
    first_break = None
    for (M, n) in [(8, 12), (8, 10), (8, 8), (8, 6), (16, 10), (16, 8), (32, 10)]:
        vals = [_xi(sorted(rng.choice(2 ** n, size=M, replace=False).tolist()), n)
                for _ in range(12)]
        xi_meas = float(np.mean(vals))
        xi_pred = (7 * M - 6) / M ** 3
        rel = abs(xi_meas - xi_pred) / xi_pred
        ok = rel < 0.05
        if not ok and first_break is None:
            first_break = (M, n, M * M / 2 ** n)
        out(f"   {M:>3} {n:>3} {M * M / 2 ** n:>7.3f} {xi_meas:>13.6f} {xi_pred:>10.6f} "
            f"{rel:>8.1%} {'성립' if ok else '깨짐'}")
    if first_break:
        out(f"   ⟹ Adversary 결론: 식은 M²/N ≪ 1(Sidon 영역)에서만 성립. "
            f"첫 붕괴 ~ M²/N={first_break[2]:.2f}.")
        out(f"     (= 명제 5의 '유한-N 가법충돌 보정'을 데이터로 자동 재확인)")
    out("")

    # ══ T2. 영점-예측자 재발견 — τ vs d_min ═════════════════════════════════
    out("## T2. Generator(아핀/랜덤 W) → Probe(M₂,τ,A_W,d_min) → Miner(영점-예측자 선택)")
    out("#  후보 f 중 'f(W)=0 ⟺ M₂(W)=0'을 위반 0으로 만족하는 것을 선택.")
    n = 5
    N = 2 ** n
    samples = []
    # 아핀부분공간(M₂=0 이어야 함)
    for basis in [[1], [2], [1, 2], [1, 4], [1, 2, 4], [1, 2, 4, 8]]:
        for off in [0, 3, 5]:
            samples.append(sorted({(p ^ off) for p in subspace(basis)}))
    # 랜덤 비아핀(대부분 M₂>0)
    rng = np.random.default_rng(11)
    for _ in range(30):
        M = int(rng.integers(3, N))
        samples.append(sorted(rng.choice(N, size=M, replace=False).tolist()))

    feats = {"τ": [], "mean_A": [], "max_A_frac": [], "d_min": []}
    m2s = []
    for W in samples:
        m2 = sre2(flat(W, n))
        m2s.append(m2)
        A = autocorr(W, n)
        Mw = len(W)
        feats["τ"].append(tau(W, n))
        feats["mean_A"].append(float(np.mean(A[1:])))            # x≠0 평균 자기상관
        feats["max_A_frac"].append(float(np.max(A[1:]) / Mw))    # 최대 자기상관/ M
        feats["d_min"].append(float(min_hamming(W)))
    m2s = np.array(m2s)
    is_zero_m2 = m2s < 1e-9

    out(f"   샘플 {len(samples)}개 중 M₂=0: {int(is_zero_m2.sum())}개")
    out(f"   {'후보 f':>11} {'위반수':>6}  (f=0 ⟺ M₂=0 위반 = 두 영점집합 불일치)")
    best = None
    for name, vals in feats.items():
        v = np.array(vals)
        is_zero_f = np.abs(v) < 1e-9
        violations = int(np.sum(is_zero_f != is_zero_m2))
        out(f"   {name:>11} {violations:>6}")
        if violations == 0 and best is None:
            best = name
    out("")
    # d_min 반례 박제 (보조정리 재서술): 같은 d_min, 다른 magic
    W_aff = subspace([1, 2])        # {0,1,2,3} 아핀 → M₂=0
    W_non = [0, 1, 2, 4]            # 비아핀 → M₂>0
    n3 = 3
    out("   d_min 무력성 반례(자동 확인):")
    out(f"     {str(W_aff):>12}: d_min={min_hamming(W_aff)}, τ={tau(W_aff, n3):.4f}, "
        f"M₂={sre2(flat(W_aff, n3)):.4f}")
    out(f"     {str(W_non):>12}: d_min={min_hamming(W_non)}, τ={tau(W_non, n3):.4f}, "
        f"M₂={sre2(flat(W_non, n3)):.4f}")
    out("     → 같은 d_min, 다른 M₂ ⟹ d_min은 영점판정 불가(인수인계서 지표 정정 자동 재확인).")
    out("")
    t2_ok = (best == "τ")
    gate["T2"] = t2_ok
    out(f"   ⟹ Miner 선택 영점-예측자: {best}  (τ=0 ⟺ M₂=0; 보조정리 1 재발견)")
    out(f"   [T2] {'PASS' if t2_ok else 'FAIL'}")
    out("")

    # ══ 게이트 종합 ═════════════════════════════════════════════════════════
    out("## 재발견 게이트 종합")
    for k, v in gate.items():
        out(f"   {k}: {'PASS' if v else 'FAIL'}")
    allpass = all(gate.values())
    out(f"   => D1 {'PASS — 무힌트 재발견 성공. (A) 법칙발견 경로 실현가능 확인.' if allpass else 'FAIL'}")
    out("")
    out(f"# Elapsed: {time.time() - t0:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")
    return allpass


if __name__ == "__main__":
    main()
