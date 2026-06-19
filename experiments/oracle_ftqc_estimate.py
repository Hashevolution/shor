"""
oracle_ftqc_estimate.py — M4/실용: 명제 6을 실제 FTQC 자원(T-count)으로 환산.

명제 6: M₂(|ψ_f⟩)>0 ⟺ f 비선형(차수≥2 ANF) ⟺ U_f 가 비클리포드(T) 게이트 요구.
여기서는 그 "비선형성"을 *실제 게이트 합성 모델*로 T-count까지 환산해, 상태 magic M₂ 와
나란히 보여준다. (T-count 는 결함허용 양자컴퓨터(FTQC)의 지배적 비용 — 각 T 는 매직상태
증류 1회를 요구.)

합성 모델 (XOR 오라클 U_f:|x⟩|y⟩→|x⟩|y⊕f(x)⟩, 출력비트별 ANF):
  • 차수 0 (상수)   → X            (클리포드, T=0)
  • 차수 1 (x_i)    → CNOT         (클리포드, T=0)
  • 차수 d≥2 단항식 → C^d(X) 다중제어-X = (2d−3) Toffoli   (Barenco et al., ancilla 사용)
      Toffoli 1개 = 7 T  (Amy–Maslov–Mosca–Roetteler, 정확 Clifford+T)
                  = 4 T  (Jones 2013, ancilla+측정)
  ⟹ T_est = Σ_{출력비트} Σ_{비선형 단항식 (차수 d)} (2d−3)·{7 또는 4}.

정직한 한계 (반드시 명시):
  이 *출력비트별 ANF* 합성은 일반적으로 최적이 아니며(부분항 공유 무시), 특히 *산술* 오라클
  (modexp)에는 매우 비효율적이다 — 실제 Shor 는 windowed 모듈러 산술로 T-count 가 훨씬 낮다
  (Gidney–Ekerå 2021: RSA-2048 ≈ 2.7×10⁹ Toffoli, ~20M 물리큐비트). 따라서 T_est 는
  *상한/지표*이지 생산용 추정이 아니다. 실용 가치는 (i) 정확한 영점판정(아핀⟺T=0),
  (ii) 정성적 magic↔T 법칙, (iii) *구조 없는* 오라클의 빠른 상한 지표.

Reproduction:
  python -u -m experiments.oracle_ftqc_estimate
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("oracle_ftqc_estimate_results.txt")

T_PER_TOFFOLI = 7   # Amy–Maslov–Mosca–Roetteler (exact Clifford+T)
T_PER_TOFFOLI_ANC = 4  # Jones 2013 (ancilla + measurement)


def graph_state(values, n_in, n_out):
    psi = np.zeros(2 ** (n_in + n_out))
    for x in range(2 ** n_in):
        psi[(x << n_out) | values[x]] = 1.0
    return psi / np.sqrt(2 ** n_in)


def anf_coeffs(truth, n):
    a = truth.astype(np.uint8).copy()
    for i in range(n):
        bit = 1 << i
        idx = np.arange(2 ** n)
        sel = (idx & bit) != 0
        a[sel] ^= a[idx[sel] ^ bit]
    return a


def synth_costs(values, n_in, n_out):
    """ANF 합성 모델로 (비선형단항식수 T_proxy, Toffoli_est, 최고차수)."""
    vals = np.asarray(values)
    deg = np.array([bin(s).count("1") for s in range(2 ** n_in)])
    nnl, toff, dmax = 0, 0, 0
    for k in range(n_out):
        truth = ((vals >> k) & 1).astype(np.uint8)
        coeff = anf_coeffs(truth, n_in)
        nz = np.where(coeff != 0)[0]
        for s in nz:
            d = int(deg[s])
            if d >= 2:
                nnl += 1
                toff += (2 * d - 3)
                dmax = max(dmax, d)
    return nnl, toff, dmax


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


def main():
    t0 = time.time()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s + "\n")

    out("# oracle_ftqc_estimate — 명제 6의 실제 FTQC 환산: 비선형 ANF → Toffoli → T-count")
    out("# 모델: 차수 d≥2 단항식 = (2d−3) Toffoli;  Toffoli = 7 T (또는 ancilla 4 T).")
    out("# 한계: 출력비트별 ANF 합성 = 상한/지표 (산술 오라클엔 비효율 — 아래 정직고지).")
    out("")

    def row(label, vals, n_in, n_out, compute_m2=True):
        nnl, toff, dmax = synth_costs(vals, n_in, n_out)
        t7, t4 = toff * T_PER_TOFFOLI, toff * T_PER_TOFFOLI_ANC
        if compute_m2 and (n_in + n_out) <= 11:
            m2 = f"{sre2(graph_state(vals, n_in, n_out)):.4f}"
        else:
            m2 = "  (skip)"
        out(f"   {label:<20} {n_in:>3} {n_out:>4} {nnl:>7} {dmax:>5} {toff:>8} {t7:>8} {t4:>8} {m2:>9}")

    # 1) 영점판정 + 동반 증가
    out("## 1. 영점판정·동반증가  (T_est=0 ⟺ M₂=0 ⟺ f 아핀)")
    out(f"   {'oracle':<20} {'in':>3} {'out':>4} {'#비선형':>7} {'dmax':>5} {'Toffoli':>8} "
        f"{'T(7)':>8} {'T(4)':>8} {'M2':>9}")
    row("Simon 선형 n=4", simon_linear(4, 3), 4, 4)
    base = simon_linear(4, 3)
    row("AND항 1개", add_and_terms(base, 4, [(0, 1, 2)]), 4, 4)
    row("AND항 2개", add_and_terms(base, 4, [(0, 1, 2), (0, 2, 3)]), 4, 4)
    row("AND항 3개", add_and_terms(base, 4, [(0, 1, 2), (0, 2, 3), (1, 3, 0)]), 4, 4)
    out("   → T_est 와 M₂ 가 함께 0에서 켜진다 (명제 6의 FTQC 환산).")
    out("")

    # 2) Shor modexp 오라클 (실제 알고리즘 오라클)
    out("## 2. Shor modexp 오라클 f(x)=aˣ mod N  (ANF-합성 상한)")
    out(f"   {'oracle':<20} {'in':>3} {'out':>4} {'#비선형':>7} {'dmax':>5} {'Toffoli':>8} "
        f"{'T(7)':>8} {'T(4)':>8} {'M2':>9}")
    for (N, a, t) in [(15, 7, 5), (15, 2, 5), (21, 2, 6), (33, 5, 7), (35, 3, 8)]:
        nout = N.bit_length()
        vals = [pow(a, x, N) for x in range(2 ** t)]
        row(f"a={a} mod {N}", vals, t, nout)
    out("   → 비선형성(따라서 T_est)이 N과 함께 급증; M₂(계산 가능 범위)도 동반 상승.")
    out("")

    # 3) FTQC 물리자원 환산 (차수 높은 예시 하나)
    out("## 3. T-count → FTQC 비용 (지배항)")
    out("   각 T ≈ 매직상태 증류 1회. 표면부호 FTQC에서 T-count 가 시간·증류공장 비용을 지배.")
    out("   예: 위 'a=5 mod 33' (ANF-합성 상한):")
    N, a, t = 33, 5, 7
    nout = N.bit_length()
    vals = [pow(a, x, N) for x in range(2 ** t)]
    nnl, toff, dmax = synth_costs(vals, t, nout)
    out(f"     Toffoli≈{toff},  T≈{toff*7} (7T) / {toff*4} (4T),  최고차수 {dmax}")
    out(f"     ⟹ 매직상태 ≈ T-count 개. (논리큐비트·코드거리는 목표 오류율에 따라 별도 산정.)")
    out("")

    # 4) 정직고지
    out("## 4. 정직한 한계 (반드시 함께 읽을 것)")
    out("   • 출력비트별 ANF 합성은 *상한*: 부분항 공유·산술 구조를 안 쓴다.")
    out("   • 특히 modexp 는 windowed 모듈러 산술이 훨씬 효율적 — 실제 Shor T-count 는 위 값보다")
    out("     수십~수천배 낮다 (Gidney–Ekerå 2021: RSA-2048 ≈ 2.7×10⁹ Toffoli, ~20M 물리큐비트).")
    out("   • 따라서 본 추정의 실용 가치 = (i) 정확한 영점판정(아핀⟺T=0), (ii) 정성적 magic↔T 법칙,")
    out("     (iii) *구조 없는/비산술* 오라클의 빠른 상한 지표. 생산용 추정기(Azure QRE 등) 대체 아님.")

    out("")
    out(f"# Elapsed: {time.time() - t0:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
