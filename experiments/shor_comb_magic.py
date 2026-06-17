"""
shor_comb_magic.py — Shor 측정후 comb 상태의 비안정자성(magic), Grover와의 대조.

동기 (magic-and-quantum-speedup.md, T3):
  Grover(2차 속도우위)의 magic은 정점 3 bit로 *포화*(밀도 M₂/n→0; experiments/grover_magic).
  그렇다면 Shor(지수 속도우위)는? — 같은 도구 magic.sre2 로 직접 계산해 대조를 자기완결화한다.
  (인용에만 의존하던 Shor 쪽을 레포 안에서 재현.)

상태 — `shor.py`가 작업 레지스터를 먼저 측정하면 계산 레지스터(t큐비트, Q=2ᵗ)는
주기적 comb로 붕괴한다:
    |ψ⟩ = (1/√m) Σ_{x ≡ x₀ (mod r)} |x⟩,   m = #{x∈[0,Q): x ≡ x₀ (mod r)}.

결과:
  1. 정정된 T1: r이 **2의 거듭제곱이면 M₂=0**(하위비트 고정 ⊗ 상위 |+⟩ = 안정자),
     **홀수 인수가 있으면 M₂>0.** (단 r·(m−1)=Q−1 류의 우연한 정렬은 |j⟩|j⟩ 꼴이 되어
     M₂=0 — 예: t=10, r=33. "M₂=0 ⟺ 2의 거듭제곱"은 거의-참이나 정확히는 아님.)
  2. comb의 magic은 **레지스터 크기 t와 함께 무한정 증가**(밀도 M₂/t ≈ 0.4–0.55 유지).
     → Shor의 magic은 *문제크기에 비례*. Grover(유한, 밀도→0)와 정반대.
  3. (회로 차원의 더 강한 결과 M₂→L, 밀도→1 은 in-circuit 상태에 대한 것: 2605.05347.
     comb은 한계분포를 품은 다른 양이라 밀도가 1엔 못 미친다 — magic-…md §6 참조.)

Reproduction:
  python -u -m experiments.shor_comb_magic
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from magic import sre2

RESULTS_FILE = Path(__file__).with_name("shor_comb_magic_results.txt")


def comb_state(t: int, r: int, x0: int = 0) -> tuple[np.ndarray, int]:
    """측정후 comb 상태 (t큐비트) 와 그 항 개수 m."""
    Q = 2 ** t
    members = [x for x in range(Q) if x % r == x0 % r]
    psi = np.zeros(Q)
    psi[members] = 1.0
    return psi / math.sqrt(len(members)), len(members)


def main() -> None:
    t0 = time.time()
    lines: list[str] = []

    def out(s: str = "") -> None:
        print(s)
        lines.append(s + "\n")

    out("# shor_comb_magic — 측정후 comb 상태의 M₂ (stabilizer 2-Rényi 엔트로피, bits)")
    out("# 측정도구: magic.sre2 (XOR-FWHT, brute-force 검증됨)")
    out("")

    # 1) 정정된 T1: 주기 구조에 따른 magic (고정 t)
    t = 10
    out(f"## 1. 정정된 T1 — comb magic vs 주기 r  (t={t}, Q={2**t})")
    out(f"   {'r':>4} {'m':>5} {'M2':>8} {'M2/t':>7}   비고")
    cases = [
        (2, "2¹  = 2의 거듭제곱"), (4, "2²"), (8, "2³"), (16, "2⁴"), (32, "2⁵"),
        (3, "홀수"), (5, "홀수"), (7, "홀수"), (9, "3²"),
        (6, "2·3"), (12, "4·3"), (24, "8·3"), (48, "16·3"),
        (15, "3·5"), (21, "3·7"),
        (33, "3·11 — 우연한 안정자(33·31=1023 ⟹ |j⟩|j⟩)"),
    ]
    for r, note in cases:
        psi, m = comb_state(t, r)
        m2 = sre2(psi)
        out(f"   {r:>4} {m:>5} {m2:>8.4f} {m2/t:>7.3f}   {note}")
    out("   → 2의 거듭제곱 ⟹ M₂=0,  홀수 인수 ⟹ M₂>0 (우연한 정렬 제외).")
    out("")

    # 2) 스케일링: comb magic 이 레지스터 크기 t 와 함께 자라는가?
    out("## 2. comb magic vs 레지스터 크기 t  (고정 홀수주기) — 문제크기에 비례?")
    for r in [3, 5, 7]:
        row = []
        for t in range(5, 15):
            psi, _ = comb_state(t, r)
            row.append(sre2(psi))
        dens = row[-1] / 14
        seg = " ".join(f"{v:4.2f}" for v in row)
        out(f"   r={r}:  M2(t=5..14) = {seg}   (밀도 M₂/t≈{dens:.2f}, 증가추세)")
    out("   → comb magic은 t와 함께 *무한정 증가* (밀도 ~0.4–0.55 유지).")
    out("")

    # 3) 속도우위 유형별 대조표
    out("## 3. 대조 — 비안정자성의 양/밀도 vs 속도우위 유형")
    out(f"   {'알고리즘':<22} {'속도우위':<10} {'magic 거동':<28} {'밀도 M₂/n':<10}")
    out(f"   {'-'*22} {'-'*10} {'-'*28} {'-'*10}")
    out(f"   {'Simon / BV':<22} {'(지수,쿼리)':<10} {'0 (클리포드)':<28} {'0':<10}")
    out(f"   {'Grover 탐색':<20} {'2차(다항)':<10} {'정점→3 bit, 답서 0':<26} {'→ 0':<10}")
    out(f"   {'Shor comb(측정후)':<18} {'지수':<10} {'∝ t (무한정 증가)':<26} {'~0.4–0.55':<10}")
    out(f"   {'Shor in-circuit':<20} {'지수':<10} {'→ L (최대; 2605.05347)':<26} {'→ 1':<10}")
    out("")
    out("   결론: 비안정자성의 *양/밀도*가 속도우위 유형을 가른다 —")
    out("   클리포드-자명(0) → 2차 속도우위(유한, 밀도→0) → 지수 속도우위(문제크기 비례→최대).")

    elapsed = time.time() - t0
    out("")
    out(f"# Elapsed: {elapsed:.1f}s")
    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
