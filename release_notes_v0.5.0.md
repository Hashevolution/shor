# Release notes — v0.5.0 (2026-06-17)

**DOI**: [10.5281/zenodo.20725965](https://doi.org/10.5281/zenodo.20725965)

## Summary

**v0.5.0** 은 v0.4.0 (DOI
[10.5281/zenodo.20688069](https://doi.org/10.5281/zenodo.20688069)) 위에 올린
**새 연구 방향(magic = 비안정자성)** 릴리스다. 기존 order-finding 정리(T1–T6)는 그대로 두고,
**양자 속도우위의 *유형*이 비안정자성(magic)의 *양/밀도*와 어떻게 연결되는지**를 이 레포의
상태벡터 시뮬레이션 위에서 직접 측정·증명한 결과를 추가한다.

## What's new — magic across quantum speedups

### 핵심 결과 (사다리)

| 알고리즘 | 속도우위 | M₂ 거동 | 밀도 M₂/n |
|---|---|---|---|
| Simon / BV | 지수(쿼리) | 0 (아핀 오라클 = 클리포드) | 0 |
| **Grover** | 2차(다항) | **정점 → 3 bit, 답에서 0** | **→ 0** |
| Shor (comb) | 지수 | ∝ t 증가 | ~0.4–0.55 |
| Shor (in-circuit) | 지수 | → L (최대) | → 1 |

> 비안정자성의 *양/밀도*가 속도우위의 유형을 가른다. magic은 회로 표면이 아니라 **문제의
> 비선형성**에 살며, 오라클·FFT 같은 블랙박스가 그것을 숨길 수 있다.

### 명제 (magic-results.md)

- **보조정리 1.** 평탄(균등진폭·동일위상) 상태 `M₂=0 ⟺ 받침이 𝔽₂-아핀부분공간`.
- **따름정리 1.** 그래프상태 `(1/√N)Σ|x⟩|f(x)⟩` 의 `M₂=0 ⟺ f 아핀`. ⟹ Simon은 선형오라클로
  `M₂≡0`(지수 쿼리속도우위), Shor modexp는 비선형으로 magic 강제 — **오라클-가림** 정식화.
- **명제 2–3.** 단일표시 Grover의 닫힌형 + 정점이 `−log₂(a⁸+(1−a²)⁴)` 의 최댓값 = **정확히
  3 bit** (a²=1/2), 밀도 `M₂/n→0`.
- **명제 2′ (일반 M).** `|ψ⟩ = b√N|+⟩^n + (a−b)|W̃⟩`: 표시집합 W가 아핀이면 magic 유한(안정자
  2개 중첩), 비아핀이면 표시상태 자체 magic이 추가.
- **정정.** "comb magic 0 ⟺ 주기가 2의 거듭제곱"은 *거의*-참(정확 판정 = 받침의 아핀성).

## What changed

### 새 파일
- `magic.py` — stabilizer 2-Rényi 엔트로피 측정도구 (XOR-FWHT, O(n·4ⁿ) 정확).
- `experiments/grover_magic.py` — Grover magic 궤적·정점 포화·일반 M.
- `experiments/shor_comb_magic.py` — comb magic, 정정된 T1, 자기완결 Shor 대조.
- `experiments/oracle_magic.py` — 오라클-가림(아핀 ⟺ M₂=0).
- `experiments/magic_proofs_check.py` — 명제 회귀검증 (42 assert).
- `magic-results.md` (명제·증명), `magic-쉬운설명.md` (초등생 수준 설명),
  `magic-and-quantum-speedup.md` (연구노트), `magic-prior-art.md` (선행조사 + DOI 검토).

### 검증
- SRE 도구를 **3개 독립 구현**(XOR-FWHT ↔ Pauli순열 ↔ kron행렬)으로 교차검증 → 10⁻¹⁵ 일치.
- 상한 0 ≤ M₂ ≤ n 확인. 명제 회귀검사 **42 assert 전부 통과**.
- 정점 수렴: n=10 → 2.815, 20 → 2.994, 30 → 3.000.
- 재현: `python -u -m experiments.magic_proofs_check` (빠름), 그 외 각
  `python -u -m experiments.<name>`.

## 선행연구와의 선 (전문/저널측 검토)

- **arXiv:2605.05347** (Paviglianiti et al.): Shor magic ↔ 주기 — 선점 결과로 인용. **아직
  프리프린트**(peer-review 전)라 본 결과의 신규성 창은 열려 있음.
- 양자걷기 magic (출판된 최근접 이웃): **Phys. Rev. Research 10.1103/7rwg-lhpv** (2506.17783),
  **Phys. Rev. B 113, 075142 (2026) 10.1103/nzrp-49mr** (2504.19750) — 둘 다 1D 격자 수송이라
  Grover/완전그래프/탐색 미포함(전문 확인). Grover 3-bit 포화는 그 빈칸.
- 알고리즘 magic 비교 (2505.17185·2507.16543): 변분/QFT 회로, Grover·Simon·Shor 사다리 미포함.

## 범위/한계

- 명제 2–3 닫힌형은 단일표시(M=1); 일반 아핀 M은 구조적 분해(명제 2′)로 유한성 보장, 명시적
  닫힌형은 미완.
- comb magic은 측정후 한계분포의 양으로, in-circuit의 → L (2605.05347)과는 구분되는 다른 양.
- 이 릴리스는 기존 order-finding 결과(T1–T6)를 변경하지 않는다 — magic은 *추가* 방향이다.
