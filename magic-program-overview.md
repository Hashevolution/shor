# PROJECT JAMES — 비안정자성(magic) 자원이론 통합본 + 로드맵

*양자 속도우위의 *유형*이 비안정자성(magic)의 *양/밀도*와 어떻게 연결되는가 — 이 레포의
상태벡터 시뮬레이션 위에서 측정·증명·부호이론화한 결과 전체와, 앞으로의 정량 과제를 한 곳에 정리.*

*기준일 2026-06-18. 세부는 `magic-results.md`(명제·증명), `magic-and-quantum-speedup.md`(연구노트),
`magic-prior-art.md`(선행조사), `experiments/`(재현 코드). 본 문서는 그 통합 인덱스 + 로드맵.*

---

## 0. 한 줄 / 프로그램 정체성

> **비안정자성의 *양/밀도*가 속도우위의 유형을 가른다.** magic은 회로 표면이 아니라 *문제의
> 비선형성*(아핀 vs 곱셈적)에 살며, 오라클·FFT 같은 블랙박스가 그것을 숨길 수 있다.
> Grover의 2차 속도우위는 *유한한*(정점 3 bit) magic을 썼다가 되돌리고, Shor의 지수 속도우위는
> magic을 문제크기와 함께 키운다(→ 최대).

핵심 측정량: $n$큐비트 순수상태의 stabilizer 2-Rényi 엔트로피
$$M_2(|\psi\rangle)=-\log_2\!\Big(\tfrac1{2^n}\sum_{P\in\mathcal P_n}\langle\psi|P|\psi\rangle^4\Big),\qquad
0\le M_2\le n,\quad M_2=0\iff |\psi\rangle\ \text{안정자}.$$

---

## Part I — 그동안 확립된 것 (DONE, v0.5.0)

### I.1 측정 인프라 — 신뢰구축 완료
- **`magic.py`**: `sre2`가 $\sum_P\langle P\rangle^4=\sum_{x,z}|\mathrm{WHT}[h_x](z)|^4$
  ($h_x(c)=\psi^*_{c\oplus x}\psi_c$) 항등식으로 **XOR-FWHT $O(n\,4^n)$ 정확** 계산.
- **3중 교차검증**(XOR-FWHT ↔ Pauli순열 ↔ kron행렬) → $10^{-15}$ 일치.
- **회귀 42 assert 전부 통과**(`experiments/magic_proofs_check.py`). 정점 수렴 n=10→2.815,
  20→2.994, 30→3.000.
- 같은 도구가 Grover/Simon/Shor comb/하이퍼그래프 등 어떤 상태벡터에도 그대로 적용됨.

### I.2 명제·정리 (`magic-results.md`)
- **보조정리 1.** 평탄(균등진폭·동일위상) 상태 $M_2=0\iff$ 받침이 $\mathbb F_2^n$의 **아핀부분공간**.
- **따름정리 1 (오라클-가림).** 그래프상태 $2^{-n/2}\sum_x|x\rangle|f(x)\rangle$의 $M_2=0\iff f$ 아핀.
  ⟹ **Simon**은 선형오라클로 $M_2\equiv0$인데도 지수 *쿼리* 우위(=magic은 쿼리우위에 불필요),
  **Shor** modexp는 곱셈적=비선형 ⟹ magic 강제. magic의 원천은 *문제의 비선형성*.
- **따름정리 2.** 측정후 comb의 $M_2=0\iff$ 잉여류가 아핀; "comb magic 0 ⟺ 주기 $2^k$"는 *거의*-참.
- **명제 2–3 (Grover M=1).** 닫힌형 $\sum_P\langle P\rangle^4$ 유도 + 정점 $\to-\log_2(a^8+(1-a^2)^4)$,
  $a^2=\tfrac12$에서 **정확히 3 bit**, 밀도 $M_2/n\to0$.
- **명제 2′ (일반 M).** $|\psi\rangle=b\sqrt N|{+}\rangle^{\otimes n}+(a-b)|\widetilde W\rangle$:
  표시집합 $W$ 아핀 ⟹ 안정자 2개 중첩(magic 유한), 비아핀 ⟹ 표시상태 자체 magic 추가.
  *일반 닫힌형은 $W$의 자기상관 $A_W(x)=|W\cap(W{\oplus}x)|$에 의존(미완).*

### I.3 속도우위 사다리
| 알고리즘 | 속도우위 | $M_2$ 거동 | 밀도 $M_2/n$ | 근거 |
|---|---|---|---|---|
| Simon / BV | 지수(쿼리) | $0$ (아핀 오라클=클리포드) | $0$ | 따름정리 1 |
| **Grover** | 2차(다항) | 정점 $\to 3$ bit, 답서 $0$ | $\to 0$ | 명제 2–3 |
| Shor (comb) | 지수 | $\propto t$ 증가 | $\sim0.4$–$0.55$ | 따름정리 2 |
| Shor (in-circuit) | 지수 | $\to L$ (최대) | $\to 1$ | 2605.05347 |

### I.4 부호이론 보정 (`experiments/marker_code_magic.py`) — 인수인계서 지표 정정
표시집합 $W$를 고전 부호로 보고 $|{\rm flat}_W\rangle$의 magic을 예측. **정확 객체 = 자기상관/Walsh 4차모멘트:**
$$M_2=-\log_2\!\Big(\tfrac1{N M^4}\sum_{x,z}\hat g_x(z)^4\Big),\quad g_x(c)=1_W(c)\,1_W(c{\oplus}x).$$
- **정정 영점판정(보조정리 1의 부호이론적 재서술):** $W$ 아핀 $\iff A_W(x)\in\{0,M\}\ \forall x$ (자기상관 2값).
- **인수인계서의 "최소 해밍 거리" 지표는 부정확:** $\{0,1,2,3\}$($M_2{=}0$)과 $\{0,1,2,4\}$($M_2{=}1.54$)는
  $d_{\min}{=}1$로 같으나 magic이 갈림 → $d_{\min}$ 무력.
- **보정 스칼라** $\tau(W)=\tfrac1N\sum_{x\ne0}\tfrac{A_W}M(1-\tfrac{A_W}M)$: $\tau=0\iff M_2=0$
  (124샘플 위반 0), Pearson $r(\tau,M_2)=0.76$, 추세 $M_2\approx7.7\,\tau^{0.32}$.

### I.5 선행연구 지형 & 신규성 방어선 (`magic-prior-art.md`)
- **주 선점:** **2605.05347**(EPFL, 프리프린트) — Shor magic↔주기 정량법칙·SRE도구·QFT손실까지.
  자기상관 객체($\Lambda$)도 여기 출처. **본 레포 잔여는 *다항우위 대조(Grover)*와 *부호이론 특수화*.**
- **부호↔magic 3갈래 점유:** (A) RM↔매직상태 증류/가중치열거자(다른 의미), (B) **Boolean
  비선형성↔하이퍼그래프 상태 magic(2308.01886 *Quantum 8,1351*; 2602.23687) — 최근접 이웃
  [전문 확인]**: 둘 다 **위상 인코딩**·**RM(2)/nonquadraticity**, *support/자기상관 미취급*,
  (C) Dicke/순열불변=구조적 support 평탄상태(대칭 특수case).
- **방어 가능한 좁은 잔여:** flat 마커상태 SRE 정확형 + 아핀⟺$A_W\in\{0,M\}$ 영점판정 +
  Grover 다중표시 적용 + 최소거리 무력성 반례 — **(i) 위상 vs support, (ii) RM(2) vs RM(1)** 로 차별화.

---

## Part II — 앞으로의 새 과제 (TODO)

### 과제 A ✅ 완료 — 부호이론 통계적 예측식 (Sidon 값)
랜덤 비아핀 $W$의 평균 magic = **Sidon 값** $\,M_2=\log_2(M^3/(7M-6))\approx 2\log_2 M-\log_2 7$
(명제 5). 랜덤 $W$는 $M\ll2^{n/2}$에서 whp Sidon → $\mathbb E[M_2]\to$ 이 값(수치 검증). 유한-$N$
보정은 가법충돌로 *하락*, $M^2/N$로 통제. `marker_code_closed_form.py` §2.
*잔여(M4로):* $M^2\gtrsim N$ 영역의 보정항 $\Delta(M,N)$ 해석 전개(현재 경험적 단조).

### 과제 B ✅ 완료 — 일반-M flat 마커상태 닫힌형
$M_2(|{\rm flat}_W\rangle)=-\log_2\!\big(\tfrac1{M^4}\sum_x E(W\cap(W{\oplus}x))\big)$, $E$=가법에너지
(명제 4). $M{=}1$·아핀 → 0 환원 확인, $\Lambda$(2605.05347)와 동형. 수치 $4{\times}10^{-15}$.
*잔여(선택):* 2-진폭 Grover(비-flat)의 명시 닫힌형(현재 명제 2′ 분해 + 명제 4로 환원 가능).

### 과제 C+D ✅ 완료 — 오라클-가림 = T-비용 가림 (명제 6)
$M_2(|\psi_f\rangle)>0\iff f$ 비선형 $\iff U_f$ 게이트분해가 비클리포드($T$) 요구 — 셋 다 $f$
아핀에서만 0(따름정리 1 + ANF 차수≥2 ⟺ Toffoli/$T$). 비선형 ANF 단항식 수 $T_{\rm proxy}$와
$M_2$ 동반 증가; Simon($0/0$) vs Shor modexp($T_{\rm proxy}{=}4{\to}156,\ M_2{=}1.54{\to}5.23$).
**2507.16543 Clifford-가림과 차별화**(오라클/추상화가 가림). `experiments/oracle_tcount_magic.py`.

### 과제 E (마감용) — (A)(C) 선행 전문 대조
정식 제출 전 RM↔증류/가중치열거자(A)와 Dicke/순열불변(C) 핵심 논문 **전문 1회 대조**
(현재 식별·요약 수준). arXiv egress 허용 환경 필요(현 환경 WebFetch 전면 차단).

---

## Part III — 로드맵

```
M1 (완료) ── 측정 인프라 + 명제 1·따름정리 1·2 + 명제 2/3/2′ + 사다리 + 42 assert
   │          산출: magic.py, magic-results.md, grover/shor_comb/oracle_magic, magic_proofs_check
   ▼
M2 (완료) ── 부호이론 보정 + 선행조사(B 전문확인)
   │          산출: marker_code_magic.py, magic-prior-art.md §5c (위상vs support / RM2 vs RM1)
   ▼
M3 (완료) ── 과제 B: 일반-M 닫힌형 = 가법에너지 (명제 4)
   │          + 과제 A: 통계 예측식 = Sidon 값 log₂(M³/(7M−6))≈2log₂M−log₂7 (명제 5)
   │          산출: marker_code_closed_form.py (닫힌형 4e-15, Sidon 1e-16, 랜덤W→Sidon 수렴)
   ▼
M4 (완료) ── 과제 C+D: 명제 6 (상태 magic ⟺ 오라클 T-비용, 둘 다 f 아핀에서만 0)
   │          + Δ(M,N) 해석 전개: 명제 5′ (E[ξ] 정확 닫힌형, MC 상대오차≤1.4e-2)
   │          산출: oracle_tcount_magic.py, marker_code_expected.py
   │          잔여: 과제 E(A·C 선행 전문 대조; 네트워크 필요)
   ▼
M5 (진행) ── 논문화 초안 작성: `magic-paper-draft.md` (명제 1–6 + 사다리 + 부호이론 통합)
              주 인용 2605.05347(Shor 선점) / 최근접이웃 (B) 차별화 / DTQW 2편 DOI
              잔여: 과제 E 전문 대조 후 정식 인용 확정 → arXiv 제출본
```

**의존성·우선순위:** A는 B와 독립 착수 가능하나, B(닫힌형)가 A(평균식)의 해석 골격을 줄 수 있어
**A·B 병행**이 효율적. C·D는 A·B와 독립(언제든). E는 제출 직전 마감 작업.

**위험·완화:**
- *A 닫힌형이 안 닫힐 위험* → 점근식 + 수치 fit(이미 $\tau^{0.32}$ 추세 확보)로 후퇴 가능.
- *시뮬레이션 규모 한계*($n\lesssim14$–16, $O(n4^n)$) → 해석적 푸리에 경로로 큰 $n$ 곡선 보강.
- *신규성 잠식*(2605.05347 출판 시) → 다항우위(Grover)·support 인코딩 각도는 (B)(2605)와 직교, 방어선 유지.
- *환경 egress 차단* → A·B·C는 로컬 계산만으로 가능; E만 외부 네트워크 필요.

---

## Part IV — 신규성·방어선 요약

| 항목 | 선행 | 본 프로그램 위치 |
|---|---|---|
| Shor magic↔주기 정량 | 2605.05347 (선점) | 재현·대조만 (신규 아님) |
| Grover(다항우위) magic 궤적·3-bit 포화 | 전용 선행 없음 | **신규(진짜 잔여)** |
| 오라클-가림(비선형성) | 2507.16543은 Clifford-가림(별개) | **신규(좁음)** |
| flat 마커상태 SRE = 자기상관, 아핀⟺$A_W\in\{0,M\}$ | (A)(B)(C) 어디에도 없음 | **신규(좁음, 방어가능)** |
| Boolean 비선형성↔magic | 2308.01886·2602.23687 [전문확인] (위상·RM2) | 차별화 대상(위상vs support·RM2 vs RM1) |
| 최소 해밍거리→magic | (인수인계서 제안) | **정정(반례로 무력화)** |

---

## 부록 — 파일 인덱스

| 파일 | 역할 |
|---|---|
| `magic.py` | SRE 측정도구 (XOR-FWHT, brute-force 검증) |
| `magic-results.md` | 명제·증명 (보조정리1·따름정리1·2·명제2/3/2′) |
| `magic-and-quantum-speedup.md` | 연구노트 (T3, 다음 과제: 부호이론) |
| `magic-prior-art.md` | 선행조사 (§5c 부호↔magic 3갈래, B 전문확인) |
| `magic-paper-draft.md` | 논문 초안 (명제 1–6 + 사다리 + 부호이론, 영문) |
| `magic-쉬운설명.md` | 초등 수준 직관 설명 |
| `experiments/marker_code_closed_form.py` | 명제 4·5 (가법에너지 닫힌형 / Sidon 법칙) |
| `experiments/marker_code_expected.py` | 명제 5′ (E[ξ] 정확 닫힌형) |
| `experiments/oracle_tcount_magic.py` | 명제 6 (상태 magic ⟺ 오라클 T-비용) |
| `experiments/grover_magic.py` | Grover magic 궤적·정점 포화·일반 M |
| `experiments/shor_comb_magic.py` | comb magic, 정정된 T1, Shor 대조 |
| `experiments/oracle_magic.py` | 오라클-가림 (아핀 ⟺ M₂=0) |
| `experiments/marker_code_magic.py` | 부호이론 보정 (자기상관·τ·d_min 반례) |
| `experiments/magic_proofs_check.py` | 명제 회귀검증 (42 assert) |
| `release_notes_v0.5.0.md` | v0.5.0 릴리스 노트 |
