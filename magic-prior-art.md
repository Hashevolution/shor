# 선행 연구 조사 — magic↔쇼어 속도우위 (magic-and-quantum-speedup.md 의 prior art)

*조사일 2026-06-16. WebSearch 기반(이 환경은 arxiv.org egress 차단 — 본문 전체는 미열람,
제목·초록 요약에 의존). 각 항목은 노트의 어느 부분과 겹치는지 표시.*

---

## 0. 한 줄 결론

> 노트의 **핵심 주장(magic의 양 ↔ 쇼어의 number-theoretic 난이도)은 이미 1개월 전
> 논문에 선점**되어 있다 — **arXiv:2605.05347, "The true cost of factoring:
> Linking magic and number-theoretic complexity in Shor's algorithm" (2026-05).**
> T2(QFT의 magic/T-count)와 측정도구(Sec 5)도 각각 기성 문헌이 있다.
> 다만 **"comb 상태 SRE = 0 ⟺ r이 2의 거듭제곱"이라는 정확한 정리 형태**가
> 2605.05347 안에 명시돼 있는지는 본문 미열람으로 미확인 — 여기만 잔여 여지일 수 있음.

---

## 1. 치명적 선점 — 노트의 중심 주장 (T1·T4·Sec 2 프레이밍)

**arXiv:2605.05347 — "The true cost of factoring: Linking magic and
number-theoretic complexity in Shor's algorithm" (2026-05-06).**
- 초록 요지: 쇼어 인수분해에서 **비안정자성(magic) 생성**을 분석, **내재적 양자 복잡도와
  바탕 수론 문제의 계산적 난이도 사이의 깊은 연결**을 드러냄. **명시적 해석(analytic)
  이론**을 세워 magic이 알고리즘 성공에 하는 근본적 역할을 보이고, **쇼어 루틴이 그
  자원을 최대로(maximally) 활용**함을 보임. 게이트수·큐비트수 같은 표준 비용 대신
  **자원 기반 지표**로 보완.
- **노트와의 겹침:** 이게 정확히 노트 §0·§2·T1·T4의 명제다 — "magic이 단순 부산물이
  아니라 *문제의 난이도(수론 구조)와 묶여 있다*"가 이 논문의 제목 그 자체.
  → **노트의 헤드라인 기여는 신규성이 거의 없음.** 반드시 주 레퍼런스로 인용하고
  재포지셔닝해야 함.
- **확인 못한 것(잔여 여지):** 이 논문이 (i) **측정후 comb/주기 상태**를 직접 다루는지,
  (ii) **stabilizer Rényi 엔트로피**를 쓰는지, (iii) **"r이 2의 거듭제곱이면 magic=0,
  홀수 인수 있으면 >0"**이라는 깔끔한 정리를 명시하는지 — 초록만으로는 불명.
  → 본문 확보 필요(아래 §6).

---

## 2. T2 (QFT의 T-count / "QFT가 얼마만큼 magic을 공급하나") — 기성 문헌 있음

- **arXiv:2409.06659 — Amortized Stabilizer Rényi Entropy of Quantum Dynamics**
  (Zhu, Chen 외, HKUST, 2024-09). **SRE를 써서 QFT의 T-count 하한을 개선.**
  → T2가 묻는 바로 그 양("QFT가 공급하는 magic")을 monotone으로 정량화. **직격.**
- **arXiv:2103.09999 — Lower bound for the T-count via unitary stabilizer nullity.**
  유니터리 자체의 비안정자성으로 T-count 하한.
- **arXiv:2306.09292 — Stabilizer Testing and Magic Entropy via Quantum Fourier
  Analysis** (Comm. Math. Phys. 2025). QFT 해석과 magic 엔트로피를 직접 연결.
- 근사 QFT의 T-count 최적화: **arXiv:2203.07739**, Nature Sci. Rep. 2025
  (s41598-025-21087-2). "QFT magic의 게이트 차원 회계"는 잘 닦인 트랙.

---

## 3. 측정도구 (노트 §5) — 신규성 없음, 단 *실현 가능성엔 호재*

- **arXiv:2106.12587 — Stabilizer Rényi entropy** (Leone–Oliviero–Hamma,
  PRL 128 050402, 2022). 노트 §5의 $M_2$ 정의 원전.
- **arXiv:2601.07824 — "Computing quantum magic of state vectors"**
  (Quantum, 2026-04-10). **상태벡터로 주어진 순수상태의 SRE를 fast Hadamard transform으로
  정확·고속 계산.** → 노트가 "amps 벡터에 즉시 적용" 하려던 그 도구가 *이미 논문화*됨.
  $t\lesssim10$ 한계를 훨씬 넘김. **T1 수치실험은 이걸 그대로 쓰면 됨.**
- **arXiv:2512.24685 — Fast exact SRE via XOR–FWHT.** 같은 FFT류 가속.
- **arXiv:2601.00761 — Exponentially Accelerated Sampling of Pauli Strings for
  Nonstabilizerness.** 큰 $Q$용 샘플링(노트가 언급한 2501.12146류의 후속).
- **arXiv:2305.19152 — Efficient quantum algorithms for stabilizer entropies.**

→ 함의: 도구는 전부 있다(좋은 소식). 그래서 **방법론은 기여가 아니다.** 기여는 *결과*여야.

---

## 4. 구조화된 상태의 magic 닫힌형 — comb 상태의 방법론적 선례

- **arXiv:2508.03534 — Stabilizer Rényi Entropy for Translation-Invariant Matrix
  Product States.** GHZ·W·Dicke·spin-coherent 등 대표 상태의 **SRE 닫힌형**을 유도,
  **큐비트 수의 홀/짝(parity) 의존성**까지 명시(예: 일반화 GHZ
  $M=\log[8/(7+\cos 8\theta)]$, 홀짝 별도식).
  → "상태의 *구조*(주기/패리티)가 magic을 결정한다"는 **T1과 똑같은 패턴**이 이미
  확립됨. comb 상태 그 자체는 아니지만, **T1을 '신규 현상'이라 부르긴 어렵게** 만든다.

---

## 5. 배경(노트 §7) — 사실 확인됨, 활발한 분야

- **arXiv:2306.14996 (Cao et al)** 비국소 magic·홀로그래픽 코드 — 노트 인용대로 실재.
- 연관 최신: **arXiv:2601.03076** Multipartite Non-local Magic & SYK.
  → "magic이 속도우위와 시공간 기하의 공통 자원"이라는 배경 서사는 정당하나, 이 또한
  독립적으로 연구 중 — 노트의 차별점은 아님.

---

## 6. 정직한 신규성 평가 + 다음 행동

| 노트 항목 | 선행연구 | 신규성 |
|---|---|---|
| §1 필요성 (Gottesman–Knill 등) | 교과서 | 없음 |
| §2·T1·T4 magic↔수론 난이도 | **2605.05347** | **거의 없음(선점)** |
| T2 QFT magic/T-count | 2409.06659, 2103.09999, 2306.09292 | 없음 |
| §5 측정도구 (state-vector SRE) | 2601.07824, 2512.24685, 2106.12587 | 없음 |
| comb=구조→magic 패턴 | 2508.03534 | 방법론 선례 존재 |
| T3 Simon vs Shor 대조 | (직접 선행 미발견) | **잠재적 잔여** |
| "comb SRE=0 ⟺ r=2의 거듭제곱" 정확 정리 | 2605.05347 본문 미확인 | **확인 필요한 잔여** |

**다음 행동(우선순위):**
1. **2605.05347 본문 확보**(arxiv.org를 egress allowlist에 추가하거나 PDF 직접). 이 논문이
   (a) 측정후 comb 상태를 다루는지, (b) 2의 거듭제곱/홀수인수 분기를 명시하는지 확인.
   - 명시돼 있으면 → 노트는 **이 논문의 *수치 재현/검증*(이 레포 시뮬레이터)** 으로 재포지션.
   - 명시 안 돼 있으면 → "**측정후 comb 상태의 SRE에 대한 정확한 period-structure 정리 +
     FFT 지름길이 magic을 '컴파일해 숨긴다'는 §3 관찰**"이 좁지만 진짜 잔여 기여.
2. T1 수치실험은 **2601.07824의 fast-Hadamard SRE**를 그대로 채택($t\gg10$까지).
3. T3(Simon vs Shor magic 궤적 대조)는 직접 선행을 못 찾았으니 **가장 신규성 높은 입구**일
   수 있음 — 우선순위 상향 고려.

---

### 인용 모음 (검증된 식별자)
- 2605.05347 — True cost of factoring: magic ↔ number-theoretic complexity in Shor (2026-05)
- 2409.06659 — Amortized Stabilizer Rényi Entropy of Quantum Dynamics (QFT T-count 하한)
- 2103.09999 — Lower bound for T-count via unitary stabilizer nullity
- 2306.09292 — Stabilizer Testing and Magic Entropy via Quantum Fourier Analysis (CMP 2025)
- 2203.07739 — T-count optimization of approximate QFT
- 2106.12587 — Stabilizer Rényi entropy (PRL 128 050402)
- 2601.07824 — Computing quantum magic of state vectors (Quantum 2026-04-10)
- 2512.24685 — Fast exact SRE via XOR–FWHT
- 2601.00761 — Exponentially accelerated Pauli-string sampling for nonstabilizerness
- 2305.19152 — Efficient quantum algorithms for stabilizer entropies
- 2508.03534 — SRE for translation-invariant MPS (GHZ/W/Dicke 닫힌형, 패리티 의존)
- 2306.14996 — Cao et al, 비국소 magic & 홀로그래픽 코드 (노트 §7 배경)
- 2601.03076 — Multipartite Non-local Magic & SYK
