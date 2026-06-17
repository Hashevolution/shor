# 선행 연구 조사 — magic↔쇼어 속도우위 (magic-and-quantum-speedup.md 의 prior art)

*조사일 2026-06-16, 본문 정독 갱신 2026-06-17(2605.05347 전문+Supplemental 확보).
각 항목은 노트의 어느 부분과 겹치는지 표시.*

---

## 0. 한 줄 결론

> 노트의 **핵심 주장(magic의 양 ↔ 쇼어의 number-theoretic 난이도)은 이미 한 달 전
> 논문에 선점**되어 있다 — **arXiv:2605.05347, "The true cost of factoring:
> Linking magic and number-theoretic complexity in Shor's algorithm"**
> (Paviglianiti, Seclì, Tirrito, Savona, EPFL, 2026-05-06).
> **전문 확인 결과**: 이 논문은 노트의 **T1·T4·§5(SRE 도구)·"QFT 손실→실패"까지
> 전부** 다룬다 — SRE $M_2$(노트 §5와 동일식)를 fast Walsh–Hadamard + MPS Pauli
> 샘플링으로 계산하고, 주기를 $r=2^k\tilde r_{\text{odd}}$로 분해해 **magic이 홀수부
> $\tilde r_{\text{odd}}$로 통제**되며 $M_2\sim\log r\to L\log2$로 포화함을 해석이론으로 증명.
> **T1의 "2의 거듭제곱 ⟹ 안정자" 직관은 부분적으로만 맞음**(§1 정정 참조).
> **남은 진짜 잔여는 T3(Simon/poly-vs-exp 대조)뿐** — 이 논문도 이를 명시적 open으로 남김.

---

## 1. 치명적 선점 — 노트의 중심 주장 (T1·T4·§5·QFT역할) [전문 확인]

**arXiv:2605.05347 — Paviglianiti, Seclì, Tirrito, Savona (EPFL), 2026-05-06.**
전문(본문 6쪽 + Supplemental 4쪽) 정독 결과, 노트와 겹치는 부분을 항목별로:

- **측정도구(노트 §5)와 동일:** SRE $M_2(|\psi\rangle)=-\log(\frac1{2^L}\sum_P\langle\psi|P|\psi\rangle^4)$
  (식 1) — 노트 §5의 그 식. 계산은 **fast Walsh–Hadamard transform**(그들의 ref [55] =
  바로 노트가 §5에서 가리킨 2512.24685)과 **Lami–Collura의 MPS Perfect Pauli Sampling**
  (ref [35])으로 $O(4^L)$ 비용을 우회. → **노트가 "쓰겠다"던 도구를 이미 그대로 사용.**

- **T1(주기구조↔magic) = 이 논문의 메인 결과:** 주기를 $r=2^k\tilde r_{\text{odd}}$로 분해.
  회로 각 단계의 상태가 균일진폭+유사난수위상 중첩 $|\psi\rangle=\frac1{\sqrt D}\sum_{m\in
  D}e^{i\theta_m}|m\rangle$ 임을 보이고, 닫힌형 $M_2=4\log D-\log(4\Lambda+6D^2-5D)$ 유도.
  핵심: **평탄역(plateau)의 magic은 홀수부 크기 $D=\tilde r_{\text{odd}}$가 통제**하고,
  $M_2\sim\log r$로 자라 $r\gtrsim2^{L/2}$에서 상한 $L\log2$(최대 magic)로 포화.
  (N=18923 등에서 수치+해석 일치, Fig 3.)

- **T4(속도우위 연결) = 이 논문의 §"Period occurrence and success rate":**
  조건부 성공률 $S=g\cdot p_{\text{succ}}\sim(r/N)^\alpha,\ \alpha\approx1$ — 큰 $r$(=큰 magic)
  realization이 전체 성공확률을 지배. $N$이 커질수록 작은 $r$(저-magic) 기여는 소멸.
  → "magic이 속도우위/문제크기/성공확률과 정량 연결"이라는 T4가 **이미 수행됨**(Fig 4).

- **노트의 "QFT 역할 / 손실"(§2-3, §6):** 회로를 짧게 잘라 magic 생성을 막으면
  ($t<\lceil\log_2 r\rceil$) 성공확률이 즉시 0 — **"magic 손실 = 알고리즘 실패"**를 증명,
  근사 QFT/게이트 절단(48–50) 위험성을 경고(Fig 5). → 노트가 던진 "QFT가 magic을
  생성/재분배/구조화하나"에 대한 직접적 답변에 해당.

- **정정 — T1의 "2의 거듭제곱 ⟹ M₂=0"은 부분적으로만 참:**
  논문은 **$r=2$ (그리고 $\tau{=}1$의 $D{=}2$)만이 회로 전 구간 정확히 안정자**($M_2=0$)임을
  증명. $r=4,8,\dots$ 같은 **일반 2-거듭제곱은 회로상 magic을 생성**한다(Fig 3에서 r=4가
  0이 아님) — 평탄 magic은 $\tilde r_{\text{odd}}$가, 마지막 $k$스텝의 추가 성장은 $2^k$가 만든다.
  - 노트의 **이상적 comb 상태**(엄밀히 $r\mid Q=2^t$ 가정) $\frac1{\sqrt m}\sum_j|x_0+jr\rangle$
    한정으로는 "$r=2^s$ ⟹ (하위 $s$비트 고정)⊗($|+\rangle$들) = 안정자, $M_2=0$"가 **맞다.**
    그러나 **실제 쇼어에선 $r\nmid 2^t$가 일반**(그래서 연분수가 필요)이라, 흥미로운 경우의
    측정후 상태는 비균일 근사-comb이고 magic은 논문의 $\tilde r_{\text{odd}}$ 분석을 따른다.
  - 즉 노트 T1의 깔끔한 정리는 **자명한 이상화($r\mid Q$, 게다가 $r=2^s$는 쉬운 경우)**에서만
    성립하고, **비자명한 실제 케이스의 magic은 이 논문이 이미 정확히 기술**한다.

→ **결론: 노트의 헤드라인(T1+T4)과 §5 도구·"QFT 손실→실패"는 신규성이 사실상 없음.**
  반드시 2605.05347을 주 레퍼런스로 인용. 차별화하려면 §6의 잔여(T3, FFT-지름길 관찰)로.

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

## 5b. T3(속도우위 유형별 magic) 정밀 조사 — "반쯤 열림"

T3가 무주공산인지 정조준 재검색한 결과, **인접 선행이 상당함**. 구분:

**점유된 것:**
- **arXiv:2505.17185 — Capecci, Santra, Bottarelli, Tirrito, Hauke, "Role of
  Nonstabilizerness in Quantum Optimization" (2025-05).** QAOA의 SRE 궤적 = **"magic
  barrier"**(깊이 따라 상승→하강), 최종 magic↔성공확률, 단열까지. 후속 **2605.01620**
  (Hypergraph QAOA). → "알고리즘 속 magic 궤적→성공" *방법론*은 신규 아님; **최적화 분야 점유.**
- **arXiv:2507.16543 — Krüger & Mauerer, "Geometric and Resource-Theoretic
  Characterisation of Non-Stabiliserness in Quantum Algorithms"** (IEEE판 "Quantum Dark
  Magic"). **여러 알고리즘에 걸친 magic 궤적 비교 프레임**(대상: VQE·QAOA·QFT 등 변분/구조
  회로)을 이미 구축. 게다가 **"Clifford 연산에 가려진 magic을 드러내는 permutation-agnostic
  거리"** 도입 → 노트의 **"블랙박스가 magic을 숨긴다"** 각도의 *일부를 선점*(단 메커니즘은
  Clifford-가림). magic이 ~75% 깊이에서 정점 후 감소. → **비교 프레임·"숨김" 개념 점유.**
- **양 끝점 기지:** Simon·BV·Deutsch–Jozsa는 (선형 오라클) **클리포드 → magic 0**
  (Gottesman–Knill; Combarro 2021). Shor는 magic-rich(2605.05347). → "Simon은 magic 없이
  지수 쿼리 우위"의 절반은 **거의 자명**, 단독으론 약함.
- (참고) **arXiv:2303.11317 "Opening the Black Box Inside Grover's Algorithm"**(PRX 2024)는
  *오라클 구조/dequantization*이지 magic이 아님 — 용어 충돌 주의용 인용.

**인접(양자걷기 magic) — Grover의 최근접 이웃, 반드시 인용:**
- **arXiv:2506.17783 — "Quantum Magic in Discrete-Time Quantum Walk"** (Phys. Rev.
  Research 게재, DOI 10.1103/7rwg-lhpv). DTQW의 SRE 계산(1D 격자, coin/walker). **Grover
  탐색은 아님.**
- **arXiv:2504.19750 — "Nonstabilizerness generation in a multiparticle quantum walk."**
  walk magic가 **시간에 로그적 성장.**
- → Grover = *완전그래프 위 양자걷기*이므로 Grover-magic 연구는 이들을 인용해야 하고,
  "walk magic ~ log(time)"라는 기성 결과가 바를 약간 올린다. **그러나 Grover 탐색 특수형
  (2D 불변부분공간 회전)을 직접 다룬 건 아직 없음.**

**비어있는 것 (방어 가능, 단 더 좁아짐):**
- **Grover의 magic/SRE 궤적 — 여전히 진짜 빈칸.** Grover 자원 분석은 거의 *coherence·
  entanglement*(trace speed, Sci. Rep. 2020)뿐; 비교논문 2편(2505.17185, 2507.16543)
  **모두 Grover 미포함**; 양자걷기 magic 2편도 **Grover 탐색은 아님**. "2차(다항) 속도우위에
  magic이 필요한가/얼마나/스케일은"은 미개척. 2605.05347도 명시적 open으로 둠. → **가장
  강한 단일 타깃.** (단 착수 시 Grover 상태가 2D 실진폭 상태라 magic이 작을 가능성 유의.)
- **쿼리/오라클 모델의 "블랙박스 가림"** — 2507.16543의 *Clifford-가림*과 **다른 메커니즘**:
  오라클/FFT 블랙박스가 *비선형성*의 magic을 숨긴다(Simon 선형⟹0 vs Shor 비선형⟹>0).
  이건 미정식화 — 단, 이제 2507.16543과 명시적으로 차별화해 서술해야 함.

**조사 범위 메모(정직):** arXiv + 저널 DOI 양쪽 검색함. 이 환경은 Crossref/저널 REST API가
egress 차단이라 *구조화된 DOI 전수질의는 불가*; 대신 WebSearch가 PRX/PRR/Quantum/npj/IEEE/
IOP/Springer 등록본을 서버측에서 읽어 출판본을 잡음(이번에 PRR·IEEE·npj판 실제 포착).
게재상태: 2605.05347·2505.17185 = 프리프린트(저널 DOI 미발견), 2507.16543 = IEEE QCE 게재,
양자걷기 2편 = PRR 게재. → 핵심 Shor-magic 선점논문은 아직 peer-review 전.

---

## 6. 정직한 신규성 평가 + 다음 행동

| 노트 항목 | 선행연구 | 신규성 |
|---|---|---|
| §1 필요성 (Gottesman–Knill 등) | 교과서 | 없음 |
| §2·T1 magic↔주기구조 | **2605.05347** [전문확인] | **없음(선점, 해석이론까지)** |
| **T4** magic↔성공률·문제크기 | **2605.05347** §success-rate (Fig 4) | **없음(선점)** |
| §2-3 "QFT가 magic을 생성/재분배" + §6 손실 | **2605.05347** Fig 5 (손실→실패) | **없음(선점)** |
| §5 측정도구 (SRE via FWHT/MPS) | **2605.05347** Methods, 2601.07824, 2512.24685, 2106.12587 | 없음 |
| T2 QFT의 T-count(게이트분해) | 2409.06659, 2103.09999, 2306.09292 | 없음 |
| "comb SRE=0 ⟺ r=2의 거듭제곱" | 2605.05347: $r{=}2$만 엄밀 안정자, 일반은 $\tilde r_{\text{odd}}$가 통제 | **노트 주장은 부분적 오류** |
| T3 방법론(algo magic 궤적→성공) | **2505.17185**(QAOA)+**2507.16543**(변분/QFT)+2605.01620 | 없음(점유) |
| T3 "magic 숨김" 개념 | **2507.16543** (Clifford-가림 permutation-agnostic 거리) | 대부분 점유 |
| T3 끝점 Simon/BV=클리포드 | Gottesman–Knill, Combarro 2021 | 자명 |
| **T3 Grover의 magic/SRE 궤적** | 전용 선행 미발견; 비교논문 2편 모두 Grover 미포함 | **진짜 잔여** |
| T3 쿼리모델 오라클-가림(비선형성) | 미정식화(2507.16543의 Clifford-가림과 별개) | **잔여(좁음)** |

**다음 행동(우선순위):**
1. **재포지셔닝 필수.** magic-and-quantum-speedup.md는 2605.05347을 **주 레퍼런스로 인용**하고,
   T1·T4·§5·"QFT 손실"을 "이 논문이 해결한 것"으로 명시. 노트를 "열린 문제 제안"이 아니라
   **(a) 이 레포 시뮬레이터에서의 독립 재현/검증** 또는 **(b) 잔여 질문**으로 다시 써야 함.
2. **T1 정정.** "2의 거듭제곱 ⟹ magic 0"은 *이상적 comb($r\mid Q$)* 에서만 참이고 그 경우는
   자명/쉬운 케이스. 실제 비자명 케이스의 magic은 홀수부 $\tilde r_{\text{odd}}$가 통제(=논문 결과).
   노트의 가설 문장을 이 구분에 맞춰 수정.
3. **T3는 "반쯤 열림" — 좁혀서 공략(§5b 참조).** 방법론(QAOA, 2505.17185)과 Clifford 끝점
   (Simon/BV)은 이미 점유/자명이므로 피하고, **(a) Grover의 magic/SRE 궤적**(진짜 빈칸)과
   **(b) 오라클/FFT가 magic을 숨긴다는 통일 관점**에 차별점을 집중.
4. 만약 굳이 comb-state 각도를 유지한다면, 차별점은 **shor.py의 FFT 지름길이 게이트 magic을
   '컴파일해 숨긴다'는 §3 관찰**(논문은 semiclassical-QFT 회로 상태를 직접 보므로 이 "숨김"
   현상은 안 다룸) 정도 — 좁고 주로 교육적 기여.

---

### 인용 모음 (검증된 식별자)
- 2605.05347 — Paviglianiti, Seclì, Tirrito, Savona (EPFL), True cost of factoring:
  magic ↔ number-theoretic complexity in Shor (2026-05-06) — **전문 확인, 주 선점 문헌**
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
- 2505.17185 — Capecci, Santra, Bottarelli, Tirrito, Hauke, Role of Nonstabilizerness
  in Quantum Optimization (QAOA "magic barrier") — T3 방법론 선행
- 2507.16543 — Krüger & Mauerer, Geometric & Resource-Theoretic Characterisation of
  Non-Stabiliserness in Quantum Algorithms ("Quantum Dark Magic") — **T3 비교프레임 +
  "Clifford-가림 magic 드러내기" 선행 (Grover/Shor 미포함)**
- 2605.01620 — Structured Parameterization & Non-Stabilizerness in Hypergraph QAOA
- 2303.11317 — Opening the Black Box Inside Grover's Algorithm (PRX 2024; 오라클 구조,
  magic 아님 — 용어 충돌 주의)
- 2506.17783 — Quantum Magic in Discrete-Time Quantum Walk (Phys. Rev. Research,
  DOI 10.1103/7rwg-lhpv) — 양자걷기 magic; Grover의 최근접 이웃(단 탐색 아님)
- 2504.19750 — Nonstabilizerness generation in a multiparticle quantum walk
  (walk magic ~ log time)
- Combarro et al. 2021 (Comput. Math. Methods) — BV/DJ를 stabilizer formalism으로 설명
  (Simon/BV가 클리포드=magic 0임의 근거)
