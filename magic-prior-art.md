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

**인접(양자걷기 magic) — Grover의 최근접 이웃, 반드시 인용 [전문 확인 2026-06-17]:**
- **arXiv:2506.17783 — Mittal & Huang, "Quantum Magic in Discrete-Time Quantum Walk"**
  (Phys. Rev. Research, DOI 10.1103/7rwg-lhpv). DTQW의 SRE(1D 격자, single/two walker,
  coin states); magic↔얽힘 상보성, 노이즈 강건성. **전문 검색: grover/complete graph/
  marked/oracle/amplitude-amp = 0회.** "saturation"은 *walk 시간축* 포화(우리와 다른 축).
- **arXiv:2504.19750 — "Nonstabilizerness generation in a multiparticle quantum walk"**
  (Moca 외; **Phys. Rev. B 113, 075142 (2026), DOI 10.1103/nzrp-49mr** — 출판됨).
  1D doublon transport; 초기 $M_2\propto(v_ft)^2$, 장시간 $\sim\log[\ln(v_ft)]$, 경계반사 후
  *시간* 포화. **전문 검색: grover/complete graph/marked/oracle = 0회.**
- → **신규성 게이트 통과:** 둘 다 *1D 격자 수송*이지 **Grover 탐색/완전그래프/2-진폭 표시상태가
  아니며, 우리의 "문제크기 $n$의 함수로 정점→3 bit 포화(밀도→0)"와 축·기하·상태구조가 모두
  다름.** Grover-magic 결과(§B 아래)의 신규성은 유지. 단 이 2편을 DTQW-magic 최근접 이웃으로
  인용하고 "완전그래프/탐색/문제크기 포화"로 차별화 서술.

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

## 5c. 부호 이론 융합(표시집합→magic) 정밀 조사 — "다층 점유, 좁은 잔여" [2026-06-18]

`marker_code_magic.py`/T3의 "다음 과제: 부호 이론" 각도(표시집합 $W$를 고전 부호로 보고
$|{\rm flat}_W\rangle$의 SRE를 $W$의 대수적 특성으로 예측)를 정조준 재검색한 결과,
**부호↔magic 공간은 세 갈래로 이미 점유**돼 있다. 인수인계서의 "RM(1,n)까지의 최소
해밍 거리" 지표는 *고전 암호 Boolean 함수론*과 *하이퍼그래프 magic*의 기성 개념과
충돌하므로, 차별화가 필수다.

**점유된 것 (세 갈래):**
- **(A) RM 부호 ↔ 매직상태 증류 / 가중치 열거자 [전문 확인 2026-06-18: 2308.05152].** "Reed–Muller
  ↔ magic"은 수십 년 트랙: RM 부호의 transversal 비클리포드 게이트로 매직상태 증류(**PRX 2, 041021**,
  Campbell–Anwar–Browne; Hastings–Haah; **2510.10852** punctured RM 서브로그 증류). 부호의
  *가중치 열거자*가 증류 성능·SRE를 통제(**1702.06990** signed quantum weight enumerators;
  **2501.10163** 불변량 이론; **2308.05152** Quantum Lego: **higher-genus 가중치열거자=SRE**(식
  III.41–47)를 텐서망으로 계산 — 전문 확인). → **여기서 "부호→magic"은 *부호화/증류되는* 상태의
  가중치열거자(SRE *계산 도구*)** 지, *평탄 마커상태의 SRE를 support 자기상관으로 닫고 Sidon/아핀
  판정으로 예측*하는 우리 각도가 아니다. **용어 충돌 주의**(2303.11317식) — 명시적 구분, 포섭 없음.
- **(B) Boolean 함수 비선형성 ↔ 하이퍼그래프 상태 magic — *최근접 이웃* [전문 확인 2026-06-18, 저자제공 PDF].**
  **2308.01886 "Magic of quantum hypergraph states"** (Chen–Yan–Zhou, **Quantum 8, 1351 (2024)**, v2),
  **2602.23687 "SRE of 3-uniform hypergraph states"** (Kagamihara–Tsuchiya, v2 2026-05-14):
  하이퍼그래프 상태 $|G\rangle=U(G)|{+}\rangle^{\otimes n}=\prod_{e\in E}CZ_e|{+}\rangle^{\otimes n}
  =\frac1{\sqrt N}\sum_x(-1)^{f_G(x)}|x\rangle$ (2308.01886 **Def 1·Eq(1)**; $f_G$의 $c$차 단항식=$c$-edge;
  2602.23687 **Eq(3)**은 CCZ만→$f$ 3차). 즉 **Boolean 함수가 *위상*에, 받침은 항상 균일 $\mathbb F_2^n$
  전체.** magic의 원천은 **2차 초과(degree $\ge3$)=nonquadraticity**: 그래프상태(2-edge=2차위상)는
  클리포드·안정자(magic 0)이고 $\ge3$-edge에서 magic 발생(2308.01886 본문, Ref[49]=Liu–Winter
  PRX Q 3,020333의 nonquadraticity 최소화와 연결; 2602.23687 **Thm 1**: SRE$=$GF(2) 대칭행렬
  $C(x){+}C(x)^T$의 rank(2차형식)로 $O(N^3 2^N)$). 2308.01886은 평균차수 상계(**Thm 2**)·랜덤상태
  최대magic 집중·순열대칭(3-complete) 상태가 $\alpha{\ge}2$에서 상수/지수적 소(小) magic을 보임.
  → 우리와 결정적으로 다른 두 축: **(i) 인코딩 — 그쪽은 $f$가 *위상*(균일 support), 우리는 $1_W$가
  *support*(균일 위상); (ii) 부호 — 그쪽 magic 원천은 *RM(2,n)*(2차)까지 거리=nonquadraticity, 우리는
  *RM(1,n)/아핀부분공간*(1차)에서의 이탈.** 그쪽 지표는 *위상함수* 비선형성, 우리 지표는 *지시함수
  $1_W$의 자기상관*. **두 논문 모두 flat/indicator(support) 상태·받침 자기상관은 전혀 다루지 않음**
  (전 본문이 위상상태 전용; 전문 대조로 확인). **반드시 (B)를 최근접 이웃으로 인용하고 "위상 vs
  support / RM(2) vs RM(1)"로 차별화.**
- **(C) 구조적 support의 평탄상태 = Dicke/순열불변 상태 magic [전문 확인 2026-06-18: 2402.08551].**
  Dicke 상태는 곧 $W=\{x:\mathrm{wt}(x)=k\}$ 위 평탄상태(=비아핀 구조적 support). **arXiv:2402.08551
  "Nonstabilizerness of Permutationally Invariant Systems"** (Passarelli–Fazio–Lucignano): *순열대칭*을
  이용해 Pauli 계산을 Dicke 기저에서 $O(N^3)$로 환원, LMG·Dicke 상태에 적용. → **대칭 support
  특수case만 점유**(전문 확인: 임의/랜덤 $W$·자기상관·Sidon·아핀판정·Grover 없음). 우리 차별점:
  **임의/랜덤 $W$ + 자기상관 영점판정 + Sidon/E[ξ] + Grover 다중표시**(대칭 가정 없음).
- **(D) flat/SMF 상태 SRE 닫힌형 — *명제 4의 직접 선점* [전문 확인 2026-06-18].**
  **2311.08463 "Magic in generalized Rokhsar–Kivelson wavefunctions"** (Tarabunga–Castelnovo,
  *Quantum* **8**, 1347 (2024)). 임의 $|\psi\rangle{=}\sum_\sigma c_\sigma|\sigma\rangle$의 $M_2$를 4-copy 공식
  **Eq (8)** $e^{-M_2}{=}\sum_{\sigma^{(1..4)}}c_{\sigma^{(1)}}c_{\sigma^{(2)}}c_{\sigma^{(3)}}c_{\sigma^{(1)}\sigma^{(2)}\sigma^{(3)}}c^*_{\sigma^{(1)}\sigma^{(2)}\sigma^{(4)}}c^*_{\sigma^{(1)}\sigma^{(3)}\sigma^{(4)}}c^*_{\sigma^{(2)}\sigma^{(3)}\sigma^{(4)}}c^*_{\sigma^{(4)}}$
  ($\sigma\sigma\sigma$=점곱=XOR)로 쓰고, SMF면 $M_2{=}-\log(Z_M/Z^4)$ (Eq 9, 4-copy 분배함수). 균일
  flat 상태($c_\sigma{=}1/\sqrt M$ on $W$)를 대입하면, 전단사
  $(\sigma^{(1)},\sigma^{(2)},\sigma^{(3)},\sigma^{(4)}){=}(a,b,c,a{\oplus}b{\oplus}c{\oplus}x)$ 로
  **그들 Eq (8) = 우리 명제 4** $e^{-M_2}M^4{=}\sum_x E(W\cap(W{\oplus}x))$ 가 정확히 일치(검산 완료).
  ⟹ **명제 4(가법에너지 닫힌형)는 신규 아님 — 2311.08463 Eq (8)의 균일-$W$ 특수화 + 재명명.**
  "$2\log_2 M$"도 그들 Eq (7) $M_2{\le}4D_{\min}$ 및 표준 $M_\alpha{\le}2\log R$의 포화로 기지.
  **그러나 그들은 이 공식을 *물리 SMF 모델*(1D/2D/3D Ising, J1-J2, 삼각격자 AFM, EA 스핀글래스)
  에만 적용**; 임의/Grover 마커·부호이론은 다루지 않음. ⟹ 명제 4는 **인용·credit**하고, 우리
  기여는 아래 *부호이론 특수화*로 재포지션.

**우리 고유 객체의 출처 정직고지:** 자기상관 $A_W$/Walsh-4차모멘트라는 핵심 양 자체는
**2605.05347의 기하항 $\Lambda$**(균일진폭·유사난수위상 comb)에서 왔고, **flat 상태 SRE 닫힌형은
2311.08463 Eq (8)이 선점**(명제 4 = 그 특수화). 따라서 신규 기여는 "닫힌형/자기상관 아이디어"가
아니라 **그 공식의 *부호이론 특수화*: 아핀⟺$A_W\in\{0,M\}$ 영점판정 + Sidon 정확법칙 + 랜덤-$W$
정확 기댓값(명제 5′) + 최소거리 무력성 반례 + Grover 다중표시 적용**이다.

**진짜로 비어있는 좁은 잔여 (방어 가능, RK 선점 반영 후):**
- **부호이론 특수화** — (i) 아핀 ⟺ 자기상관 2값($A_W\in\{0,M\}$) 영점판정, (ii) Sidon 정확법칙
  $\log_2\frac{M^3}{7M-6}$(상수·유한형), (iii) 랜덤-$W$ 정확 기댓값 $\mathbb E[\xi]$(명제 5′ 4중쌍 분류),
  (iv) **Grover 다중표시(명제 2′)에 적용** — 은 (A)(B)(C)(D) 어디에도 없다(RK는 물리 SMF 모델만).
- **인수인계서 지표 정정의 박제값:** 최소 해밍 거리(=1차 비선형성류)는 magic을 결정하지 못함
  (`marker_code_magic.py` §1: $\{0,1,2,3\}$ vs $\{0,1,2,4\}$ 동일 $d_{\min}{=}1$, magic 0 vs 1.54).
  cf. bent 함수(아핀까지 최대거리)는 *위상* 함수 개념이라 support 마커와는 다른 대상.
- **명제 6(오라클-가림 = T-비용)** 은 RK와 무관 — 상대적으로 안전(단 ANF↔Toffoli·따름정리 1은 기지).

**조사 범위 메모(정직):** 이 환경은 **WebFetch가 전면 403**(arXiv·Quantum 저널·NASA ADS·
Semantic Scholar 본문 미확보) — WebSearch 스니펫 + HF 논문검색만 가능. **단 (B) 2편은 사용자가
PDF를 제공하여 전문 정독 완료**(2026-06-18): 결정적 2점(① 위상 인코딩 — 2308.01886 Def 1·Eq(1),
2602.23687 Eq(3); ② magic 원천=nonquadraticity/degree$\ge3$=RM(2)까지 거리 — 2308.01886 Thm 2·
Ref[49], 2602.23687 Thm 1의 2차형식 rank)을 **본문 식·정리로 확정**, 또한 **둘 다 support/indicator·
받침 자기상관은 전혀 다루지 않음**을 확인. → "위상 vs support / RM(2) vs RM(1)" 차별화 **최종 확정.**
**(A)(C)(D)도 전문 확인 완료(2026-06-18, 사용자 제공 PDF):** (A) 2308.05152=가중치열거자→SRE 계산도구,
(C) 2402.08551=순열대칭 환원(대칭 support만), (D) 2311.08463 Eq(8)=flat-state 닫힌형 선점. ⟹ **부호↔magic
4갈래(A·B·C·D) 전문 대조 종료.** 결론: **명제 4·"2logM"은 선점(D·표준한계)**, 살아남는 신규는
**부호이론 특수화**(아핀⟺$A_W{\in}\{0,M\}$·Sidon 상수·랜덤 E[ξ]·$d_{\min}$ 반례·Grover 적용) **+ 명제 6
+ Grover 3-bit 사다리**. (잔여 전문 미대조: 2510.10852·1702.06990·2501.10163 — 모두 *증류* 쪽이라 영향 없음.)

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
| 부호: RM ↔ magic상태 증류/가중치열거자 | PRX 2,041021, 2510.10852, 1702.06990, 2501.10163, 2308.05152 | 없음(다른 의미·용어충돌) |
| 부호: Boolean 비선형성 ↔ magic | **2308.01886·2602.23687** [전문확인] (위상인코딩, RM(2)/nonquadraticity, support 미취급) | 없음 — 단 위상 vs support·RM(2) vs RM(1)로 차별 |
| 부호: 구조적 support 평탄상태 magic | 2402.08551(순열불변), Dicke/대칭 SRE | 대칭case 점유 |
| flat/SMF 상태 SRE 닫힌형 (=명제 4) | **2311.08463 Eq (8)** [전문확인] (RK, 임의 $c_\sigma$) | **없음(선점)** — 명제 4는 그 특수화 |
| "$2\log_2 M$" 성장 | $M_2{\le}2\log R$ 한계 (표준), 2311.08463 Eq (7) | 없음(한계 포화) |
| **부호이론 특수화: 아핀⟺$A_W{\in}\{0,M\}$ + Sidon 정확법칙 + 랜덤 $\mathbb E[\xi]$ + Grover 적용** | (A)(B)(C)(D) 어디에도 없음(RK는 물리모델만) | **좁은 잔여(방어가능)** |
| 인수인계서 "최소 해밍거리→magic" | 반례로 무력화(`marker_code_magic.py` §1) | **정정(지표 부정확)** |

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
- **부호 이론 융합(§5c) 관련:**
  - PRX 2, 041021 — Campbell–Anwar–Browne, Magic-state distillation in all prime dimensions
    using quantum Reed–Muller codes (RM ↔ 증류)
  - 2510.10852 — Sublogarithmic Distillation in all Prime Dimensions using Punctured RM Codes
  - 1702.06990 — Signed quantum weight enumerators characterize qubit magic state distillation
  - 2501.10163 — Invariant Theory, Magic State Distillation, and Bounds on Classical Codes
  - 2308.05152 — Quantum Lego Expansion Pack: Enumerators from Tensor Networks (가중치열거자→SRE)
  - **2308.01886 — Chen, Yan, Zhou, Magic of quantum hypergraph states, Quantum 8, 1351 (2024)**
    [**전문 확인**] Def 1·Eq(1) 위상 인코딩 $\prod CZ_e|{+}\rangle^n$; Thm 2 평균차수 상계;
    nonquadraticity 연결(Ref[49]=Liu–Winter PRX Q 3,020333); support상태 미취급 — 최근접 이웃
  - **2602.23687 — Kagamihara, Tsuchiya, SRE of 3-uniform hypergraph states (v2 2026-05-14)**
    [**전문 확인**] Eq(3) ∏CCZ 위상; Thm 1 SRE=GF(2)행렬 $C{+}C^T$ rank, $O(N^3 2^N)$;
    support상태 미취급 — 최근접 이웃
  - 2402.08551 — Passarelli, Fazio, Lucignano, Nonstabilizerness of Permutationally Invariant
    Systems [전문 확인] 순열대칭→Dicke 기저 $O(N^3)$ 환원, LMG/Dicke만; 임의-W 미취급
  - 2308.05152 — Cao, Gullans, Lackey, Wang, Quantum Lego: Enumerators from Tensor Networks
    [전문 확인] higher-genus 가중치열거자=SRE(식 III.41–47) 텐서망 계산도구; 마커집합 부호이론 아님
  - **2311.08463 — Tarabunga, Castelnovo, Magic in generalized Rokhsar–Kivelson wavefunctions,
    Quantum 8, 1347 (2024)** [전문 확인] Eq (8) 임의 상태 4-copy SRE 공식 = **명제 4의 선점**
    (flat 상태로 특수화 시 우리 가법에너지 닫힌형과 전단사로 일치); 물리 SMF 모델에만 적용
  - (참고) 2510.01380 — Non-stabilizerness in quantum-enhanced metrological protocols (대칭상태 SRE 환원)
  - 2512.19657 — Extremizing Measures of Magic on Pure States by Clifford-stabilizer States (검색상 인접; 전문 미확보)

### 게재 상태/DOI (저널측 검토 2026-06-17)
*이 환경은 Crossref/저널 API egress 차단 — WebSearch로 저널 페이지를 읽어 확인. 구조화된
전수질의는 아님.*

**출판됨(DOI 확인):**
- 2106.12587 — Stabilizer Rényi entropy → Phys. Rev. Lett. 128, 050402 (2022),
  10.1103/PhysRevLett.128.050402
- Lami–Collura (2303.05536) → Phys. Rev. Lett. 131, 180401 (2023),
  10.1103/PhysRevLett.131.180401
- 2504.19750 → Phys. Rev. B 113, 075142 (2026), 10.1103/nzrp-49mr  *(이전엔 프리프린트로 오기)*
- 2506.17783 → Phys. Rev. Research, 10.1103/7rwg-lhpv
- 2306.09292 → Commun. Math. Phys. (2025), 10.1007/s00220-025-05421-3
- 2601.07824 — Computing quantum magic of state vectors → Quantum (2026), q-2026-04-10-2059

**프리프린트(저널 DOI 미발견, peer-review 전):**
- 2605.05347 (주 선점, Shor-magic), 2505.17185 (QAOA-magic), 2409.06659 (Amortized SRE),
  2512.24685 (XOR-FWHT)
- 2507.16543 (Krüger–Mauerer) — arXiv; IEEE QCE 2025 conference 판("Quantum Dark Magic") 존재.

→ **함의:** 주 선점논문 2605.05347은 아직 미출판이라, 본 레포 결과의 신규성 창은 열려 있음.
양자걷기 인접 2편(2504.19750·2506.17783)은 *출판된* 최근접 이웃이므로 인용 시 DOI 표기.
