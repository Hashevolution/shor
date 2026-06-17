# Magic(비안정자성)과 양자 속도우위 — 쇼어에서 무엇이 풀렸고 무엇이 남았나

*연구 과제 노트. shor 레포의 상태벡터 시뮬레이션과 직접 연결.*
*2026-06 개정: 선행연구 조사(→ `magic-prior-art.md`)로 핵심 질문 대부분이 이미
[Paviglianiti–Seclì–Tirrito–Savona, arXiv:2605.05347 (EPFL, 2026-05)]에서 해결됨을
확인. 본 노트를 "열린 문제 제안"에서 "**해결된 것의 정리 + 남은 잔여(T3)**"로 재포지션.*

---

## 0. 한 줄

> 쇼어가 magic을 *필요로 한다*는 것도, **magic의 *양*이 주기 $r$의 수론적 난이도와
> 정량적으로 묶인다**는 것도 이제 **확정**됐다(2605.05347, 해석이론 + 수치).
> 이 레포 고유의 잔여는 **다항 속도우위와의 대조(T3)** — 그 첫 결과로 **Grover의 magic은
> 정점 3 bit로 포화(밀도→0)하고 정답에서 0으로 되돌려짐**을 확인했다(§3, `grover_magic.py`).
> 같은 도구로 Shor 측 comb의 magic이 **문제크기 $t$에 비례해 증가**함도 확인(`shor_comb_magic.py`).
> → magic의 *양/밀도*가 속도우위 유형을 가른다: 클리포드-자명(0) → 2차(유한,밀도→0) →
> 지수(문제크기 비례→최대).

---

## 1. 선행연구가 이미 해결한 것 — 반드시 인용

**주 레퍼런스: A. Paviglianiti, M. Seclì, E. Tirrito, V. Savona, "The true cost of
factoring: Linking magic and number-theoretic complexity in Shor's algorithm,"
arXiv:2605.05347 (2026-05-06).** 이 논문은 본 노트가 던졌던 질문 대부분을 *이미* 답했다:

| 한때 "열린 문제"로 본 것 | 2605.05347의 결과 |
|---|---|
| magic의 *양* ↔ 주기 $r$ (T1) | **메인 결과.** $r=2^k\tilde r_{\text{odd}}$ 분해, 닫힌형 $M_2=4\log D-\log(4\Lambda+6D^2-5D)$. magic은 **홀수부 $\tilde r_{\text{odd}}$가 통제**, $M_2\sim\log r\to L\log2$ 포화 (Fig 3) |
| magic ↔ 성공확률·문제크기 (T4) | 조건부 성공률 $S\sim(r/N)^\alpha,\ \alpha\approx1$; 큰 $r$(=큰 magic)이 전체 성공을 지배 (Fig 4) |
| QFT가 magic을 생성/재분배/구조화? (§2-3) | 회로를 잘라 magic 생성을 막으면 ($t<\lceil\log_2 r\rceil$) 성공확률 즉시 0 — **magic 손실 = 실패** (Fig 5). 근사 QFT/게이트 절단 경고 |
| comb/주기 상태 SRE 측정도구 (§5) | **SRE $M_2$를 fast Walsh–Hadamard(2512.24685) + MPS Perfect Pauli Sampling(Lami–Collura)으로 계산** — 본 노트 §5가 "쓰겠다"던 바로 그 도구 |

→ **함의:** 아래 §2의 옛 T1·T2·T4와 §4 측정도구는 *신규 기여가 아니다.* 재현·검증
가치는 있으나, 새 결과를 원한다면 §3(T3)로 가야 한다.

---

## 2. 정정 — T1의 "2의 거듭제곱 ⟹ magic 0" 가설은 부분적 오류

원래 가설은 *"$r$이 2의 거듭제곱이면 comb 상태는 안정자($M_2=0$), 홀수 인수가 있으면
$M_2>0$"* 였다. 2605.05347 전문 대조 결과:

- **이상적 comb 상태**($r\mid Q=2^t$ 가정) $\frac1{\sqrt m}\sum_{j}|x_0+jr\rangle$ 한정으로는
  맞다: $r=2^s$면 (하위 $s$비트 고정) $\otimes$ ($|+\rangle$들) = 안정자, $M_2=0$.
  — 그러나 이 케이스는 **자명하고 쉬운 경우**다($r\mid 2^t$이면 연분수도 필요 없다).
- **실제 회로에선 $r=2$(및 $\tau{=}1$의 $D{=}2$)만이** 전 구간 엄밀히 안정자다.
  $r=4,8,\dots$ 같은 일반 2-거듭제곱은 **magic을 생성**한다(Fig 3에서 $r{=}4$는 0이 아님).
- 실제 쇼어는 $r\nmid 2^t$가 일반이고, 그때 측정후 상태는 비균일 근사-comb이며 그
  magic은 **홀수부 $\tilde r_{\text{odd}}$**가 통제한다 — 즉 "magic을 켜는 것은 *2의 거듭제곱
  여부*가 아니라 *주기의 홀수부*"다.

> **정정된 명제:** comb의 magic을 통제하는 것은 주기의 **홀수 부분 $\tilde r_{\text{odd}}$**이며,
> $2^k$ 인자는 회로 말미 $k$스텝의 추가 성장만 만든다. "인수분해를 비자명하게 만드는
> 주기"와 magic이 함께 켜진다는 직관 자체는 옳았으나, 그 메커니즘은 *2-거듭제곱 vs 아님*이
> 아니라 *$\tilde r_{\text{odd}}$의 크기*다. (이 모두 2605.05347이 이미 증명.)

---

## 3. 남은 잔여 — T3: 속도우위 유형별 magic (단, "반쯤 열린" 문제)

T3는 **무주공산이 아니다.** 정밀 조사 결과 인접 선행이 상당하다 — 정직하게 구분하면:

### 이미 점유된 것 (T3가 피해야 할 곳)
- **"알고리즘 속 magic 궤적 → 성공확률" 방법론**: Capecci–Santra–Bottarelli–Tirrito–Hauke,
  *"Role of Nonstabilizerness in Quantum Optimization,"* arXiv:2505.17185 (2025-05). QAOA의
  **"magic barrier"**(상승 후 하강) + 최종 magic↔성공확률, 단열까지. 후속 2605.01620.
- **비교 프레임 + "magic 숨김" 개념**: Krüger & Mauerer, *"Geometric and Resource-Theoretic
  Characterisation of Non-Stabiliserness in Quantum Algorithms,"* arXiv:2507.16543
  ("Quantum Dark Magic"). **여러 알고리즘에 걸친 magic 궤적 비교**(VQE·QAOA·QFT)를 이미
  구축하고, **"Clifford 연산에 가려진 magic을 드러내는 permutation-agnostic 거리"**를 도입 —
  노트의 "블랙박스가 magic을 숨긴다" 각도의 *일부를 선점*(메커니즘은 Clifford-가림).
- **양 끝점은 자명**: Simon·BV·Deutsch–Jozsa는 (선형 오라클) **클리포드 → magic 0**
  (Gottesman–Knill; Combarro 2021). Shor는 magic-rich(2605.05347). → "Simon은 magic 없이
  지수 쿼리 우위"의 절반은 **거의 자명** — Simon 단독은 약하다.

### 진짜로 비어있는 것 (더 좁아진 표적)
1. **Grover의 magic/SRE 궤적 — 여전히 빈칸.** Grover 자원 분석은 대부분 *coherence·
   entanglement*(trace speed, Sci. Rep. 2020)이고, 위 두 비교논문(2505.17185, 2507.16543)
   **모두 Grover를 안 다룬다.** "2차(다항) 속도우위에 magic이 *얼마나* 필요한가, 궤적·스케일은"은
   미개척. 2605.05347도 Grover/다항우위를 명시적 open으로 남김. → **가장 강한 단일 타깃.**
2. **쿼리/오라클 모델의 블랙박스 가림 — 2507.16543과 *다른 메커니즘*.** 그쪽은 Clifford 켤레가
   magic을 가리는 것이고, 여기서 말하는 건 **오라클/FFT 블랙박스가 *비선형성*의 magic을
   숨긴다**는 것: Simon(선형⟹분해해도 0) vs Shor(비선형 모듈러 곱⟹Toffoli⟹분해하면 >0).
   이는 §6의 FFT 지름길과 같은 현상 — 단 **반드시 2507.16543과 차별화해 서술.**

### T3 가설 (좁혀서)
- (a) **Grover**: $M_2$ 궤적을 반복수의 함수로. 2차 속도우위가 magic을 얼마나 요구하나, 그 양이
  탐색공간 $N$·해 개수와 어떻게 스케일하나 (Shor의 $M_2\sim\log r$ 포화와 대조).
- (b) **오라클-가림 가설**: "보이는 magic=0이어도 오라클을 게이트분해하면 magic이 드러난다"를
  Simon(선형⟹0)과 Shor(비선형⟹>0)로 정량화 → "magic은 회로 표면이 아니라 *문제의 비선형성*에
  산다." (2507.16543의 Clifford-가림과 구분: 이쪽은 *오라클 추상화*가 가리는 것.)

### T3(a) 결과 — 구현·확인됨 (`experiments/grover_magic.py`, `magic.py`)
단일 표시($M{=}1$) Grover 상태의 $M_2$를 직접 계산(SRE via XOR-FWHT, brute-force 검증):
- **궤적:** $M_2$는 0(균일 중첩=안정자)에서 올라 *탐색 중간*에서 정점, **정답 $k^*$에서 다시
  ≈0**(정답=계산기저=안정자). 즉 Grover는 magic을 **썼다가 되돌린다.**
- **닫힌형:** $\Sigma_P\langle P\rangle^4 = 1+(N{-}1)(a^2{-}b^2)^4+(N{-}1)(b^2(N{-}2){+}2ab)^4
  +(N{-}1)(\tfrac N2{-}1)(2b(a{-}b))^4$ (진폭 $a$ 1개, $b$ $N{-}1$개). 수치와 $10^{-13}$ 일치.
- **점근(핵심):** $N\to\infty$에서 $M_2\to-\log_2(a^8+(1-a^2)^4)$, $a^2=\tfrac12$에서 **정점
  $=-\log_2(1/8)=3$ bit**로 *포화*. 따라서 **magic 밀도 $M_2/n\to0$.**
- **대조:** Shor는 $M_2\to L$(밀도 $\to1$; 2605.05347). → **비안정자성의 양/밀도가 속도우위
  유형(다항 vs 지수)을 가른다**는 직접 증거. (단 양자걷기 magic 2506.17783/2504.19750을
  최근접 이웃으로 인용하고, Grover=완전그래프 walk라는 점에서 차별화해 서술.)
- **남은 일:** Grover 회로의 게이트분해 T-count 대조, 일반 $M$·다중해, 명제의 해석적 증명.
  (b) 오라클-가림 정식화는 아래에서 완료.

### T3 자기완결 대조 — `experiments/shor_comb_magic.py`
인용에만 의존하던 Shor 쪽을 같은 `magic.sre2`로 직접 계산:
- **정정된 T1 확인:** 측정후 comb $\frac1{\sqrt m}\sum_{x\equiv x_0\,(r)}|x\rangle$의 $M_2$는
  $r$이 **2의 거듭제곱이면 0**, **홀수 인수가 있으면 >0**. (단 $r(m{-}1)=Q{-}1$ 류의 우연한
  정렬은 $|j\rangle|j\rangle$ 꼴 → $M_2{=}0$; 예 $t{=}10,r{=}33$. "0 ⟺ 2의 거듭제곱"은 거의-참.)
- **스케일링:** comb magic은 **레지스터 크기 $t$와 함께 무한정 증가**(밀도 $M_2/t\approx0.4$–$0.55$).
  → Shor의 magic은 문제크기에 비례. (in-circuit의 더 강한 $M_2\to L$, 밀도$\to1$은 2605.05347.)

| 알고리즘 | 속도우위 | magic 거동 | 밀도 $M_2/n$ |
|---|---|---|---|
| Simon / BV | 지수(쿼리) | 0 (클리포드) | 0 |
| **Grover 탐색** | 2차(다항) | 정점→3 bit, 답서 0 | **→ 0** |
| Shor comb(측정후) | 지수 | $\propto t$ (무한정 증가) | ~0.4–0.55 |
| Shor in-circuit | 지수 | → L (최대; 2605.05347) | → 1 |

### T3(b) 결과 — 오라클-가림 정식화·확인 (`experiments/oracle_magic.py`)
함수 그래프상태 $|\psi_f\rangle=\frac1{\sqrt N}\sum_x|x\rangle|f(x)\rangle$ ($=H^n$ 후 오라클 친
상태; Simon/Shor 공통 중간상태)의 magic을 `magic.sre2`로 계산:
- **명제(수치 확인):** $M_2(|\psi_f\rangle)=0 \iff f$ 가 $\mathbb F_2$-아핀. 비선형이면 $M_2>0$,
  **비선형(AND)항 수에 따라 증가**(0→1.54→2.48→3.70).
- **Simon:** 숨은문자열 $s$의 2-대-1 함수를 *선형*으로 잡을 수 있어($f(x)=x\oplus(x_p\cdot s)$)
  오라클이 클리포드 → **회로 전체 $M_2\equiv0$인데도 지수 쿼리 속도우위.** ⟹ magic은
  (쿼리모델) 지수 속도우위에 *필요하지 않다.*
- **Shor:** $f(x)=a^x\bmod N$ 은 곱셈적=비선형 → 오라클이 magic을 강제, $M_2>0$이 문제크기와
  함께 증가($N{=}15{\to}1.54,\ 21{\to}5.80,\ 33{\to}7.12$).
- ⟹ **magic의 원천은 회로 표면이 아니라 *문제의 비선형성*이고, 쿼리 오라클(과 shor.py의
  `np.fft` 지름길)이 그 magic을 블랙박스로 *숨긴다*.** (2507.16543의 Clifford-가림과 구분.)

### 이 레포에서의 실현
- `experiments/simon_sigma_curve.py`로 Simon 측 $M_2$가 (오라클을 클리포드로 두면) 0임을
  대조군으로 박제하고, **Grover를 추가 구현**해 $M_2$ 궤적을 Shor(2605.05347 셋업 재현)와 한
  그래프에 겹친다. SRE 파이프라인은 §5 공통.

---

## 4. 보조 잔여 (작고 주로 교육적)

- **FFT 지름길의 "magic 숨김":** `shor.py`의 `simulate_period_finding`은 (작업레지스터
  선측정 → comb) + (역QFT를 `np.fft.fft`로 통째 적용)이라, 게이트 차원의 비-Clifford
  ($R_k$, Toffoli)가 FFT 블랙박스로 **컴파일되어 보이지 않는다.** 2605.05347은 *반고전
  QFT 회로의 단계별 상태*를 직접 보므로 이 "숨김" 현상 자체는 다루지 않는다 → 이 레포만의
  관찰로 짧게 문서화할 가치는 있으나 새 물리는 아니다.
- **이상적 comb($r\mid Q$)의 정확한 안정자 판정**을 §2의 정정과 함께 예제로 박제하면
  교육적. (자명 케이스임을 분명히 표시할 것.)

---

## 5. 측정 도구 (그대로 유효, 단 신규 아님)

순수상태 $|\psi\rangle$($L$큐비트)의 stabilizer 2-Rényi 엔트로피:
$$M_2 = -\log_2\!\Big(\frac1{2^L}\sum_{P\in\mathcal P_L}\langle\psi|P|\psi\rangle^4\Big),\qquad M_2(\text{안정자})=0.$$
- **이 레포 구현:** `magic.py` — `sre2(psi)`가 위 식을 **XOR-FWHT로 $O(n\,4^n)$ 정확 계산**
  ($\Sigma_P\langle P\rangle^4=\sum_{x,z}|\mathrm{WHT}[h_x](z)|^4,\ h_x(c)=\psi^*_{c\oplus x}\psi_c$).
  `sre2_bruteforce`(Pauli 행렬 직접, $n\le5$)로 검증 완료(랜덤·GHZ·$T$텐서 모두 일치).
- 더 큰 경우는 **MPS Perfect Pauli Sampling**(Lami–Collura, PRL 131 180401) — 2605.05347
  사용. state-vector용 고속 SRE는 2512.24685·2601.07824 참고.

---

## 6. 정직한 경계

- magic은 속도우위에 **필요**하나 단독으로 **충분**하진 않다(얽힘과 함께 작동).
- 이 레포의 지름길(작업레지스터 선측정 + FFT)은 **전체 유니터리가 아니라 한계 분포**를
  재현한다 — "회로가 *소비*하는 magic"(2605.05347이 다룬 것)과 "comb 상태가 *품은*
  magic"은 다른 양이므로 구분해 해석할 것.
- **magic↔속도우위의 정량 법칙은 *쇼어에 한해* 이제 알려졌다**(2605.05347). 일반 알고리즘으로의
  일반화(특히 다항 속도우위)는 미정 — 그래서 T3가 유효한 연구다.

---

## 7. 배경 — 왜 이 질문이 생겼나

분리된 한 작업(창발 시공간/홀로그래픽 코드, 별도 레포로 이관)에서, **비국소 magic이
홀로그래픽 코드의 넓이연산자를 상태의존으로 만든다**(중력 역반응; Cao et al.
arXiv:2306.14996)는 효과를 수치로 재현했다. 그 과정에서 **magic이 (i) 쇼어의 양자
속도우위와 (ii) 시공간 기하 둘 다의 공통 자원**임이 드러났다. 그 통찰이 "쇼어에서 magic의
*양*은 속도우위와 어떻게 연결되나"라는 질문으로 이어졌고 — 그 질문의 쇼어-특수적 답은
2605.05347이 주었다. 본 노트의 잔여 기여는 그 답을 **다른 속도우위 구조(Simon/다항)와
대조**(T3)하는 데 있다.

---

## 부록 — 명제·증명, 선행연구
- **명제·증명:** `magic-results.md` (보조정리 1: 평탄상태 안정자성 ⟺ 받침 아핀; 따름정리 1:
  그래프상태 $M_2{=}0\iff f$ 아핀; 명제 2–3: Grover 닫힌형·정점→3 bit). 재현 검증
  `experiments/magic_proofs_check.py` (35 assert).
- **선행연구 전체:** `magic-prior-art.md` (2605.05347 전문 분석 + T-count·SRE·양자걷기·
  알고리즘 magic 문헌, 신규성 게이트).
- **실험:** `magic.py`(SRE 도구), `experiments/{grover_magic,shor_comb_magic,oracle_magic}.py`.
