# Magic(비안정자성)과 양자 속도우위 — 쇼어에서 무엇이 풀렸고 무엇이 남았나

*연구 과제 노트. shor 레포의 상태벡터 시뮬레이션과 직접 연결.*
*2026-06 개정: 선행연구 조사(→ `magic-prior-art.md`)로 핵심 질문 대부분이 이미
[Paviglianiti–Seclì–Tirrito–Savona, arXiv:2605.05347 (EPFL, 2026-05)]에서 해결됨을
확인. 본 노트를 "열린 문제 제안"에서 "**해결된 것의 정리 + 남은 잔여(T3)**"로 재포지션.*

---

## 0. 한 줄

> 쇼어가 magic을 *필요로 한다*는 것도, **magic의 *양*이 주기 $r$의 수론적 난이도와
> 정량적으로 묶인다**는 것도 이제 **확정**됐다(2605.05347, 해석이론 + 수치).
> 이 레포 고유의 잔여는 단 하나 — **Simon/다항 속도우위와의 대조(T3)**.

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

## 3. 남은 잔여 — T3: 선형(Simon) vs 비선형(Shor), 다항 vs 지수 속도우위

직접적 선행을 못 찾았고, **2605.05347 스스로 Discussion에서 명시적 open으로 남긴** 유일한
실질 입구:

> *"it would be interesting to study whether the behavior of magic changes in the
> presence of a polynomial quantum advantage, rather than an exponential one as in
> the case of Shor."* — 2605.05347, Discussion

### T3 가설
- **Simon 알고리즘**(주로 Hadamard + 오라클, $\mathbb F_2$-선형 숨은부분군)은 구조가
  "Clifford스럽다". 만약 Simon의 핵심 상태가 (오라클 표를 안정자로 둘 때) **낮은/0에 가까운
  $M_2$**를 보인다면, "선형(안정자스러운) 문제 vs 비선형(magic 필요) 문제"의 magic 차이가
  직접 드러난다.
- 대조축: **(i) 속도우위 유형**(Simon: 쿼리모델 지수분리; Grover: 다항 √N; Shor: 지수),
  **(ii) 문제의 대수 구조**(선형 vs 곱셈적 주기). magic 궤적이 (i)을 따르는지 (ii)를
  따르는지가 핵심 질문.
- 검증/반증 모두 명확: Simon의 $M_2$ 궤적이 Shor처럼 $\log(\text{구조크기})$로 자라며
  포화하면 "magic은 속도우위의 보편 통화", 거의 평탄하면 "magic은 *비선형* 구조의 표지".

### 이 레포에서의 실현
- 레포에 `experiments/simon_sigma_curve.py`가 이미 있다 → Simon 측 상태의 $M_2$ 궤적을
  Shor와 **동일한 SRE 파이프라인**(2605.05347이 쓴 FWHT/MPS, 또는 작은 $t$는 정확합)으로
  계산해 한 그래프에 겹쳐 그린다.
- Shor 측은 2605.05347의 셋업(반고전 QFT, 각 측정 직전 $M_2$)을 그대로 재현하면 비교가
  공정하다.

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
- 작은 경우($L\lesssim10$)는 $4^L$ Pauli 합으로 정확 계산.
- 큰 경우는 **fast Walsh–Hadamard(2512.24685)** 또는 **MPS Perfect Pauli Sampling
  (Lami–Collura, PRL 131 180401)** — 둘 다 2605.05347이 실제로 사용. state-vector용
  정확·고속 SRE는 2601.07824(Quantum 2026)도 참고.

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

## 부록 — 선행연구 전체 목록
→ `magic-prior-art.md` (2605.05347 전문 분석 + T-count·SRE·구조화상태 magic 관련 문헌).
