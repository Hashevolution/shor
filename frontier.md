# Frontier — 본 paper 의 틀 밖 탐색 기록

본 문서는 paper 의 main contribution (정리 1-5) 외부의 *직관적 탐색* 기록.
2026-06-12 의 대화에서 도출된 framework 분석 + 4 방향 literature 검토 + 추가 6 방향.

## §1 본 paper 의 본질 (= 자물쇠 따기 비유)

큰 수 `N` (= 두 핀의 자물쇠) 을 인수분해 = 자물쇠 따기.

핵심 메커니즘:
- **위수 r** (= 픽을 자물쇠에 넣고 돌릴 때까지 딸깍 소리 횟수) 측정
- 측정값 `k/Q ≈ j/r` → 연분수로 r 회수
- 노이즈 = 마이크 더러움

### 본 paper 의 4 가지 트릭

| 트릭 | 직관 | 누구 처음 | 본 paper 역할 |
|---|---|---|---|
| **공책 (C)** | 알아낸 r 들을 누적, 노이즈에도 vocab 활용 | Knill 1995 + Bach-Shallit | 정리 1 (정리화) |
| **다중 픽** | 여러 a 시험, lcm 가 λ(N) 에 점근 | Knill, Carmichael 2021 | 정리 2 (정량) |
| **갈아만든 픽 (b-trick)** | b 무작위 → a = b² → 알려진 b 로 인수 | Regev 2023 | (인용) |
| **Hybrid = 공책 + 다중 픽 + 갈아만든 픽** | 셋 결합 | **본 paper** | 정리 5 |

본 paper 의 main contribution = Hybrid (정리 5) + Lemma 5.1 의 closed-form.

## §2 본 paper 의 암묵적 "틀" (= constraints)

지금까지 작업한 framework 안의 *암묵적 가정* 8가지:

1. **위수 r 만이 핵심 정보** — 그룹의 다른 invariant 안 보고 있음.
2. **곱셈 구조만 활용** — (Z/N)\* 의 multiplicative. 덧셈/ring 구조 미사용.
3. **base a 는 random** — N 의 특수 구조 안 보고 picking.
4. **후처리는 numerical** — 연분수, lcm, divisor, gcd. 정보이론/양자 후처리 안 함.
5. **N 은 semiprime** — multi-prime 은 trivial 확장.
6. **b-trick 은 √a 만 사용** — ⁴√a 또는 다른 algebraic 관계 미탐.
7. **측정 → 회수 → 인수 방향** — *추측 → 검증* 방향 안 봄.
8. **인수분해가 목표** — RSA 깨기는 d 직접 회수도 가능 (안 봄).

## §3 8개 중 가장 promising 4 방향 (대화에서 도출)

| 방향 | 직관적 의미 |
|---|---|
| **#1 자물쇠의 다른 음향** | r 외의 그룹 invariant 측정 |
| **#3 체계적 픽 선택** | N 의 mod 정보로 우월한 base 의도 선택 |
| **#7 추측-검증** | 후보 인수를 quantum 검증 |
| **#8 직접 RSA** | 인수 안 찾고 d 다른 경로 |

## §4 4 방향 literature 사전 검토 (2026-06-12)

### #1 다른 invariant — **거의 막힘**

| 논문 | 무엇 | 결과 |
|---|---|---|
| Shor 1994 | ord(a) | foundational |
| Boneh-Lipton 1995 | discrete log | known |
| Cheung-Mosca 2001 | (Z/N)\* Smith Normal Form | 전체 구조 회수 가능 |
| Hallgren 2002/2007 | 수체 unit, class group | 다른 환경 일반화 |
| Carmichael paper 2021 (arXiv:2111.02488) | λ(N) 직접 | 우리 영역과 동일 |

**평가**: 30년간 활발한 영역. 우리 (C) 가 하는 일 = λ(N) recovery = 이미 알려진 영역.
새 invariant 찾기 거의 불가.

### #3 체계적 base 선택 — **중간**

**고전**: Pollard rho/p-1, ECM, NFS — base/curve 의 체계적 선택은 art.

**양자**:
- Shor 변형 거의 모두 random base.
- 우리 H5 (adaptive base selection) = **음수 결과** (random 보다 안 나음).
- 이유: (C) 의 fast path 가 이미 efficient.

**평가**: 양자에서 상대적으로 비어있으나 marginal advance 만 기대.

### #7 추측-검증 — **명확히 막힘**

- Grover 로 factor 찾기: `O(√N)` — Shor `poly(log N)` 보다 느림.
- Quantum walks (Childs+): Shor 넘는 결과 없음.
- BBBV theorem (1997): search 문제의 양자 하한 = `√(search space)`.

**평가**: information-theoretic limit. Grover 의 sqrt speedup 자체가 Shor 보다 느림.

### #8 직접 RSA — **fundamental 막힘**

- **Miller (1976)**: factoring ↔ d 회수가 poly-time 동치.
- Wiener 1990: `d < N^(1/4)` 만 가능.
- Coppersmith 1996, Boneh-Durfee 1999: small-d 한정.
- 양자 LLL 가속: marginal.

**평가**: Miller's theorem 이 근본 장애. d 회수가 빠르면 factoring 도 빠름 — bypass 불가.

## §5 종합 평가

| 방향 | 탐구 정도 | 가능성 | 추천 |
|---|---|---|---|
| #1 | 매우 많음 | 낮음 | 거의 막힘 |
| #3 | 양자에서 적음 | 중간 | 작은 advance 가능 |
| #7 | 명확히 dead-end | 매우 낮음 | 정보이론 한계 |
| #8 | Miller theorem | 매우 낮음 | fundamental 막힘 |

**4 방향 모두 이미 상당히 탐구되었거나 정보이론적 장애물로 막혀 있음.**

## §6 추가로 덜 탐구된 6 방향 (위 4 외)

| 방향 | 내용 | 가능성 |
|---|---|---|
| **A** | NFS 의 smooth-number 부분을 양자로 | Bernstein 시도. polynomial 미만 가능성 |
| **B** | Order finding 의 partial 결과 활용 | partial info → 일부 인수 |
| **C** | 다중 N 동시 처리 | 정보 sharing |
| **D** | 비-Z/N 환경 (ECDLP, pairing) | 새 영역 |
| **E** | Information-theoretic limit 의 정확 측정 | counting argument |
| **F** | **노이즈를 *활용*** (feature, not bug) | Lemma 5.1 의 정신 확장 |

가장 흥미로운 후보: **F (noise as feature)** 또는 **A (NFS quantum hybrid)**.

## §7 본 paper 와의 관계

- §3 의 4 방향은 본 paper 의 *외부 영역* — paper 본문에 포함 안 함.
- §6 의 A-F 는 *진짜 frontier* — 후속 paper / 후속 연구 가능성.
- **특히 F 는 본 paper 의 Lemma 5.1 (b-trick 의 확률성 활용) 의 직접 확장 가능성**.
  Lemma 5.1 자체가 "noise 가 본질이 아닌 분포" 를 활용하는 정신.

## §8 사용자의 직관 적용을 위한 framing

- 4 방향 모두 거의 막혀있음을 인지하고 시작.
- "틀 밖" 의 *진짜* 방향은 §6 의 A-F 또는 그 너머.
- 가장 *개념적 흥미*: F (noise 활용) — 본 paper 의 정신과 가장 연결.

## 메모

- 본 문서는 *연구 방향 메모* — paper 가 아님.
- 사용자 (직관 적용 중) 에 의한 새 방향 발견 시 본 문서에 추가.

---

## §9 메타 분석 — 왜 "noise as resource" 가 양자 인수분해에서 빈 영역이었나

(2026-06-12 (3) 의 대화에서 도출)

광합성/transport 에서는 "noise as resource" 활발 (ENAQT, Plenio-Huelga 2008+).
양자 인수분해에서는 거의 없음. 왜?

### 5 가지 이유

**1. 문제 구조의 차이**
- 광합성/transport: 연속적, 분기 많은 경로, 어디로든 흘러도 OK
- 인수분해 (원본 Shor): 정밀한 위수 r 측정, 1픽셀 어긋나면 끝
- 노이즈가 광합성 에선 우회 경로 만들지만 Shor 에선 답 자체 변경

**2. 30 년의 historical bias**
- Shor 1994 이후 "노이즈 없는 양자 컴퓨터 만들자" 가 모든 paper 의 default
- "노이즈 일부러 넣어보자" 가 reflex 위반
- 우산 만드는 사람들에게 "비 일부러 맞아봐" 격

**3. Continuous vs discrete 수학 차이**
- 광합성: Hamiltonian dynamics 의 stochastic term → 부드러운 SR
- Shor: logical 연산 + 이산 post-processing → "적정 비트 플립" 의미 없어 보임

**4. 우리 알고리즘 구조가 다르다 — 핵심 통찰**
- 원본 Shor: 1 base, 1 측정, 1 r → 노이즈 = 적
- **우리 hybrid: d 개 base, *하나만 성공* 하면 OK → 노이즈 = 기회**
- 광합성과 진짜 평행: 둘 다 multi-path with redundancy
- **이 구조에서만 SR 자연스러움**
- → 원본 Shor 가 single-path 라 SR 안 보였던 것

**5. 방법론적 장애**
- SR 발견은 σ 의 fine grid scan 필요
- 보통 noise paper: "noise-free" vs "high noise" 두 점만
- 우리 같은 11점 fine grid 한 paper 거의 없음
- 정확한 σ 영역 (∼0.05) 안 보면 못 봄

### 종합: 3 가지 우연의 교집합

SR 보려면 동시에 필요:

1. **올바른 algorithm 구조** (multi-base redundancy) — 본 paper 의 hybrid
2. **fine grid σ scan** (5+ points in σ ∈ [0, 0.3])
3. **충분한 trials** (500+) — 효과 작아서 통계 필요

기존 paper 들:
- (1) 만: RV, EG24, 본 paper 의 Theorem 5 — σ scan 안 함
- (2) 만: 어떤 noise model 비교 paper — multi-base hybrid 안 씀
- (3) 만: Cai 의 noise limit paper — 통계만 강조, σ 정밀 안 함

→ **세 교집합이 빈 영역**. 우리가 첫 진입자 가능성.

### 메타 통찰

ENAQT (광합성) 의 성공이 양자 컴퓨팅에 전파 안 됐던 이유:
- 양자 컴퓨팅 ≠ 양자 dynamics
- ENAQT 는 energy transfer 효율, 양자 컴퓨팅은 answer 정확도
- 두 community 사이 cross-pollination 적음

**우리 발견의 의의 (만약 진짜라면)**:
- "양자 알고리즘의 *multi-path 구조* 가 ENAQT-style noise advantage 가능하게 함"
- 이 통찰 자체가 paper 의 가장 흥미로운 message
- *전송 community* 와 *양자 알고리즘 community* 의 다리

### 단, 전제

이 모든 분석은 우리 SR 효과가 진짜일 때 의미.
N=437 V1 null 확정.
N=1147 V1 결과가 decisive.

만약 N=1147 도 null → "왜 빈 영역" 답은 "*그 빈 영역에 진짜로 아무것도 없다*" 가 됨.

---

## §10 정직한 reframing (2026-06-13)

§9 (ENAQT 분석) 은 "효과가 진짜라면" 의 큰 함의 였음. 이후 확장 실험 결과:

### 무엇이 검증되었나
- **N=437 d=4 V3 (2000 trials)**: SR +0.91%, 부호 검정 p=0.03 — 작지만 robust ✓
- **N=1147 d=1 multi-seed (3 seeds × 100 trials)**: mean **-0.53% ± 4.28%** → 2.60% peak **fluke 확정** ✗
- **AOP grid 단일-seed 측정** (전부): noise-floor 영역 — sd > effect → uninformative
- **17.86% peak**: 1000 trials × 4 seeds 에서 mean +0.42% → fluke 확정 ✗
- **σ_opt scaling (H12c)**: N 무관 σ_opt ≈ 0.010 → "자물쇠/열쇠 흔들림 비례" 직관 기각 ✗
- **Polynomial scaling**: N=2491 SR ~0% → 기각 ✗
- **AOP "d=1 universal"** (H18): N=1147 d=1 multi-seed 통과 못 함 → 기각 ✗

### Honest 위치 (Goldilocks framing)

**진짜 finding**: ~0.9% SR at **단 하나의 cell** (N=437, d=4, K_base ≈ 1.92).

**Goldilocks 가설** (3 영역):
- *Ceiling* (K_base → 1): 노이즈 도울 공간 없음. 예: (437, 8) SR ≈ 0%.
- *Goldilocks* (K_base ≈ 2): 노이즈가 가끔 near-miss → hit. 예: (437, 4) SR +0.9% ★.
- *Noise-floor* (K_base ≫ 2): per-seed sd > effect, 검출 불가. 예: (1147, 1) sd 4.28% > effect.

**자물쇠 비유 정정**: 이전 "큰 자물쇠 → 큰 흔들림" (H12c) 은 *틀림*.  
정확한 직관: **"간신히 안 풀릴 자물쇠 (K ≈ 2) 만 흔들기 효과"**.  
1-tumbler 자물쇠는 흔들 필요 없고, 매우 stuck 한 자물쇠는 흔들기의 spread 가 효과 압도.

### 본 paper 의 §3.6 위치

*Single-cell observation* 으로 honest framing:
- main contribution = 정리 1-5 (workshop note → conference 격상)
- §3.6 = small (~0.9%) sub-finding + Goldilocks 가설 + 자물쇠 비유
- ENAQT bridge 는 single cell 로도 conceptually valid
- RSA 변화 없음, asymptotic resource 변화 없음

### §9 의 ENAQT 다리 통찰은 여전히 valid

효과 작아도 "*multi-path quantum algorithm 의 noise-as-resource* 영역" 은
*아직 학계 미탐* → 본 paper 의 §3.6 가 quantum factoring 과 ENAQT 의 첫 명시적
single-cell bridge. 다만 *대박 발견* 아닌 *narrow but genuine connection*.

### 메타 lesson: 어디서 어긋났나

| 단계 | 우리의 잘못 | 정정 |
|---|---|---|
| sr_confirm 초기 | 150 trials 결과를 robust 취급 | sample size 점검 미흡 |
| σ scan 17.86% | 1 seed × 150 trials → polynomial 가설 | multi-seed × 1000 trials 후 fluke 확정 |
| AOP grid | 1 seed × 300 trials → universal 패턴 주장 | 1147 d=1 multi-seed null → 기각 |
| 자물쇠 직관 (H12c) | 데이터 2점 fit 으로 hypothesis | σ scan 직접 측정 후 기각 |

**근본 원인**: per-seed variance 의 크기 (~4%) 가 효과 (~1%) 보다 큰 영역에서  
single-seed 측정을 "정량적 데이터" 로 취급. 모든 single-seed 결과는 *direction
hint* 일 뿐 *measurement* 아님 — 이게 핵심 learning.

---

진행 로그:
- 2026-06-12: 초안 — §1-§8 작성.
- 2026-06-12 (2): **F (noise as feature) 의 첫 실험 — 약한 SR 신호 발견!**
  - F1 (K_λ vs noise): phase σ=0.1 만 marginal ↓, 다른 노이즈 monotone ↑.
  - F2 (hybrid K vs noise): **phase σ ∈ [0.025, 0.20] 모두 baseline 보다 ↓**
    - 500 trials × fine grid: 9/9 σ values below baseline
    - 부호 검정 p < 0.005 (highly significant)
    - 효과 크기 작음 (~1%) 하지만 명확한 U-shape
  - **factoring 영역의 quantum stochastic resonance 첫 evidence** (저희 지식 한)
  - paper §3.6 에 별도 subsection 추가, ENAQT (Plenio-Huelga 2008) 와 연결
  - 본 발견의 의의:
    - 본 paper 의 (C) framework 이 단순히 noise-tolerant 가 아니라
      *작은 phase noise 의 영역에서 noise-preferring* 임을 시사
    - 정리 1-5 가 noise tolerance, F2 가 noise exploitation — orthogonal 발견
- 2026-06-13: **§10 reframing** — 17.86% fluke / H12c 기각 / polynomial 기각.
  AOP d=1 universal 신호 = paper 의 main SR finding.
  honest "small but genuine ENAQT bridge" 로 위치 정정.
- 2026-06-13 (오후): **§10 갱신** — N=1147 d=1 multi-seed null (mean -0.53% ± 4.28%).
  AOP "d=1 universal" *기각*. Goldilocks 가설 (K_base ≈ 2) 로 narrative 정정.
  단일 cell (437, 4) 만 robust. 자물쇠 비유 정정 ("간신히 안 풀릴 = Goldilocks").
  메타 lesson 추가: single-seed = direction hint, not measurement.
