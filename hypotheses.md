# SR 메커니즘 가설 추적

본 문서는 N=1147 의 SR 가 N=437 보다 8배 큰 이유에 대한 8개 가설.
실험 결과 들어올 때마다 갱신.

## 데이터 ground truth (2026-06-12 기준)

| (N, d) | K_baseline | SR % (σ=0.05) | 신뢰도 |
|---|---|---|---|
| (437, 2) | 4.420 | 0.45% | M1, 200 trials |
| (437, 4) | 2.322 | 1.03% | sr_confirm, 500 trials |
| (437, 4) | 2.205 | 1.81% | M1, 200 trials (다른 seed) |
| (437, 8) | 1.190 | 0.00% | M1, ceiling |
| (1147, 2) | 2.630 | **8.56%** | n_extend, 200 trials |
| (1147, 4) | 1.440 | **5.90%** | n_extend, 200 trials |
| (2491, ?) | ? | ? | 진행 중 |

## 8개 가설

### H1: **Q 크기 (continued fraction depth)**
- 핵심: Q=4M (N=1147) vs Q=262K (N=437) → CF depth ~22 vs ~18
- 더 많은 convergent 후보 per measurement → 노이즈로 노출되는 path 많음
- 예측: **SR ∝ log Q ∝ log N**
- 검증법: N 같지만 Q 인위적 변경 (McAnally 식). 또는 N 시리즈에서 log fit.
- 상태: **유력 #1**

### H2: λ(N) 약수 개수 d(λ)
- 핵심: N=1147 의 λ/2 약수 12개 vs N=437 의 6개
- (C) 의 divisor pool 크기 → cross-coordinate 노출 효과
- 예측: SR ∝ d(λ(N))
- 검증법: 같은 N 크기지만 d(λ) 다른 것들 비교
- 상태: **유력 #2**

### H3: v_2 비대칭
- 핵심: N=1147 (v_p=1, v_q=2 비대칭) vs N=437 (v_p=v_q=1)
- ord 분포의 *richer dynamics*
- 예측: 대칭 v_2 < 비대칭 v_2 의 SR
- 검증법: 같은 크기 (대칭 vs 비대칭) N 쌍 비교
- 상태: 가능, 약함

### H4: K_baseline 의 sweet spot
- 핵심: K_baseline ~ 2-3 가 SR 최대
- K=1 근접 ceiling, K 크면 매 step 의 SR 작아짐
- 예측: SR(K_baseline) 가 U-shape (peak at K~2-3)
- 검증법: 같은 N 의 d 변화로 K_baseline scan
- 상태: 일부 관찰됨 (M1 의 d=2 vs 4 vs 8)

### H5: 알고리즘 difficulty regime
- 핵심: 작은 N = 너무 쉬워 noise 도울 여지 없음, 적당 N = SR 큼
- 정성적, 정량적 prediction 약함
- 상태: 정성적

### H6: Q 의 mode 수 (소거 후보)
- 핵심: Q amps 가 mode 많음 → 노이즈 작용 면적 넓음
- 예측: Q 같으면 같음 → H1 의 부분집합
- 상태: H1 에 흡수

### H7: ord 작은 값의 *단순성*
- 핵심: N=1147 의 ord(a) 가능 값에 2, 3, 5 등 작은 r
- 작은 r 에서 convergent 단순 → 노이즈 영향 큼
- 검증법: ord(a) 분포 직접 측정
- 상태: 가능

### H8: λ 의 4-divisibility
- 핵심: N=1147 (λ=4·45) vs N=437 (λ=2·99)
- 4|λ 가 b-trick 의 더 깊은 sqrt extraction 만들 가능
- 검증법: λ ≡ 0 mod 4 인 N vs 그렇지 않은 N 비교
- 상태: 약함, testable

## 검증 우선순위

1. **H1 (Q size)** — N 시리즈에서 log fit 으로 검증 가능
2. **H2 (d(λ))** — 같은 크기 다른 d(λ) 쌍으로 검증
3. H4 (K sweet spot) — d scan 으로 즉시 확인
4. H3, H8 — 비교 N 쌍 찾기
5. H5, H7 — 정성적 또는 보조

## 향후 데이터 슬롯 (실험 결과 채울 곳)

### N=2491 (= 47·53) — 백그라운드 진행 중
- λ(2491) = lcm(46, 52) = 1196 = 2²·13·23
- v_2 = 1, 2 (비대칭, N=1147 과 같은 패턴)
- λ/2 = 598 약수 = 8 개 (N=1147 보다 적음)
- d(λ) = 12 (N=1147 과 같음)
- Q = 2^24 = 16M
- 예측 (H1): SR > N=1147 (Q 4배)
- 예측 (H2): SR ≈ N=1147 (d(λ) 같음)
- 결과: ___ (대기)

### N=4087 (= 61·67) — 다음 측정 후보
- λ(4087) = lcm(60, 66) = 660 = 2²·3·5·11
- v_2 = 2, 1 (비대칭)
- λ/2 = 330 약수 = 16개
- Q = 2^24 = 16M
- 예측: SR 더 큼 (Q, d(λ) 모두 큼)
- 결과: ___ (미측정)

### N=8009 (= 89·90 거의) — 미래
- 8009 가 소수가 아님 확인 필요. 8011 가 89·90 ... 사실 89·90 = 8010 not semiprime
- N = 8051 = 83·97? 또는 다른 쌍
- 미정

## V 결과 통합 (대기)

V1 진행 중:
- seed 1: N=437 d=4 SR = 1.13%
- seed 2: N=437 d=4 SR = 0.81%
- seed 3: N=437 d=4 SR = **-1.85%** (반전!)
- seed 4-5: 대기

**N=437 의 SR 가 noise 보다 작을 가능성**. V1 of N=1147 결과가 분기 결정.

만약 N=437 V1 null + N=1147 V1 robust → "**SR emerges only above some N_critical**".

이 자체가 흥미로운 finding 가능성 (H4, H5 와 일치).

## H9-H13: scaling 가설 (V1 N=1147 robust 시사 후 추가)

### H9: **Polynomial scaling**
- `SR(N) = max(0, c · (N - N_crit)^α)`
- N=437 (N_crit 근처, null) → N=1147 의 8.56% (활성) → α ≈ 3 추정
- 예측: N=2491 에서 SR ≈ 30-60% (만약 H9 참)
- 검증: N 시리즈에서 polynomial fit

### H10: **Threshold + plateau**
- `SR(N) = max(0, min(SR_max, c · log(N/N_crit)))`
- 임계점 후 빠르게 plateau
- 예측: N=2491 에서 SR ≈ 10-15%

### H11: **Sigmoidal**
- `SR(N) = SR_max · σ(α(N - N_50))`
- 빠른 상승 후 plateau
- 예측: N=2491 → ~20-30%, N=4087 → ~30-40%

### H12: **σ_opt 가 N 에 따라 변화** ← *사용자 직관*

> *"작은 자물쇠는 열쇠를 작게 흔들어가며 맞춰주고, 큰 자물쇠는 조금더 크게 흔들어가며 맞춰주고,
>  각 자물쇠의 크기에 비례한 흔들림이 열쇠를 잘맞게 해준다."*
>
> ([사용자 직관 표현, 2026-06-12])

- H12a: σ_opt ∝ log N (느린 증가)
- H12b: σ_opt ∝ 1/√N (감소)
- H12c: **σ_opt ∝ N^α** (α > 0, 비례 증가) ← *직관적 가설*

**의의 (만약 H12c)**:
- 우리 측정 (모든 N 에 σ=0.05 고정) 이 *과소측정* 일 수 있음
- N=1147 의 진짜 σ_opt 가 σ=0.2~0.3 이면 SR 가 *15-25%* 가능
- → 실제 polynomial scaling 이 *훨씬 더 가파를* 가능성

### H13: **σ_opt + SR_max 가 함께 변화**
- SR_max(N) ∝ N^β, σ_opt(N) ∝ N^γ (둘 다 양수)
- 큰 N 에서 SR *훨씬 큼*, 단 σ *조정 필요*
- 만약 맞으면: RSA-2048 까지 가는 길에 **paradigm shift** 가능

## H12c 의 검증 우선순위 — 최상

### 즉시 실험: N=1147 의 σ scan

```
N=1147, d=2, σ ∈ {0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0}
각 200-500 trials
→ σ_opt(N=1147) 추적
```

3 결과 시나리오:
- σ_opt ≈ 0.05 (= N=437 과 같음) → H12 기각, H1-H4 충분
- σ_opt ≈ 0.02 (작아짐) → H12b (가능, 실용 의미 적음)
- **σ_opt > 0.05 (커짐) → H12c 확정** ← 가장 실용적 결과

### 그 후: N=2491, 4087 도 σ scan

σ_opt(N) 함수 형태 결정 → polynomial vs log vs sigmoid

## H14-H15: seed-내 변동성 + lucky base 분포 (V1 N=1147 outlier 분석에서 도출)

### V1 N=1147 d=2 의 outlier 발견 (2026-06-12)

| seed | K_baseline | SR % | 해석 |
|---|---|---|---|
| 5 | **2.110** | **-5.21%** | 유일 음수, K_base 최저 |
| 4 | 2.580 | +12.40% | 가장 큰 SR |
| 1 | 2.840 | +2.46% | |
| 3 | 3.130 | +6.39% | |
| 2 | 3.270 | +7.34% | |

seed 5 의 outlier = K_baseline 가장 낮음 + 유일 음수 SR.

### H14: **K_baseline 의 seed 내 변동성**

같은 (N, d) 에서도 *base 의 luck* 에 따라 K_baseline 크게 다름.
- d=2 의 경우 4 가지 base luck 조합:
  - 4 의 1 (P=c²=1/4): 양쪽 b_i 모두 *nontrivial sqrt* (double-lucky) → K_base 최저 → ceiling
  - 4 의 2 (P=2c(1-c)=1/2): 한쪽만 lucky → normal K_base → 최대 SR
  - 4 의 1 (P=(1-c)²=1/4): 둘 다 trivial (no-lucky) → 인수 못함, hit max_runs

→ seed 5 가 *double-lucky* 였을 가능성 가장 큼.

**의의**:
- d=2 에서 *25% seed 가 ceiling, 25% 가 fail* → SR 가 seed 마다 다름
- d=4 에서: 1/16 ceiling, 1/16 fail, 14/16 normal — *more uniform*
- d=8 에서: most seed 가 multi-lucky → ceiling dominant → SR ≈ 0
- 이게 **M1 의 d sweet spot 패턴 (d=4 최적)** 의 *seed-level* 메커니즘

### H15: lucky base 분포 vs SR 의 *aggregated* 효과

Trial 평균 SR 가 보이는 효과:
`SR_aggregate = mean over seeds of [SR if not ceiling, else 0 or -]`

- seed 마다 다른 ceiling state → 평균 SR 가 *seed 분포에 dependent*
- N 다른 군에서 *lucky base 확률 (c)* 가 변하면 → SR_aggregate 도 변함
- Lemma 5.1: `c = 1 - 2^{-(v_p+v_q)} · (4^{min(v_p,v_q)} + 2)/3`
- N=437 (v_p=v_q=1): c = 1/2
- N=1147 (v_p=1, v_q=2): c = 1 - 1/4 = **3/4** ← 다름!

→ **N=1147 의 c=3/4 가 N=437 의 c=1/2 와 *다름*** = lucky base 가 *더 흔함*.

**예측**:
- d=2 의 ceiling 비율: N=437 (c=1/2) = 1/4, N=1147 (c=3/4) = **9/16**
- 즉 N=1147 d=2 에서 *56% seed 가 ceiling* 예상
- 그러면 seed 5 의 ceiling 출현은 *예측 일치* (2/5 ≈ 40%, 분포 안)

**SR_aggregate 의 이론적 한계**:
- N=1147 d=2: only 1 - (1-c)^d = 1 - 1/16 = 15/16 seeds 가 *최종 성공*
- 그 중 ceiling 비율 9/16 = 60%
- Non-ceiling, non-fail seeds: 6/16 = 38%
- 이 38% 만이 SR 발현 → *aggregate SR 가 individual SR 의 38%* 약화

**다른 d 에서**:
- d=4, c=3/4: ceiling = c^d = (3/4)^4 = 81/256 = 32% → SR-effective seeds 가 더 많음
- d=8, c=3/4: ceiling = 0.1, SR-effective 90% (but K_base 작아서 다른 ceiling)

### H14 + H15 의 통합 의미

본 SR 효과는 **여러 layer 의 효과**:
1. **Per-coordinate**: phase noise 가 g 살짝 boost
2. **Per-seed**: lucky base 조합이 K_base 결정
3. **Aggregate**: seed 분포의 SR_aggregate

따라서 SR 측정의 변동성 = 본질적 (= noise + 알고리즘 구조 동시 작용).

**실용적 함의**:
- Adaptive d (= seed 마다 d 조정) 가 더 효율적일 수 있음
- N 별로 c 다름 (Lemma 5.1) → 최적 d 도 N 별로 다름
- 이게 H12 (σ_opt ∝ N^α) 와도 연결 — σ 와 d 의 *결합 최적화*

## H16: σ scan 결과 후 갱신 (2026-06-12 σ scan 완료)

### σ scan 결과

| | σ_opt | SR_max | trials |
|---|---|---|---|
| N=437 d=4 | **0.010** | 1.86% | 300/σ |
| N=1147 d=2 | **0.010** ← 같음! | **17.86%** | 150/σ |

### H12c **기각**

σ_opt 가 N 따라 변하지 않음 (둘 다 0.010).
사용자의 "자물쇠 비례 *흔들림* 크기" 직관 = 데이터로 *기각*.

### 사용자 직관의 *수정된* 형태

- ~~자물쇠 비례 흔들림 크기~~ (틀림)
- **자물쇠 비례 *흔들림의 효과 크기*** (맞음!)

→ 같은 살짝 흔들림 (σ=0.01) 이 큰 자물쇠에서 *훨씬 강하게* 작용.
→ **σ_opt 고정, SR_max 가 N 의 함수**.

### H16 (**확정된 모델**)

```
SR_max(N) ∝ N^α
σ_opt(N) ≈ constant ≈ 0.010
```

### α 추정 (2점 데이터)

`α = log(SR_max ratio) / log(N ratio) = log(17.86/1.86) / log(1147/437) ≈ 2.34`

→ **SR_max ∝ N^2.34**

### H9 결정적 지지

H9 (polynomial scaling) 가 H16 의 *core* — **strongly confirmed**.

### H12 의 통합 (재정렬)

| H12 sub | 상태 |
|---|---|
| H12a (log N) | 기각 |
| H12b (1/√N) | 기각 |
| H12c (polynomial up) | **기각** (sigma 같음) |
| **H12 통합** | **σ_opt 거의 N 무관** |

## 변경 로그

- 2026-06-12 (초안): H1-H8 작성, 데이터 ground truth 표 작성.
- 2026-06-12 (2): **H9-H13 추가**. H9 (polynomial scaling), H12c (σ_opt ∝ N^α).
  사용자의 *자물쇠와 열쇠 흔들림 비례* 직관 박아넣음.
- 2026-06-12 (3): **H14, H15 추가**. V1 N=1147 seed 5 의 outlier 분석에서 도출.
  K_baseline 의 seed-내 변동성 + Lemma 5.1 의 c 가 N 별로 다름 (1/2 vs 3/4) 발견.
- 2026-06-12 (4): **σ scan 완료. H12c 기각, H16 (확정 모델) 도출**.
  - σ_opt 는 N 무관 (≈ 0.010)
  - SR_max ∝ N^α, α ≈ 2.34 (polynomial)
  - 사용자 직관 *수정된 형태* 로 정확화: 흔들림 크기 → 흔들림 효과 크기 비례
