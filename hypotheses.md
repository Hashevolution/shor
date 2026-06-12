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

## 변경 로그

- 2026-06-12 (초안): H1-H8 작성, 데이터 ground truth 표 작성.
