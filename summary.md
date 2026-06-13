# SR 연구 종합 도표 (2026-06-12 ~ 13 누적)

본 문서는 모든 실험 결과 + 가설 평가 종합. 도중 진위 / fluke / 진짜 finding 의
구분 명확화.

## 1. 모든 SR 측정 (대표 N별)

### N=437 (= 19 × 23, λ=198, c=1/2)

| d | K_base | K(σ=.05) | SR % | trials | source |
|---|---|---|---|---|---|
| 1 | 7.277 | 7.220 | +0.78% | 300 | AOP |
| 2 | 4.473 | 4.483 | -0.22% | 300 | AOP |
| 2 | 4.420 | 4.400 | +0.45% | 200 | M1 |
| 3 | 2.707 | 2.673 | +1.23% | 300 | AOP |
| 4 | 2.322 | 2.298 | +1.03% | 500 | sr_confirm |
| 4 | 2.205 | 2.165 | +1.81% | 200 | M1 |
| **4** | **1.9160** | **1.8985** | **+0.91%** | **2000** | **V3 ✓★** |
| 8 | 1.190 | 1.195 | +0.00% | 200 | M1 (ceiling) |

### N=1147 (= 31 × 37, λ=180, c=3/4)

| d | K_base | K(σ=.05) | SR % | trials | source | 신뢰 |
|---|---|---|---|---|---|---|
| 1 | 5.780 | 5.630 | **+2.60%** | 300 | AOP | 단일 |
| 2 | 2.630 | 2.405 | +8.56% | 200 | n_extend | ⚠ |
| 2 | 2.800 | 2.300 | +17.86% | 150 | σ scan | **fluke 확정** |
| 2 | 2.427 | 2.397 | +1.24% | 300 | AOP | 단일 |
| **2** | **~2.78** | **~2.77** | **+0.42%** | **4×1000** | **confirm ✗** |
| 4 | 1.440 | 1.355 | +5.90% | 200 | n_extend | 변동 |
| 8 | 1.055 | 1.050 | +0.47% | 200 | ceiling | |

### N=2491 (= 47 × 53, λ=1196, c=3/4)

| d | K_base | K(σ=.05) | SR % | trials | source |
|---|---|---|---|---|---|
| 1 | ?? | ?? | **진행 중** | 300 | AOP |
| 2 | 2.250 | 2.360 | -4.89% | 200 | n_extend |
| 4 | 1.065 | 1.085 | -1.88% | 200 | ceiling |
| 8 | 1.000 | 1.015 | -1.50% | 200 | floor |

## 2. V Validation 결과

### V1 N=437 d=4 (5 seeds × 300 trials)

| seed | K(σ=0) | K(σ=.05) | SR % |
|---|---|---|---|
| 1 | 1.770 | 1.750 | +1.13% |
| 2 | 1.643 | 1.630 | +0.81% |
| 3 | 1.800 | 1.833 | -1.85% |
| 4 | 1.747 | 1.740 | +0.38% |
| 5 | 1.703 | 1.683 | +1.17% |

**mean = +0.33%, sd = 1.26, t = 0.58 → borderline**

### V1 N=1147 d=2 (5 seeds × 100 trials)

| seed | K(σ=0) | K(σ=.05) | SR % |
|---|---|---|---|
| 1 | 2.840 | 2.770 | +2.46% |
| 2 | 3.270 | 3.030 | +7.34% |
| 3 | 3.130 | 2.930 | +6.39% |
| 4 | 2.580 | 2.260 | +12.40% |
| 5 | 2.110 | 2.220 | -5.21% |

**mean = +4.68%, sd = 6.57, t = 1.59 → borderline (high variance)**

### V3 N=437 d=4 σ scan (2000 trials per σ) ★

| σ | hybrid K | Δ | SR % |
|---|---|---|---|
| 0.000 | 1.9160 | - | baseline |
| 0.025 | 1.8985 | -0.018 | +0.91% |
| 0.050 | 1.8985 | -0.018 | +0.91% |
| 0.075 | 1.8990 | -0.017 | +0.89% |
| 0.100 | 1.9005 | -0.016 | +0.81% |
| 0.200 | 1.9070 | -0.009 | +0.47% |

**sign test: 5/5 negative → p = 1/32 = 3.1% (유의) ✓**

### 17.86% Final Confirm ✗

| seed | K(σ=0) | K(σ=.01) | SR % |
|---|---|---|---|
| 1 | 2.705 | 2.648 | +2.11% |
| 2 | 2.818 | 2.868 | -1.77% |
| 3 | 2.739 | 2.707 | +1.17% |
| 4 | 2.852 | 2.847 | +0.18% |

**mean = +0.42%, sd = 1.66, t = 0.50 → NOT significant**

**→ 17.86% scan 결과 = fluke 확정** ✗

## 3. σ Scan (N=437 d=4, 300 trials per σ)

| σ | K | SR % | 영역 |
|---|---|---|---|
| 0.000 | 2.153 | baseline | |
| 0.010 | 2.113 | +1.86% | σ_opt |
| 0.025 | 2.113 | +1.86% | plateau |
| 0.050 | 2.113 | +1.86% | plateau |
| 0.075 | 2.113 | +1.86% | plateau |
| 0.100 | 2.120 | +1.55% | declining |
| 0.150 | 2.120 | +1.55% | declining |
| 0.200 | 2.120 | +1.55% | declining |
| 0.300 | 2.130 | +1.08% | declining |
| 0.500 | 2.173 | -0.93% | anti-SR |

**σ_opt = 0.010 (가장 작음), peak σ range = [0.010, 0.075]**

## 4. K_baseline vs SR % 산점도

```
K_base │  (N, d)      │ SR %    │ 영역
───────┼──────────────┼─────────┼───────────────
1.000  │ (2491, 8)    │ -1.50%  │ floor
1.055  │ (1147, 8)    │ +0.47%  │ ceiling
1.065  │ (2491, 4)    │ -1.88%  │ ceiling
1.190  │ (437, 8)     │ +0.00%  │ ceiling
1.440  │ (1147, 4)    │ +5.90%  │ surprising (variable)
1.916  │ (437, 4)     │ +0.91%  │ V3 robust ✓
2.250  │ (2491, 2)    │ -4.89%  │ low (N effect)
2.322  │ (437, 4)     │ +1.03%  │ sr_confirm
2.427  │ (1147, 2)    │ +1.24%  │ AOP
2.630  │ (1147, 2)    │ +8.56%  │ n_extend (variable)
2.707  │ (437, 3)     │ +1.23%  │ AOP
2.800  │ (1147, 2)    │ +17.86% │ FLUKE
4.420  │ (437, 2)     │ +0.45%  │ M1
4.473  │ (437, 2)     │ -0.22%  │ AOP
5.780  │ (1147, 1)    │ +2.60%  │ AOP (단일)
7.277  │ (437, 1)     │ +0.78%  │ AOP
```

**관찰**: K_base 와 SR 의 관계는 단순 monotonic 아님. Sweet spot ~2-3 있으나 변동 큼.

## 5. 가설 평가 종합

| 가설 | 설명 | 상태 |
|---|---|---|
| H1 | SR ∝ log Q | 부분 (정량 부족) |
| H2 | SR ∝ d(λ) | 부분 (정량 부족) |
| H3 | v_2 비대칭 | H15 로 흡수 |
| **H4** | **K_baseline sweet spot** | **지지** ✓ |
| H5 | difficulty regime | 정성적 |
| H8 | λ 의 4-divisibility | 모호 |
| H9 | Polynomial SR ∝ N^α | **기각** ✗ |
| H12c | σ_opt ∝ N^α (자물쇠 직관) | **기각** ✗ |
| **H14** | **K_base seed 변동** | **지지** ✓ |
| **H15** | **Lemma 5.1 c = f(N)** | **수학적 확정** ✓ |
| H16 | σ_opt 고정 + SR_max polynomial | **기각** (17.86% fluke) |
| **H17** | **Narrow N window** | **지지** ✓ |
| H18 | d 고정시 SR ∝ N (AOP 시사) | **기각** (N=1147 d=1 multi-seed null) ✗ |
| AOP | Slack → SR 단조 | **기각** (N=1147 d=1 multi-seed -0.53% ± 4.28%) ✗ |
| **H19** | **Goldilocks K_base ≈ 2 = sweet spot** | **단일 cell 지지** (N=437 d=4 V3) |

## 6. 최종 mechanism observation (Plan A')

### 13-seed σ scan at (437, 4): mechanism universality 확정

| 측정 | 값 |
|---|---|
| Cell | (N=437, d=4) |
| Setup | 13 seeds × 200 trials × 12 σ = 31,200 trial-measurements |
| Mean SR (σ=0.050) | **+0.144%** |
| sd (between-seed) | 1.016% |
| SE | 0.282% |
| t | +0.51 |
| **p (1-sided)** | **0.31** (NOT significant) |
| Sign test | 8/13 positive (p=0.29) |

### Mechanism level (★ 진짜 finding)

```
Boundary-flip 분포 (13/13 universal):
  K=1/K=2 boundary: 10 seeds (76.9%)  ← primary
  K=2/K=3 boundary:  2 seeds (15.4%)  ← secondary
  K=3↔K=1 long-jump: 1 seed (7.7%)    ← rare

σ-curve direction asymmetry:
  Positive seeds (8): plateau + decline (high σ 에서 baseline 회귀)
  Negative seeds (5): monotonic worsening (K=1 reservoir 무한)

Sanity 검증: 13/13 histogram Δ = K_mean Δ × n (정확 일치)
```

### Direction 의 *base-set 결정성* 직접 증거

```
seed 3, 4, 13: 모두 K_base = 1.720
  seed 3: SR = -0.872%  (K=1→K=2)
  seed 4: SR = +0.581%  (K=2→K=1)
  seed 13: SR = -0.291% (K=1→K=2)
```

→ **같은 K_mean, 다른 direction** — direction 은 base set 의 *internal K-distribution* 결정.

### Cross-cell verification (기존 데이터)

| (N, d) | K_base | regime | SR % | match |
|---|---|---|---|---|
| (437, 8) | 1.19 | ceiling | +0.00% | ✓ exact |
| (1147, 8) | 1.06 | ceiling | +0.47% | ✓ small |
| (2491, 4) | 1.07 | ceiling | -1.88% | ✓ small |
| (2491, 8) | 1.00 | floor | -1.50% | ✓ noise |
| (1147, 1) | 5.78 | noise floor | -0.53% ± 4.28% | ✓ var > effect |
| (437, 4) | 1.92 | active boundary | +0.144% (13-seed) | ✓ small, no net |
| (1147, 2) | 2.43 | active boundary | +0.42% (4×1000) | ✓ regresses |

→ **모든 regime 의 prediction *일치*** (mechanism universal across regimes).

## 7. 시간순 narrative

```
초기 (sr_confirm, σ scan): 17.86% 발견 → 흥분
+ V validation:           N=437 의 ~1% 진짜 확인
+ N extend N=2491:        SR 사라짐 → polynomial 기각
+ 17.86% confirm:         **fluke 확정** → 진짜 ~0.5%
+ AOP:                    slack 가설 검증, non-monotone 패턴
                          → K_base sweet spot 의 미세 구조 발견
```

## 8. 최종 평가 (Plan A')

### 확정된 mechanism facts

- **Boundary-flip mechanism universal** (13/13 seeds at (437,4))
- **Multi-boundary distribution**: 77% K=1/K=2, 15% K=2/K=3, 8% long-jump
- **σ-curve classical SR shape** (sub-threshold + plateau + decline + overload)
- **σ-curve direction asymmetry** (positive seeds 회귀, negative seeds monotonic)
- **Direction 의 base-set 결정성** (seed 3, 4, 13 같은 K_base 다른 direction)
- **σ_opt ≈ 0.01-0.05** (N 무관)
- **Phase noise specific** (depol/amp 안 함)
- **K=20 trials** (max_runs ceiling failures) 은 mechanism 영향 없음 (dilution 만)

### 통계적 미확정 (Plan A' caveat)

- **Net SR direction**: mean +0.144%, t=0.51, p=0.31 — **NOT significant**
- Sign test 8/13 (p=0.29) — NOT significant
- 30+ seeds 면 작은 net bias 검출 가능성 있음

### 기각된 주장 (6개)

1. Polynomial scaling (H9, H16)
2. "자물쇠 비례 흔들림" 직관 (H12c)
3. 17.86% peak (fluke)
4. 2.60% peak (N=1147 d=1 multi-seed null)
5. AOP "d=1 universal" (H18)
6. **V3 "p=0.03" significance** — within-seed σ correlation 으로 *허위* (★ new retraction)

### Paper 등급 (Plan A')

- 정리 1-5 (folklore + hybrid + Lemma 5.1) = main contribution
- **§3.6 universal trial-level mechanism observation**
  - 13/13 boundary flip universality
  - σ-curve asymmetry
  - Direction base-set 결정성
- **Mid-tier journal 가능 등급**: 워크숍 격상 + mechanism depth
- **RSA 보안에 변화 없음**

## 9. 진행 중 / 다음 후보

### 완료
- σ scan (437, 4) 13 seeds × 200 trials × 12 σ → mechanism universal
- K-histogram backfill + flip analysis (13/13 universal boundary flip)

### 진행 중 (2026-06-13 18:19~)
```
σ scan (1147, 2) compact: 5 seeds × 100 trials × 5 σ (background)
  - σ = {0.000, 0.025, 0.050, 0.150, 0.200}
  - 예상 ~1.6 시간
  - 검증:
    * mechanism universality cross-N
    * K_base=2.43 → boundary distribution 차이 (K=2/K=3 더 활성?)
    * direction at 다른 N (positive vs negative split)
    * σ-curve 동일 shape?
```

### 미진행 후보 (mechanism universality 추가 확인)
```
- (4087, 4) K_base=1.43 — 큰 N, ~5 시간 (1 seed × 200 trials per σ × 5 σ ~ 1 시간 가능)
- (437, 3) K_base=2.71 — boundary distribution shift 검증
- (1147, 3) K_base=1.67 — Goldilocks edge

각 cell 의 marginal value: cross-cell universality 강화. 시간 대비 모두 *marginal*.
```

## 10. 메타 lesson

- **작은 표본의 화려한 결과 = fluke 가능성 매우 높음**
- 가설 세우기 전 high-stat confirm 필수
- "자물쇠 직관" 같은 정량적 모델은 데이터로 직접 검증
- 정직한 평가가 결국 paper 의 신뢰도 결정

## 11. 추가 발견 (2026-06-13 저녁)

### Cross-cell (1147, 2) — High-K rescue 발견 ★

5 seeds × 100 trials × 5 σ at (N=1147, d=2), K_base mean = 2.92:

| seed | K_base | SR (σ=0.05) | dominant flip |
|---|---|---|---|
| 1 | 2.80 | +1.43% | K=2→K=1 (classical) |
| 2 | 3.25 | 0.00% | K=3→K=1 |
| 3 | 3.39 | **+9.44%** ★ | **K=15→K=5 + K=11→K=5 + K=20→K=6** |
| 4 | 2.92 | **+8.56%** ★ | **K=8→K=4 (3 trials!)** |
| 5 | 2.24 | -2.68% | K=1→K=2 (classical neg) |

**Cross-seed**: mean +3.35%, sd 5.37%, t=1.39, p≈0.12 (t-dist) / 0.082 (normal approx) — marginal

**핵심 발견 — High-K rescue**:
- 큰 K (8, 11, 15, 20) trials 가 *moderate K* (4, 5) 로 jump
- Per-seed |SR| 5x amplified vs (437, 4) max
- K_base 큰 cell 에서 *새 mechanism channel* 활성

### Engineered amplification ★

Mild thinned hybrid (ALL convergents + NO (C) augmentation):
```
Full hybrid:    K_base=2.08, per-seed |SR| 0-1.16%
Over-thinned:   K_base=19.87, SR=0% (no borderlines)
Mild thinned:   K_base=2.92, per-seed |SR| 4-5% (★ 5x amplification)
```

→ (C) augmentation 가 *noise buffer* 역할 — 제거 시 mechanism 노출.  
→ Direction 여전히 stochastic (1+, 2-).

### Algorithm-structure regime map (testable conjecture)

| Algorithm | SR | source |
|---|---|---|
| Shor (단일 base) | small (≤1%) | **PREDICTED** |
| Regev (LLL) | negative | **PREDICTED** |
| Hybrid full | +0.14% | **MEASURED** |
| Hybrid mild-thinned | 5x amp | **MEASURED** |
| Hybrid over-thinned | 0 | **MEASURED** |

→ Multi-base + per-coord + no buffer = amplification sweet spot

## 12. Paper / Zenodo 준비 상태 (2026-06-13)

```
paper.md / paper.tex:
  ✓ §3.6 SR observation (honest framing)
  ✓ §3.7 joint interpretation (T5 + SR orthogonal)
  ✓ §7 reproducibility (sr_aop.py 추가)
  ✓ §8 conclusion (T4-5 + SR mention)
  ✓ Abstract 마지막 한 줄 (SR over-claim 없이)

Zenodo 메타데이터:
  ✓ CITATION.cff (GitHub citation 위젯)
  ✓ .zenodo.json (Zenodo 자동 메타데이터)
  ✓ LICENSE (MIT)

남은 작업:
  - N=1147 d=1 multi-seed 결과 → §3.6 단일 셀 confirm 통합
  - 선택: N=4087 d=1 single 측정
  - GitHub release 생성 → Zenodo 자동 sync → DOI 획득
  - arXiv 제출 (endorsement 후)
```

---

진행 로그:
- 2026-06-13 (오전): V + confirm 완료, 17.86% fluke 확정. AOP partial. paper Abstract/§3.6/§3.7/§7/§8 polish 완료.
- 2026-06-13: Zenodo 메타데이터 (CITATION.cff/.zenodo.json/LICENSE) 추가.
- 2026-06-13 (오후): **N=1147 d=1 multi-seed confirm 결과 = noise floor**. 2.60% peak
  도 fluke 확정 (mean -0.53% ± 4.28%). AOP "d=1 universal" 무너짐. paper §3.6
  **Goldilocks 재작성** (K_base ≈ 2 cell 단일 robust, ceiling/Goldilocks/noise-floor
  3 영역 가설). Abstract / §3.7 / §8 정합성 갱신. H18 / AOP 기각, H19 (Goldilocks)
  단일 cell 지지.
- 2026-06-13 (저녁): **σ scan (437, 4) 13 seeds × 200 trials × 12 σ 완료**. K-histogram
  backfill 후 *multi-boundary flip mechanism* 발견 (13/13 universal). Plan A' 확정:
  mechanism universal, net direction 통계적 미확정. Goldilocks 재정정 (boundary trial
  분포 의 K_base 함수). V3 "p=0.03" *허위 significance* 발견 (within-seed σ correlation).
  paper.md/paper.tex §3.6 multi-boundary mechanism 으로 finalize. H19 약화, H22 신규
  (multi-boundary universal).
- 2026-06-13 (밤): (1147, 2) compact scan background 진행. Plan A' cross-cell verification.
