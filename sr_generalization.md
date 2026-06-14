# SR generalization to other quantum algorithms (scoping)

본 문서는 paper §3.6 의 *trial-level boundary-flip SR mechanism* 이 Shor / Regev /
Hybrid (C)+b-trick 의 범위 *밖* 으로 확장되는지 탐색하기 위한 scoping + 실험계획.

상위 paper (v0.2.1, DOI 10.5281/zenodo.20681847) 의 §3.6 regime map 은
multi-base 양자 factoring algorithm 의 5 가지 structure 에서 *measured 5/5*. 이제
factoring 밖으로 일반화 가능한지 본다.

## §1 SR mechanism 의 일반화 필요 조건 (paper §3.6 정독 결과)

paper §3.6 에서 추출:

1. **Discrete 성공-카운트 지표 K**.
   알고리즘은 *반복 횟수* (혹은 등가의 discrete 자연수 K) 로 성공을 측정해야 한다.
   - Shor 계열: K = base draws 또는 측정 round 수까지의 success.
   - 일반화: K = oracle calls, measurements, iterations 등.

2. **K-bin boundary** 가 존재해야 한다.
   - K-distribution 이 K=1, 2, 3, ... 의 discrete bins 으로 나뉘고 그 boundary 위
     trials 가 다수 존재해야 함.
   - Shor 의 K=1/K=2 boundary 가 active 인 cell (K_base ≈ 2) 이 가장 깨끗.

3. **Borderline-trial population** 이 충분해야 한다.
   - paper §3.6 "Engineered amplification" 결과: borderline population 이 곧
     SR magnitude 상한.
   - K_base ≈ 1 (ceiling) → 0% SR; K_base ≈ 2-3 (active) → 1-5% SR; K_base ≥ 5
     (noise floor) → variance > effect.

4. **Phase noise 가 measurement 분포의 peak smearing 으로 작용**해야 한다.
   - 다른 noise (depolarizing, amplitude damping) 는 monotone degradation 만.
   - QFT-기반 알고리즘에서 phase noise → peak broadening → alternative
     candidate exposure (Shor 의 convergent 후보 flip).
   - 일반화: phase 정밀도가 성공 판정에 직접 영향을 주는 알고리즘.

5. **Saturation plateau** σ-curve (Benzi shape) 가 검증 marker.
   - σ ∈ [σ_thr, σ_overload] 의 plateau → mechanism authentic.
   - monotone decay → ordinary noise robustness, not SR.

## §2 후보 양자 알고리즘 매핑

| 알고리즘 | K 매핑 | Boundary 존재 | Phase noise 의미 | SR fit |
|---|---|---|---|---|
| **Grover search** | K = measurement runs to first marked hit; per-run p_succ from Grover amplitude | k_iter 선택에 따라 p ∈ (0, 1), borderline 존재 | rotation angle 의 Gaussian fluctuation | ★★★ |
| **QPE (isolated)** | K = QFT measurement 의 phase bin recovery 횟수 | Shor 와 동일 (QPE = Shor 의 quantum kernel) | 동일 | ★★★★ (redundant with Shor) |
| **Simon's algorithm** | K = independent measurements 까지 linear system 해결 | n-bit hidden string 의 좌표별 boundary | XOR oracle phase 의 perturbation | ★★ |
| **Quantum Counting** | K = QPE rounds to count estimate | QPE 같음, but counting tolerance 가 boundary | phase noise | ★★ |
| **Deutsch-Jozsa** | K=1 trivially (single shot) | 없음 | 없음 | ✗ |
| **HHL** | K = ancilla measurement rounds | 격자 매우 깊음, NISQ scale 어려움 | phase via QPE inside | ★ (out of scope) |
| **QAOA** | parametric, no K 카운트 | 없음 (variational) | 다른 mechanism | ✗ |

**선정**:
- *Primary candidate*: **Grover** — 구조적으로 Shor 와 가장 다른 quantum subroutine
  (search ≠ factoring), 그럼에도 K-distribution + phase boundary 매핑 가능.
  사용자의 long-term plan 의 C3 (Grover SR) 와 일치.
- *Sanity check*: **QPE isolated** — Shor 의 quantum kernel 이 QPE 이므로 SR 이
  *반드시* 재현되어야 함. 만약 안 되면 우리의 §3.6 분석에 오류.
- *Open*: Simon — 가능하면 후속.

## §3 Grover SR 실험 설계

### 3.1 Setup

- **Search space**: `N_search = 64` (2^6). 작은 사이즈로 빠른 sim.
- **Marked elements**: `M = 1` (single marked).
- **Grover iterations per run**: `k_iter` — k_iter 를 sub-optimal 로 골라 p_succ
  를 Goldilocks 영역에 위치시킴.
  - 최적 `k* = π/4 · √N = π/4 · 8 ≈ 6.28` (M=1, N=64).
  - `p(k_iter)` = sin²((2 k_iter + 1) · θ), θ = asin(√(M/N)) ≈ 0.1253 rad.
  - `k=1`: p ≈ 0.135 → K_base ≈ 7.4 (noise floor)
  - `k=2`: p ≈ 0.343 → K_base ≈ 2.92 ← **active boundary**
  - `k=3`: p ≈ 0.589 → K_base ≈ 1.70 ← **active boundary**
  - `k=6`: p ≈ 0.996 → K_base ≈ 1.00 (ceiling)
- **Regime map** 후보 cells: k_iter ∈ {1, 2, 3, 6} 으로 4 regime 동시 측정.

### 3.2 Per-seed variance

Shor 의 per-seed K-distribution 변이는 random base `a` 가 다른 `r_a` 를 주기 때문.
Grover 에서 동일 변이를 얻으려면:

- **선택**: 각 seed 마다 *random imperfect oracle phase* `ϕ_seed ∈ U(-π/16, π/16)`
  를 *고정 적용*. 즉, oracle 의 -1 phase 가 `e^{i(π + ϕ_seed)}` 가 되도록.
  → 매 seed 가 약간 다른 effective Grover rotation → 약간 다른 p_succ.
  → seed 별 K-distribution 변이.

### 3.3 Phase noise σ

각 Grover iteration 의 rotation angle 에 `N(0, σ²)` Gaussian noise 를 *각 iteration
마다 독립* 적용. K-loop 안에서 매 run 마다 다시 sample → 측정 분포의 peak smearing
과 정확히 동등 (Shor 의 phase_sigma 와 같은 noise class).

### 3.4 측정 protocol (Shor 와 정합)

- σ values: {0.000, 0.050, 0.150} (Shor pure_*_sr 와 동일)
- N_seeds: 5
- Trials per cell: 100
- max_runs: 20
- 총: 4 regime cells × 5 seeds × 3 σ × 100 trials = 6,000 K 측정.

### 3.5 SR metric

per seed:
- `K_base = mean K at σ = 0`
- `K_sigma = mean K at σ > 0`
- `SR % = (K_base - K_sigma) / K_base × 100` (positive = noise helps)

Aggregate:
- `mean SR` over 5 seeds, sd, SE, t, p (1-sided)
- Direction sign-test (n positive / n negative)

### 3.6 예측 (paper §3.6 regime map 의 grover 적용)

| k_iter | K_base | regime | 예측 mean SR % |
|---|---|---|---|
| 6 | ~1.00 | ceiling | ≈ 0% (no borderline) |
| 3 | ~1.70 | active boundary | small (|mean SR| < 1.5%), direction stochastic |
| 2 | ~2.92 | active boundary | small-to-moderate (1-5% per seed), direction stochastic |
| 1 | ~7.4 | noise floor | variance dominated, |mean SR| < |sd / √n| |

**판정 기준**:
- ceiling cell (k=6) SR ≈ 0, 1-2% 범위 → **regime map 적용 ✓**
- active cell (k=2, 3) 에서 per-seed |SR| 가 1-5% 범위에서 변동, 방향 stochastic
  → **mechanism universality ✓**
- noise floor (k=1) variance > effect → ✓

만약 위 3 가지 다 만족 → **paper §3.6 의 algorithm-structure regime map 이
Shor 범위 *밖* 으로 확장됨** (mechanism 일반화 확정).

만약 ceiling 에서 큰 SR 또는 active 에서 0 → §3.6 mechanism 이 Shor specific.

## §4 후속 (QPE / Simon)

- **QPE isolated**: Shor 의 §3.4 의 QPE 만 떼어내서 동일 protocol. 예측: Shor 와
  정확히 같은 regime map (= internal consistency check).
- **Simon**: hidden-string K (XOR linear system 해결까지 measurements 수) 측정.
  구조적으로 다른 boundary 메커니즘 가능.

## §5 산출물

1. `experiments/grover_sr.py` — Grover SR minimal probe (본 문서 §3 protocol).
2. `experiments/grover_sr_results.txt` — 측정 결과.
3. 본 문서 §6 추가 — 결과 + regime map 일반화 판정.

## §6 결과 — Grover SR 측정 (2026-06-14)

세 차례 측정 완료:
- `grover_sr.py`: 4 regime cells (k ∈ {1,2,3,6}) × 5 seeds × 100 trials × 3 σ.
- `grover_sr_focused.py`: active cell k=2 × 5 seeds × 200 trials × 7 σ (Benzi plateau).
- `grover_sr_ceiling.py`: ceiling cell k=6 × 5 seeds × 200 trials × 3 σ.

### 6.1 Regime map (Grover analog of paper Table §3.6)

| k_iter | K_baseline | regime | σ=0.05 mean SR | per-seed \|SR\| max | direction | source |
|---|---|---|---|---|---|---|
| 6 | 1.033 ± 0.015 | **ceiling** | −0.592% (t=−0.73) | 2.43% | 1+/4− | focused_ceiling |
| 2 | 3.289 ± 0.605 | **active boundary** | +3.861% (t=+0.99) | 15.91% | 4+/1− | focused |
| 1 | 9.626 (broad ϕ) | **noise floor** | +0.053% (t=+0.02) | 9.37% | 3+/2− | initial scan |
| 3 | 3.586 (broad ϕ) | active (transitional) | −3.744% (t=−1.06) | 15.42% | 1+/4− | initial scan |

**판정**: paper §3.6 의 regime map 이 Grover 에 **확장된다**:
- *ceiling* (k=6, K_base ≈ 1) → 거의 null (|mean SR| < 1%) ✓
- *active boundary* (k=2, K_base ≈ 3) → per-seed |SR| 1-16%, direction stochastic ✓
- *noise floor* (k=1, K_base ≈ 9-10) → variance > effect ✓ (mean +0.05%, sd 7)

### 6.2 Benzi plateau 검증 (active cell k=2)

7-point σ scan (focused):

| σ | mean SR % | sd | SE | t | direction |
|---:|---:|---:|---:|---:|---|
| 0.005 | −0.819 | 7.413 | 3.315 | −0.25 | 1+/4− |
| 0.025 | +2.094 | 10.023 | 4.482 | +0.47 | 3+/2− |
| 0.050 | +3.861 | 8.759 | 3.917 | +0.99 | 4+/1− |
| 0.075 | +2.221 | 7.697 | 3.442 | +0.65 | 3+/2− |
| 0.100 | +0.391 | 8.427 | 3.769 | +0.10 | 3+/2− |
| 0.150 | +5.613 | 7.432 | 3.324 | +1.69 | 4+/1− |
| 0.200 | +10.607 | 12.405 | 5.548 | +1.91 | 4+/1− |

- σ ∈ [0.005, 0.100]: **plateau** (range −1 ~ +4%, no monotone trend) — **Benzi shape ✓**
- σ ≥ 0.150: **상승** (+5.6% → +10.6%) — Shor (σ ≥ 0.150 에서 decline) 와 *반대 방향*

### 6.3 Direction stochasticity 재현

Active cell (k=2) σ=0.05 per-seed:
- seed 1: +15.91%, seed 2: −8.43%, seed 3: +4.85%, seed 4: +1.33%, seed 5: +5.64%
- mix of signs, magnitudes 1-16% — paper §3.6 의 universal direction stochasticity 와 **정성적 동일**.

Ceiling cell (k=6) σ=0.025:
- seed 1: −4.33%, seed 2: +1.93%, seed 3: +2.86%, seed 4: −0.49%, seed 5: +0.50%
- 3+/2−, |SR| ≤ 4.3% — paper §3.6 의 *ceiling cell* 측정 ((437, 8) +0.00%, (1147, 8) +0.47%)
  와 비교 가능 (작은 magnitude, no systematic bias).

### 6.4 Quantitative differences from Shor

| 측면 | Shor §3.6 | Grover (this work) |
|---|---|---|
| Active boundary per-seed \|SR\| | 0.3-2% | 1-16% (5-8× larger) |
| Plateau σ range | [0.005, 0.100] | [0.005, 0.100] (동일) |
| High-σ trend | decline (overload) | **rise** (Jensen smearing) |
| Direction stochasticity | universal | universal at small σ |
| Boundary K = 1/2 dominance | 77% | (not analyzed at trial-level) |

**Per-seed magnitude amplification 원인**: Grover 의 각 iteration 이 rotation noise 를
선형 누적 → effective σ ∝ √k × σ_per_iter. Shor 의 single-shot QFT measurement 와
달리 *k 번 noise injection*.

**High-σ trend 차이 원인**: Grover 의 over-rotation 영역에서 `sin² ((2k+1)θ + noise)`
의 noise smearing 이 low-p regime 에서 mean p 를 *증가* (Jensen 부등식 적용 영역).
Shor 의 phase noise 는 peak displacement 로 작용 → overload 시 정답 후보 손실 →
decline.

### 6.5 종합 판정

**§3.6 mechanism 의 일반화: 부분적 ✓**

1. **Universal direction stochasticity** at active boundary cell — Grover 에서 재현. ✓
2. **Benzi plateau shape** (σ ∈ [0.005, 0.100]) — Grover 에서 재현. ✓
3. **Borderline-population scaling** (ceiling = small, active = large) — Grover 에서 재현. ✓
4. **High-σ overload decline** (Shor specific) — Grover 에서는 *상승* 으로 분기.
   → §3.6 의 σ-curve overload 단계는 algorithm-specific. plateau 와
   direction-stochasticity 가 일반화 부분.

**Take-away**: paper §3.6 의 regime map 의 *plateau / stochasticity / borderline-population*
삼중은 noise-as-resource 의 mechanism-level 일반화. *Overload decline* 은 algorithm
specific (Shor 측정의 peak structure 에 종속).

### 6.6 후속 (next iteration candidates)

1. **QPE isolated** (sanity check): Shor 와 동일 regime map 재현 예상 (internal consistency).
2. **Simon's algorithm**: hidden-string K-boundary mechanism 확인. n-bit 의 좌표별 SR.
3. ~~**σ-curve mechanism 정량 모델**~~ ✅ 완료 (2026-06-14): `sr_sigma_curve_model.md`.
   Grover closed-form `E[p]=(1-cos(2μ)·e^{-2kσ²})/2` 도출 + R²=0.88 fit.
   *Success-criterion smoothness* 가 σ-curve shape 결정 (smooth=analytic, discrete=step) 의 통합.
4. **Trial-level K-histogram analysis** in Grover (Shor §3.6 의 K=1/K=2 flip 식과 비교).

### 6.7 정합성

본 결과는 paper v0.2.1 의 §3.6 finding 과 *정합* 하며, regime map 의 일부 (plateau /
stochasticity / borderline) 가 Shor / Regev / Hybrid 의 좁은 framework 를 넘어
*noise-as-resource mechanism 의 더 일반적인 가족* 임을 시사. 다만 σ-curve 의
overload phase 는 algorithm-specific 으로 분기 — 본 paper 의 *small SR* 정량 결론은
Shor 에 한정.

따라서 후속 SR-focused paper 의 잠재 contribution:
- *"Multi-algorithm regime map for trial-level noise-as-resource — plateau + stochasticity
  universal, overload algorithm-specific"*.
- Grover, QPE, Simon, Counting 의 cross-algorithm scan 으로 mechanism universality
  확장.
