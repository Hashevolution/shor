# Roadmap: Workshop note → Conference paper

현재 `paper.md` 는 workshop note 등급. 본 문서는 conference 등급으로 격상하기 위한
다단계 로드맵과, Phase 0 (E 스코핑) 의 결론을 담는다.

## 전체 단계

| Phase | 내용 | 기간 (예상) | 산출물 |
|---|---|---|---|
| 0 | E 스코핑 — Regev 2023 와의 통합 가능성 | (완료) | 본 문서 §1 |
| 1 | A — 정량적 K_λ 정리 (반소수 N=pq) | 3-4 주 | 정리 + 증명 + 실험 검증 |
| 2 | D — noise→covered 점근식 | 2-3 주 | 정리 + curve fit |
| 3 | B — HSP/이산로그 확장 | 4-6 주 | 정리 또는 음수결과 |
| 4 | C — hardware 데모 (IBM Q) | 1-2 주 | 실측 데이터 |
| 5 | E — Regev 통합 본실행 (Phase 0 결과 따라) | 미정 | 통합 결과 또는 음수결과 |

---

## §1 Phase 0 — E (Regev 2023 다리) 스코핑

**질문**: (C)-식 다중 base lcm 후처리가 Regev 2023 의 격자 후처리와 결합 가능한가?
가능하다면 무엇이 얻어지는가?

### 1.1 Regev 2023 알고리즘 구조 요약

Regev (arXiv:2308.06572) 의 핵심 차이 (Shor 와 비교):

- Shor: 단일 base `a`, 위수 `r_a` 회수, 후처리 = `k/Q` 의 연분수
- Regev: `d ≈ √(log N)` 개 base `a_1, …, a_d` 를 **동시에** 사용, 단일 QFT 측정에서 격자점 (lattice point) 회수, 후처리 = LLL 격자 환원으로 joint period 회수

회로 측면 이득:
- gates: O(n²) (Shor) → **O(n^{3/2})** (Regev)
- 측정 횟수: O(log r) → **O(√(log N))**

대가:
- 보조비트 더 많음, ancilla 구조 다름
- 후처리가 **격자 환원 (LLL)** — `k/Q` 연분수보다 무거움
- 측정 노이즈가 격자점에 누적되므로 **noise tolerance 가 좁다** (격자 환원은 정확한 점이 필요)

### 1.2 (C) 의 핵심 자산 재확인

(C) 는 측정 후처리 (classical) 단계의 다음 두 성질을 보유:
1. **다중 base 누적**: `L = lcm(r_{a_1}, r_{a_2}, …)`. 매번 새 정보 누적.
2. **노이즈 불변성**: `r_a | L` 이면 측정 분포 무관 (정리 1).

Regev 와 비교:
- Regev 는 다중 base 를 *병렬* 로 처리 (한 측정에 묶음); (C) 는 *순차* 로 누적
- Regev 의 격자 환원은 노이즈에 민감; (C) 의 divisor search 는 노이즈 불변

### 1.3 통합 가능성 — 4가지 후보

#### 후보 A. 직렬 합성 (Regev → (C))

가장 명백한 합성: Regev 로 초기 `L_0 = lcm(r_{a_1}, …, r_{a_d})` 를 한 번에 회수 →
이후 base 는 (C) 의 fast path (`a^{L_0} ≡ 1` 검사) 또는 slow path (측정 + 분모/약수 후보).

**플라우저블 ★★★**. 단순 합성. 동작 보장.

**의미** — 그러나 작음. Regev 한 번이면 인수분해 끝. (C) 의 *후속 사용* 가 의미를
가지는 시나리오는 한정적 (같은 N 에 대해 여러 base 의 위수를 다 알아야 하는 cryptanalysis
같은 경우; 표준 인수분해 위협 모델에서는 한 번이면 충분).

→ paper 한 줄 코멘트 (`Regev + (C)` 합성 가능) 정도. 큰 결과 아님.

---

#### 후보 B. Regev 의 격자 측정에서 (C)-식 누적 추출

Regev 가 격자에서 회수하는 것은 *각 base 의 위수* 가 아니라 *joint period vector*.
하지만 그 벡터에서 lcm 정보가 추출 가능한가?

좀 더 구체적: Regev 의 measurement 가 `(k_1, …, k_d) ≈ (j_1 Q / r_{a_1}, …, j_d Q / r_{a_d})`
로 보이도록 framing 가능. 그러면 각 좌표에 (C) 의 연분수+divisor 후처리를 *독립적으로* 적용 가능.

**플라우저블 ★★**. Regev 의 격자 framing 이 좌표별 (C) 처리와 호환되는지 확인 필요.
호환되면 → Regev 의 격자 환원 일부를 (C) 로 *대체* 가능.

**의미** — 중간. 노이즈 환경에서 Regev 의 약점 (격자 환원 실패) 을 (C) 로 보완.
"노이즈 robust Regev variant" 같은 결과 가능.

---

#### 후보 C. (C) 의 격자 일반화

(C) 의 본질을 "베이스 a → 분류기 a^d ≡ 1 + lcm 누적" 으로 추상화. 격자 setting 에서
동일 구조를 정의할 수 있는가?

- 격자 측정 → 격자 sub-lattice 의 누적 `Λ`
- 새 측정 → `Λ` 안에 들어가면 (`Λ` 가 v 를 cover) 노이즈 무관하게 회수

**플라우저블 ★**. 격자 sub-lattice 구조와 lcm 구조 사이의 대응이 자연스럽지 않음.
다만 성공시 → "lattice 다중-base 후처리의 noise-invariance 정리" 같은 새 정리.

**의미** — 성공시 큼 (격자 framework 의 새 후처리 기법). 실패 확률도 큼.

---

#### 후보 D. Hybrid noise-aware 알고리즘

실용 알고리즘:
- 노이즈 추정값이 낮으면 → Regev (회로 효율적)
- 노이즈 추정값이 높으면 → (C) 다중 base (노이즈 robust)

**플라우저블 ★★★**. 동작 보장. 엔지니어링.

**의미** — paper 의 contribution 으로는 작음 (단순 분기). 다만 hardware 실측 (Phase 4)
과 결합하면 "현실적 NISQ 알고리즘 선택 가이드" 라는 응용 contribution 가능.

---

### 1.4 종합 판정

| 후보 | 플라우저블 | 임팩트 | 우선순위 |
|---|---|---|---|
| A. 직렬 합성 | ★★★ | ★ | (이미 paper 코멘트로 충분) |
| B. Regev 격자에서 (C) 추출 | ★★ | ★★ | **1순위** — 본 phase 의 표적 |
| C. (C) 의 격자 일반화 | ★ | ★★★ | 2순위 — 욕심 옵션 |
| D. Hybrid | ★★★ | ★ | Phase 4 와 결합 |

**Phase 5 의 본실행 표적**: 후보 B. 우선 (C) 를 격자 측정에서 좌표별로 적용 가능한지
formalize. 가능하면 → noise-robust Regev variant. 작업 중 후보 C 의 단서가 보이면 확장.

**Go/no-go 결정**: **Go**. 후보 B 가 viable 으로 보임. 다만 실행 전 Regev 2023 의 정확한
격자 framing 을 읽어야 (인터넷 접근시) — 본 스코핑은 사전 지식 기반의 상위 평가.

### 1.5 E 본실행 (Phase 5) 의 첫 단계

1. Regev 2023 paper (arXiv:2308.06572) 의 §3-4 (격자 환원 후처리) 정독.
2. 격자 measurement 의 좌표별 `k_i` 가 `(j_i Q / r_{a_i})` 형태로 해석되는지 확인.
3. 가능하면 좌표별 (C) 후처리 시뮬레이션 (numpy 격자 환원 또는 LLL 라이브러리).
4. 노이즈 시나리오에서 Regev 단독 vs Regev + 좌표별 (C) 의 성공률 비교.

---

## §2 Phase 1 — A: 정량적 K_λ 정리 (도출 완료)

### 2.1 정리 진술

반소수 `N = pq` (`p, q` 서로 다른 홀소수) 와 `(Z/N)*` 의 균등 무작위 `K` 개 원소
`a_1, …, a_K` 에 대해, `L_K := lcm(ord_N(a_1), …, ord_N(a_K))`. 각 prime `ℓ | λ(N)` 에 대해

`s_ℓ := |{ξ ∈ {p, q} : v_ℓ(ξ-1) = v_ℓ(λ(N))}| ∈ {1, 2}`

를 정의한다 (`v_ℓ` 는 `ℓ`-adic valuation).

**정리 2 (K_λ 분포).**

> (a) **꼬리 상한 (sharp).** `K ≥ 1` 에 대해
> `P[L_K < λ(N)] ≤ Σ_{ℓ | λ(N)} ℓ^{−K · s_ℓ}`.
>
> (b) **꼬리 상한 (simple).**
> `P[L_K < λ(N)] ≤ ω(λ(N)) · 2^{−K}`.
>
> (c) **평균.** `K_λ := min{K : L_K = λ(N)}` 에 대해
> `E[K_λ] ≤ 1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} − 1)`.
>
> (d) **고확률.** 임의 `ε ∈ (0, 1)` 에 대해
> `K_λ ≤ ⌈log₂(ω(λ(N))/ε)⌉ + 1` 가 확률 `≥ 1 - ε`.

**전형 점근**: ω(λ(N)) ≤ ω(p-1) + ω(q-1) = O(log log N) (Hardy-Ramanujan), 따라서
`E[K_λ] = O(log log log N)` (전형 반소수).

### 2.2 증명

**Step 1 (군 분해).** N = pq 에서 CRT 로 `(Z/N)* ≅ (Z/p)* × (Z/q)* ≅ C_{p-1} × C_{q-1}`.
균등 a ∈ (Z/N)* 은 (x, y) ∈ C_{p-1} × C_{q-1} 의 균등 쌍에 대응하며 두 좌표는 **독립**.

**Step 2 (cyclic group 위수 분포).** Cyclic group `C_n` 에 대해 `v := v_ℓ(n)`, 균등 무작위
`x ∈ C_n` 와 `1 ≤ k ≤ v`:

`P[v_ℓ(ord(x)) ≥ k] = 1 − 1/ℓ^{v−k+1}`.

*증명.* `C_n` 의 부분군은 `n` 의 약수 `m` 에 대응하는 `H_m = {x : ord(x) | m}`, `|H_m| = m`.
조건 `v_ℓ(ord(x)) < k` 는 `ord(x) | m` 인 `m` 으로 `v_ℓ(m) ≤ k - 1`. 그러한 `m | n` 의 최대값은
`m* = n · ℓ^{k-1-v}` (l 의 지수만 v→k-1 로 낮춤). 따라서 `{x : v_ℓ(ord(x)) < k} = H_{m*}`,
`|H_{m*}| = m* = n/ℓ^{v-k+1}`. `P = 1 − (n/ℓ^{v-k+1})/n = 1 − 1/ℓ^{v-k+1}`. ∎

(수치 검증: `python -c "..."` 로 N=120, 96, 1000 에서 정확히 일치 확인.)

**Step 3 (단일 base 당 prime ℓ cover 확률).** `v := v_ℓ(λ(N)) = max(v_ℓ(p-1), v_ℓ(q-1))`,
`v_p := v_ℓ(p-1)`, `v_q := v_ℓ(q-1)`. 균등 base `a = (x, y)` 에 대해:

`ord_N(a) = lcm(ord(x), ord(y))`,
`v_ℓ(ord_N(a)) = max(v_ℓ(ord(x)), v_ℓ(ord(y)))`.

이 값이 `< v` 일 확률:
- `v_p < v`: x 는 항상 `v_ℓ(ord(x)) ≤ v_p < v`. 조건은 자동 만족. P[x 의 조건] = 1.
- `v_p = v`: Step 2 (k=v) 로 `P[v_ℓ(ord(x)) < v] = 1/ℓ`.
- 같은 분석 q.

따라서 (독립성):

`P[v_ℓ(ord_N(a)) < v] = (1/ℓ)^{s_ℓ}` (where s_ℓ ∈ {1, 2} as defined).

**Step 4 (K-fold tail bound).** `K` 개 독립 base 에서 `ℓ`-성분이 한 번도 cover 되지 않을 확률:

`P[L_K 의 ℓ-성분 < v_ℓ(λ(N))] = (1/ℓ)^{K · s_ℓ} = ℓ^{−K s_ℓ}`.

`L_K = λ(N)` ⇔ 모든 `ℓ | λ(N)` 의 ℓ-성분이 cover. 여집합에 union bound:

`P[L_K < λ(N)] ≤ Σ_{ℓ | λ(N)} ℓ^{−K s_ℓ}`. → **(a) 증명.**

`s_ℓ ≥ 1`, `ℓ ≥ 2` 로 `ℓ^{−K s_ℓ} ≤ 2^{−K}`, 합은 `ω(λ(N)) · 2^{−K}`. → **(b) 증명.**

**Step 5 (Expectation).**

`E[K_λ] = Σ_{K=0}^∞ P[K_λ > K] = Σ_{K=0}^∞ P[L_K < λ(N)]`.

`K=0` 에서 `L_0 = 1 < λ(N)`, P = 1. `K ≥ 1` 에서 (a) 적용:

`E[K_λ] ≤ 1 + Σ_{K=1}^∞ Σ_ℓ ℓ^{−K s_ℓ}` (interchange)
`= 1 + Σ_ℓ ℓ^{−s_ℓ} / (1 − ℓ^{−s_ℓ})`
`= 1 + Σ_ℓ 1/(ℓ^{s_ℓ} − 1)`. → **(c) 증명.**

**Step 6 (고확률).** (b) 에서 `P[L_K < λ(N)] ≤ ω · 2^{−K} ≤ ε` 풀면
`K ≥ log₂(ω/ε)`. ⌈⌉ 처리. → **(d) 증명.** ∎

### 2.3 선행연구 비교

- **Pomerance et al. (2017)**: 일반 abelian `G` 에 대해 `e(G) ≤ d + 2.752` (생성에 필요한
  기댓값). semiprime `(Z/N)* = C_{p-1} × C_{q-1}` 에서 각 Sylow 가 2-cyclic 이하 → `d ≤ 2`
  → `e ≤ 4.75`. *생성* 은 *exponent 도달* 보다 강하므로 `E[K_λ] ≤ e(G) ≤ 4.75` 가
  Pomerance et al. 의 따름정리. 정리 2(c) 는 이를 다음 두 측면에서 강화:
  - (Z/N)\* 의 구조를 사용한 ℓ-Sylow-별 분석.
  - "exponent 도달" 의 의미만큼만 비용 부담.

- **Knill (1995), Bach-Shallit (1996)**: 다중 base lcm 이 λ(N) 으로 수렴하는 것을 언급하지만
  **정량적 K 분포 분석은 없음**.

- **Carmichael paper (2021, arXiv 2111.02488)**: λ(N) 추정 알고리즘이지만 base 수의 분포
  에 대한 명시적 정리 없음. Algorithm 1 의 분석은 high-prob bound 에 가깝지만 본 정리 형태와
  다름.

**위치**: 정리 2 는 Phase 1 의 새 기여. 이전 작업에 없던 정량적 상한.

### 2.4 실험 검증 계획

- `experiments/k_lambda_dist.py` 를 확장: 각 N 에 대해 정리 2(c) 의 예측값 `1 + Σ 1/(ℓ^{s_ℓ}-1)`
  계산, 경험 평균과 비교 컬럼 추가.
- 큰 N (10⁴ 이상) 에서 ω(λ) 분포 수집 → (d) 형태의 고확률 보증 확인.

### 2.5 다음 세션 시작점

- §2.4 의 확장 스크립트 구현 + 실행.
- paper.md/tex 의 §3.5 에 정리 2 삽입 — Theorem 1 (Determinism) 옆에 자매 정리로.
- §5 (Related Work) 에 Pomerance et al. 비교 한 문단 추가.

---

## §3 Phase 2 — D: Noise→covered 점근식 (착수)

### 3.1 문제

`(C)` 알고리즘이 누적하는 `L` 은 매 trial 마다 잠재적으로 확장된다. 노이즈는 매 trial 의
**확장 성공 확률** 만 줄인다. 정리 2 (노이즈-free) 에서 `E[K_λ]` 이 도출됐다 — 노이즈 하
의 일반화를 원함.

### 3.2 핵심 통찰 (sketch)

매 trial 의 사건은 두 단계로 분리:

1. **Sylow 진전 (coupling).** 새 base 가 새 ℓ-Sylow 정보를 가지는가? 노이즈 무관.
   P = 1 - (1/ℓ)^{s_ℓ} (정리 2 Step 2 그대로).
2. **추출 성공 (noise-dependent).** (C) 가 measurement k 로부터 그 정보를 회복하는가?
   이 단계가 noise 의존.

매 trial 의 L 확장 확률 = P_progress · g(η) where g(η) := P[(C) succeeds | r_a ∤ L_before, noise η].

### 3.3 정리 3 (안)

> **정리 3 (noise-adapted K_λ).** 노이즈 모델 M 의 effective extension prob `g_M ∈ (0, 1]`
> 에 대해, 노이즈 하의 `K_λ^M` 은
>
> `E[K_λ^M] ≤ E[K_λ] / g_M  ≤  (1 + Σ 1/(ℓ^{s_ℓ}-1)) / g_M`
>
> P[L_K^M < λ(N)] ≤ ω(λ(N)) · 2^{-K · g_M / log₂ e}.

### 3.4 g_M 의 계산 (모델별)

- depolarizing p: `g = (1-p) · g_0` (g_0 ≈ 0.9 in noise-free regime)
- bias_zero p: 같음
- modexp q: `g = (1-q) · g_0` (구조 파괴 시 회복 0)
- phase σ: g(σ) = ?(peak 분산도 의존). closed form derivation 필요 (gaussian peak smearing).
- amp_damp γ: g(γ) ≈ exp(-γ·t) · g_0 ? 검증 필요.

### 3.5 1차 측정 결과 (2026-06-12)

`experiments/g_eta.py` 로 N=437, L_before=1 에서 5종 노이즈 측정 (500 trials/condition).

`g_0 = 0.38` (noise-free, L=1 의 어려운 조건. 이는 *위수 회수 확률* 이지 L 확장
확률이 아님; 새 base 가 새 Sylow info 를 추가할 확률 α 와 별개).

| 모델 | p/q/σ/γ → g/g_0 거동 |
|---|---|
| depolarizing | p=0.1→0.89, 0.5→0.59, 0.9→0.24. (1-p) 식보다 항상 큼. |
| bias_zero | p=0.1→0.85, 0.5→0.49, 0.9→0.12. (1-p) 식과 거의 일치. |
| modexp_error | q=0.1→0.84, 0.5→0.34, 0.9→0.18. (1-q) 식과 일치-근사. |
| phase_sigma | σ=0.5→0.85, 1.0→0.51, 2.0→0.17. 비선형 (peak 분산). |
| amplitude_damp | γ=0.001→0.44. 매우 민감 (Q=2^18 에서 exp(-γQ) 효과). |

### 3.6 핵심 모델 (2-parameter form)

데이터가 다음 형태로 깔끔히 정리됨:

`g_M(η) ≈ (1 − η_eff(η)) · g_0  +  η_eff(η) · g_unif_M`

- `η_eff(η)`: 모델의 효과적 "destruction rate" (depol/bias/modexp 는 η_eff=η,
  phase/amp_damp 는 모델별 함수).
- `g_unif_M`: 노이즈의 stationary 분포에서 (C) 가 회수할 확률.
  · depol: `g_unif ≈ 0.06` (uniform k 의 lucky rate)
  · bias_zero: `g_unif = 0` (k=0 → 분모=1 만 후보)
  · modexp: `g_unif ≈ 0.075`
  · phase/amp_damp: ~0.05

### 3.7 정리 3 (수정안)

> **정리 3 (algorithm K_λ under noise).** 노이즈 M 의 effective extraction rate `g_M(η)` 가
> 정의됐을 때, (C) 알고리즘이 `L = λ(N)` 에 도달하는 base 수 `K_λ^alg`:
>
> `E[K_λ^alg(η)] ≤ E[K_λ^ideal] / g_M(η)  ≤  (1 + Σ 1/(ℓ^{s_ℓ}-1)) / g_M(η)`
>
> 여기서 `E[K_λ^ideal]` 은 정리 2(c) 의 상한이다.

증명 sketch: 각 trial 의 L-확장 확률 = (Sylow-새-정보 확률 α_s) × (C 회수 확률 g_M).
독립성으로 단계별 기댓값 / g_M.

### 3.8 다음 세션 시작점

- 정리 3 의 엄밀 증명 작성 (coupling: 각 ideal trial 은 평균 1/g_M 개의 algorithm trial 에 대응).
- η_eff 의 모델별 정확한 표현 정리 — depol/bias/modexp 의 (1-η) 식 정당화, phase/amp 의
  근사식 도출.
- N 다양화: N=437 외에 N=1147, 4087 에서 동일 측정 → g_M 의 N-의존성 확인.
- paper.md/tex 의 §3.4 또는 §4.1 에 정리 3 + Phase 2 검증 표 추가.

## §4 진행 로그

- **2026-06-12 (1)**: Phase 0 스코핑 완료. E 의 후보 B 를 Phase 5 표적으로 확정.
  본 문서 §1 추가. Phase 1 (A) 착수 — §2 의 정리 윤곽 작성.
- **2026-06-12 (3)**: **Phase 2 착수.** `g_eta.py` 작성 + N=437 에서 5종 노이즈 측정.
  핵심 발견: g_M(η) ≈ (1-η_eff)·g_0 + η_eff·g_unif_M 의 2-parameter form.
  depol/bias/modexp 는 (1-η) 식 잘 맞음. phase/amp_damp 는 비선형 의존.
  정리 3 (수정안) 작성 — `E[K_λ^alg] ≤ E[K_λ^ideal] / g_M`. 본실행은 다음 세션.

- **2026-06-12 (2)**: **Phase 1 완료.**
  - 정리 2 (K_λ 분포) 완전 도출 — 4 형태 (tail sharp, tail simple, expectation, high-prob).
    Step 2 (cyclic group 위수 분포) 를 N=120, 96, 1000 수치 검증으로 확인.
  - `experiments/k_lambda_dist.py` 확장: 정리 2(c) 상한 계산 + 경험 평균과 비교.
    17 semiprimes × 1,000 trials 에서 **모든 N 에서 상한 통과**.
    N=4087 (ω(λ)=4): 경험 평균 2.26 vs 상한 2.475 vs Pomerance 보장 4.752.
  - `paper.md` / `paper.tex` 통합:
    · §3 → "Main theorems" (3.1 정리 1, 3.2 정리 2, 3.3 joint interpretation).
    · §3.2 에 완전 증명 (5 step) + Pomerance et al. 비교.
    · Abstract, Conclusion, §6 limitations, §7 reproducibility 업데이트.
    · Appendix D 신설 (17행 정리 2 검증표).
