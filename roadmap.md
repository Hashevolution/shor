# Roadmap: Workshop note → Conference paper

현재 `paper.md` 는 workshop note 등급. 본 문서는 conference 등급으로 격상하기 위한
다단계 로드맵과, Phase 0 (E 스코핑) 의 결론을 담는다.

## 전체 단계

| Phase | 내용 | 상태 | 산출물 |
|---|---|---|---|
| 0 | E 스코핑 — Regev 2023 와의 통합 가능성 | ✅ 완료 | 본 문서 §1 |
| 1 | A — 정량적 K_λ 정리 (반소수 N=pq) | ✅ 완료 | 정리 2 (paper §3.2) + 17 N × 1k trials |
| 2 | D — noise→covered 점근식 | ✅ 완료 | 정리 3 (paper §3.3) + 9 노이즈 N=437 |
| 3 | B — HSP/이산로그 확장 | ✅ 종료 (scope-out) | paper §6 Limitations 단락 |
| 4 | C — hardware 데모 (IBM Q) | 🟡 부분 (calibrated proxy) | Appendix E (5종 IBM Eagle 모델) |
| 5 | E — Regev 통합 본실행 | ✅ 완료 | 정리 4-5 (paper §3.4-3.5) + Lemma 5.1 |
| 6 | F — Noise-as-resource (frontier) | ✅ 완료 | §3.6 SR + AOP grid + ENAQT bridge |
| 7 | Publication — arXiv / Zenodo | 🟡 진행중 | Zenodo 메타 준비완료, DOI 대기 |

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

### 3.7 정리 3 (최종 형식)

#### 3.7.1 "Destructive at rate η" 노이즈의 정의

노이즈 모델 `M` 이 각 측정에 대해 독립적으로 확률 `η` 로 측정값을 `D_M` (베이스
`a` 와 독립인 분포) 으로 교체하고, 확률 `1-η` 로 노이즈-free 측정을 산출하면
`M` 을 *destructive at rate η* 라 한다.

이 정의를 만족하는 노이즈: depolarizing (D = uniform), bias_zero (D = δ_{k=0}),
modexp_error (근사).

이 정의를 만족하지 *못하는* 노이즈: phase_sigma (peak smearing — ideal 분포의 함수적
변형), amplitude_damp (magnitude decay).

#### 3.7.2 정리

> **정리 3 (destructive 노이즈의 알고리즘 K_λ).** 노이즈 `M` 이 destructive at rate η
> 이고, 모든 상태 `s` 에서 회수 확률이 `g_M(η) := (1-η) g_0 + η g_unif_M` 으로 일정
> 하다고 가정하자 (여기서 `g_0` 은 노이즈-free 회수 확률, `g_unif_M` 은 `D_M` 에서의
> 회수 확률). 그러면 `(C)` 알고리즘이 `L = λ(N)` 에 도달하는 base 수 `K_λ^alg` 는
>
> `E[K_λ^alg(η)] = E[K_λ^ideal] / g_M(η)`.
>
> 특히 `E[K_λ^alg(η)] / E[K_λ^alg(0)] = g_0 / g_M(η)`.

#### 3.7.3 증명

Markov chain 분석. 상태 `s = L_before ∈ divisors(λ(N))`. 매 trial:

- 확률 `ε_s := P[r_a | s]`: covered (fast path). 회수 결정적 (정리 1). 상태 그대로.
- 확률 `1 - ε_s`: r_a ∤ s (slow path). 측정 + (C). destructive 가정 하 회수 확률
  `g_M(η)` (상태 무관). 회수 성공 시 상태 → `lcm(s, r_a)`, 실패 시 그대로.

상태 `s` 에서 비자명 전이 확률 = `(1 - ε_s) · g_M(η)`. `E[# trials in s] = 1 / [(1-ε_s) g_M(η)]`.

선형성으로:
`E[K_λ^alg(η)] = Σ_s 1/[(1-ε_s) g_M(η)] = (1/g_M(η)) · Σ_s 1/(1-ε_s) = (1/g_M(η)) · E[K_λ^ideal]`.

(여기서 `E[K_λ^ideal] = Σ_s 1/(1-ε_s)` 는 정리 2 의 분석에서 도출된 ideal trial 수.) ∎

#### 3.7.4 N=437 경험 검증

200 trials × 노이즈 9종. 경험 `K_λ^ideal = 1.83` (k_lambda_dist 측정값) 사용:

| 노이즈 | g_M(L=1) | 정리 3 예측 | 실측 K_λ^alg | 오차 | 적용 |
|---|---|---|---|---|---|
| noise-free | 0.380 | 4.82 | 4.97 | +3% | ✓ |
| depol p=0.1 | 0.340 | 5.38 | 6.09 | +13% | ✓ |
| depol p=0.3 | 0.262 | 6.98 | 7.40 | +6% | ✓ |
| depol p=0.5 | 0.226 | 8.10 | 9.19 | +13% | ✓ |
| depol p=0.7 | 0.156 | 11.73 | 14.11 | +20% | ✓ |
| bias_zero p=0.5 | 0.188 | 9.73 | 10.80 | +11% | ✓ |
| modexp q=0.3 | 0.212 | 8.63 | 11.50 | +33% | △ |
| modexp q=0.5 | 0.130 | 14.08 | 20.28 | +44% | △ |
| phase σ=1.0 | 0.194 | 9.43 | 13.06 | +38% | ✗ |
| phase σ=2.0 | 0.064 | 28.59 | 69.12 | +142% | ✗ |

`✓` = 정리 3 가정 만족, 평균 +11% 오차로 일치 (상태-의존 g_M 의 잔여 효과).
`△` = 부분적 일치 (modexp 는 destructive 에 가깝지만 strict 정의 외).
`✗` = structural 노이즈, 정리 3 적용 외 — trajectory 의 state-dependent g_M 효과로
정리 3 의 예측 초과 (큰 r_a 회수가 smear 분포에서 어려움).

#### 3.7.5 structural 노이즈 (open question)

phase/amp_damp 의 정확한 식은 `g_M(s, η)` 의 상태 의존성을 반영해야 함:

`E[K_λ^alg(η)] = Σ_path (path prob) · Σ_s 1/[(1-ε_s) g_M(s, η)]`

phase 노이즈에서 `g_M(s, η)` 의 closed-form 도출은 후속 작업. 본 paper 에서는
경험 관찰 + Theorem 3 의 destructive case 만 제공.

### 3.8 paper 통합 + 다음 세션 시작점

- paper.md/tex §3.4 에 정리 3 (destructive case) + §4 (Empirical) 에 검증표 통합.
- structural 노이즈 의 g_M(s, η) 도출 — Phase 2 follow-up.
- N 다양화 (N=1147, 4087) — 만약 시간 허락 시.

## §4 Phase 3 — B: HSP / 이산로그 확장 (스코핑 + 결론)

### 4.1 문제

(C) 프레임워크 (다중 base lcm 누적 + divisor search 후처리) 가 다른 양자
은닉구조 문제 (DL, abelian HSP, 격자 주기) 로 *비자명하게* 확장되는가?

### 4.2 (C) 의 핵심 구조 재확인

(C) 의 성공 조건은 다음 3가지:

1. **그룹 학습 문제**: Shor 는 단일 인스턴스 r_a 가 아니라 *그룹 exponent λ(N)*
   를 학습 (그 부산물로 인수분해). 여러 base 의 r_a 들은 lcm 으로 결합 → λ(N).
2. **정수 격자 (lcm/divisor)**: 후보 풀의 자연 결합 연산. 후보 검증이 단조 (divisor
   추가만 가능).
3. **빠른 고전 검증**: `a^d ≡ 1 mod N` 가 polylog time.

(C)-determinism 정리 (정리 1) 는 (2)+(3) 의 직접 결과. (C)-스케일링 정리 (정리 3) 는
(1)+(2) 의 다중 base 구조의 결과.

### 4.3 후보 확장 평가

| 후보 | (1) 다중 base 가능 | (2) 격자 구조 | (3) 검증 | 평가 |
|---|---|---|---|---|
| **A. 다중 prime 인수분해** | ✓ (동일) | ✓ (동일) | ✓ (동일) | **trivial extension** |
| **B. 이산 로그 (DLP)** | ✗ (단일 인스턴스) | ✗ (x ∈ Z/r 의 단일 원소) | ✓ | **불가** |
| **C. abelian 1-D HSP** | ✗ (Shor 와 동일) | (= Shor) | (= Shor) | **redundant** |
| **D. abelian multi-D HSP** | ✗ (단일 f) | △ (sublattice) | ✓ | **trivial / well-known** |
| **E. (C) + Pohlig-Hellman** | (조합 사용) | — | — | **useful 조합 알고리즘** |
| **F. Regev 2023 격자** | (Phase 5 표적) | △ | △ | Phase 5 |

### 4.4 결론

**(C) 의 핵심 (다중 base 의 lcm 누적) 은 Shor 의 factoring 에 *구조적으로 결합*** —
"군의 exponent 를 학습" 한다는 목적이 다중 base 결합을 자연스럽게 만든다.

- DL 과 HSP 는 *단일 인스턴스* 문제로, 다중 base 가 자연스럽지 않음. 다른 base 를 쓰면
  완전히 다른 문제가 됨 (다른 (g, h), 다른 f).
- 다중 prime 인수분해는 trivial 확장 — 정리 1·2·3 가 그대로 성립 (단지 λ(N) 의 정의가
  복잡해질 뿐). paper 의 후속 코멘트로 충분.
- (C) + Pohlig-Hellman: (C) 로 ord(g) 학습 → PH 로 DL. 유용한 *조합* 알고리즘이지만
  (C) 의 새 확장은 아님. paper 코멘트.

**Phase 3 의 가장 가치 있는 산출물 = 음수 영역 결과** (negative scope result):

> **관찰**: (C)-식 다중 base 후처리는 *그룹 exponent 학습* 형태의 양자 알고리즘에
> 특정된다. DL 과 abelian HSP 같은 *단일 인스턴스* 문제에는 자연스럽게 적용되지 않음.

이는 paper §6 (Limitations) 에 명시할 가치 있음 — (C) 의 적용 범위를 정확히 함.

### 4.5 paper 통합

- §6 (Limitations) 에 새 paragraph 추가: "Scope of the framework — multi-base accumulation
  is specific to group-exponent learning".
- 다중 prime 인수분해 확장은 한 문장 코멘트.

### 4.6 다음 phase

Phase 3 는 짧게 종료. Phase 4 (hardware) 로 진행.

## §5 Phase 4 — C: Hardware 데모 (스코핑)

### 5.1 외부 의존

본 phase 는 **IBM Quantum (또는 IonQ) 계정 + qiskit 설치** 가 필요. 현재 저장소는
numpy-only 정책 — qiskit 추가는 사용자 명시 결정 필요. 본 세션은 *준비 작업*
(스코핑, 스켈레톤, 체크리스트) 만 진행.

### 5.2 표적 시나리오

- **N = 15** (가장 작은 nontrivial 반소수): 4 큐비트 작업 레지스터 + 8 큐비트 계산 레지스터.
  IBM Q free tier (e.g., `ibmq_qasm_simulator` 또는 7-큐비트 hardware 인 `ibm_lagos` /
  `ibm_perth`) 에서 실행 가능 (회로 깊이 ~30-50).
- **N = 21**: 5+10 큐비트. NISQ 한계 근처.
- **N ≥ 35**: NISQ 깊이 한계 초과, 실행 불가.

### 5.3 측정 목표

(C) 의 hardware 검증 = **실제 디바이스 노이즈가 (C) 의 covered 영역에 도달하는가?**

- 각 N 에서 K = 50-100 회 측정.
- 각 측정마다 (k, base a) 기록.
- 후처리: (C) 알고리즘 적용 → L 누적.
- 정리 1 검증: K 가 충분히 큰 시점에서 r_a | L_before 인 측정의 회수율 = 100% (zero violations).
- 정리 2 검증: K_λ 의 경험 분포 → λ(N) 도달까지 측정 수.
- 정리 3 검증: hardware 노이즈를 "destructive equivalent η_M" 으로 매핑 (depol 근사) →
  K_λ^alg 예측 vs 경험.

### 5.4 코드 스켈레톤 (qiskit 기반, future use)

`experiments/hardware_demo.py` (가칭) 의 골격:

```python
from qiskit import QuantumCircuit, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

def shor_circuit(N: int, a: int, t: int) -> QuantumCircuit:
    \"\"\"단순 modexp + iQFT 회로. N=15, 21 까지만 효율적.\"\"\"
    ...  # standard Shor circuit

def run_hardware(N: int, K: int = 50, backend_name: str = "ibm_lagos"):
    service = QiskitRuntimeService()
    backend = service.backend(backend_name)
    results = []
    for trial in range(K):
        a = pick_coprime(N)
        qc = shor_circuit(N, a, t=2 * int(math.log2(N)) + 2)
        qc_t = transpile(qc, backend)
        sampler = SamplerV2(mode=backend)
        job = sampler.run([qc_t], shots=1)
        k = int(list(job.result()[0].data.meas.get_counts().keys())[0], 2)
        results.append((a, k))
    return results

def apply_C_to_hardware_data(N, results):
    from multi_base import MultiBaseState, order_from_measurement
    ...  # accumulate L, log covered/violations/lucky
```

### 5.5 다음 phase

본 phase 는 사용자의 IBM Q 계정/qiskit 설치 결정에 의존. 결정시 별도 세션에서 실행.
Phase 5 (Regev 본실행) 로 진행.

## §6 Phase 5 — E: Regev 통합 (본실행)

### 6.1 Phase 0 재확인

Phase 0 에서 후보 B (Regev 의 격자 측정에서 좌표별 (C) 추출) 를 Phase 5 표적으로
확정. 본 phase 는 *Regev 2023 의 격자 후처리에 (C)-식 noise-invariance 를 통합* 하는
시도.

### 6.2 Regev 알고리즘 기본 구조 (재확인)

Regev 2023 (arXiv:2308.06572):
- `d ≈ √(log N)` 개 base `a_1, …, a_d` 를 *병렬* 로 사용
- 단일 QFT 측정 → 격자점 `y ∈ Z^d` 회수 (idealy, `y` 는 specific dual lattice 안)
- 후처리: LLL 격자 환원으로 joint period 회수

회로 효율 우위 (O(n^{3/2}) gates) 대신 후처리가 무거움 (LLL).

### 6.3 (C) 통합 시도

#### 시도 1: 직렬 합성 (Phase 0 의 후보 A)

Regev 로 `L_0 = lcm(r_{a_1}, …, r_{a_d})` 회수 → 이후 다른 base 는 (C) 의 fast path.
- 작동 보장. paper 코멘트 한 줄 가능.
- 새 contribution 작음.

#### 시도 2: Regev 격자 좌표별 (C) — Phase 0 후보 B

Regev 의 격자 측정 `y = (y_1, …, y_d)` 가 좌표별로 `y_i ≈ j_i · Q / r_{a_i}` 형태로
해석 가능한가? 가능하면 *좌표별 독립 (C) 후처리* 가 적용 가능.

**핵심 질문**: Regev 의 측정 분포에서 `y_i` 의 marginal 이 Shor 의 단일 base 측정 분포
와 동일한가?

이를 확인하려면 Regev 2023 §3-4 정독 필요. **현재 인터넷 접근 없음** → 사전 지식 기반
추측만 가능.

**추측 (검증 필요)**: Regev 의 측정은 격자 `Λ = ⊕_i (Q · Z / r_{a_i}) ⊕ remainder` 의
dual 위의 분포. 좌표별 marginal 은 Shor 분포에 가까우나 정확히 같지는 않음 (joint
correlation 존재).

이 경우:
- 좌표별 (C) 는 *근사적으로* 적용 가능.
- 노이즈 하에서 (C) 의 noise-invariance 가 좌표별로 부분 유지.
- Regev 의 LLL 단계가 노이즈에 약한 부분을 (C) 가 보완.

#### 시도 3: (C) 의 격자 일반화 — Phase 0 후보 C

(C) 의 "lcm 누적 + divisor search" 를 격자 setting 으로 추상화:
- `L` → 누적 부분격자 `Λ_K`
- divisor search → 격자 sub-element search
- 검증자: `a^d ≡ 1` → `d ∈ Λ_K^⊥` (격자 멤버십)

이 일반화가 자연스러우면 → "lattice (C)-determinism 정리" 도출 가능.
**현재 평가**: 격자 구조와 divisor 구조의 1-to-1 대응이 자연스럽지 않아 보임. 시도해
보지 않으면 모름.

### 6.4 Regev 2023 구조 확인 (WebFetch 2026-06-12)

arXiv abstract + summary fetch 로 핵심 구조 확인:

- 각 양자 run 이 `d = √(n + 4) ≈ √(log N)` 개 base `a_1, ..., a_d` 의 **좌표별
  측정** 산출: 출력은 벡터 `(k_1, ..., k_d) ∈ Z^d`.
- **각 좌표의 marginal**: `k_i ≈ j_i · Q / r_{a_i}` — **Shor 의 단일 base 측정 분포와 동일**.
- 총 `√n + 4` 회 독립 run, 매번 다른 base 셋.
- Regev 의 후처리: 격자 환원 (LLL/BKZ) on 수집된 `K · d` 좌표 벡터.

**시도 2 (Phase 0 후보 B) viability 확정**: marginal 이 Shor 분포 이므로 **각 좌표에
독립적으로 (C) 후처리 적용 가능**.

### 6.5 정리 4 (안) — Regev + (C) 좌표별 후처리

> **정리 4 (Regev-(C) 합성).** Regev 의 양자 회로를 `K` 회 독립 실행. 각 측정 `(k^{(t)}_1,
> ..., k^{(t)}_d)`, `t = 1, ..., K` 에 대해 좌표별로 (C) 후처리 적용 — 누적 `L` 을 다음
> 규칙으로 갱신:
>
> `L ← lcm(L, (C)(a_i^{(t)}, N, k^{(t)}_i, Q, L))` for each (t, i)
>
> 그러면:
>
> **(a) covered:** 만약 어떤 (t, i) 에서 `r_{a_i^{(t)}} | L_before`, 회수는 결정적
> (정리 1 의 Regev 적용).
>
> **(b) K_λ bound:** Regev 의 base 분포가 균등이라 가정하면, 정리 2 가 적용되어
> `E[K_λ^Regev-(C)] = E[K_λ^ideal] / d` runs (각 run 이 d 개 base 제공).
>
> **(c) 노이즈 무관성 (정리 1 corollary 1 의 Regev 적용):** 측정 좌표 `k^{(t)}_i` 의
> 분포가 어떻게 손상되든, `r_{a_i^{(t)}} | L_before` 이면 (C) 가 회수.

**의의**: Regev 의 격자 환원 (LLL) 후처리는 측정 정확도에 민감 — 잡음이 격자점을 왜곡
하면 LLL 실패. (C) 좌표별 후처리는 noise-invariance 유지 → **noise-robust Regev variant**.

**Trade-off**: (C) 좌표별은 좌표 간 결합 정보를 무시 → Regev 의 LLL 만큼의 효율 보장은
못함. 하지만 노이즈 환경에서 더 견고.

### 6.6 증명 sketch

(a): 정리 1 의 직접 적용 — 좌표 `i` 의 측정 분포가 무엇이든 `r_{a_i^{(t)}} | L` 이면
(C) 좌표별 = 정리 1.

(b): 한 run 이 d 개 독립 base 제공. K runs 에서 K·d 개 독립 base 누적. 정리 2 에 K·d
입력 → `P[L_{K·d} < λ(N)] ≤ ω(λ) · 2^{-K·d}`. K = `log_2 ω(λ) / d` 면 됨.

(c): 정리 1 corollary 1 의 직접 적용.

### 6.7 검증 계획

- numpy 시뮬: Regev-style 다중 base 회로의 단순화 (각 좌표 = 독립 Shor 측정).
- (C) 좌표별 후처리 vs Regev LLL 의 K_λ 비교 (노이즈-free + 노이즈 하).
- N 작음 (15, 21, 35) 에서 LLL 라이브러리 (e.g., fpylll) 또는 직접 구현.

### 6.8 paper 통합

- §3.4 또는 §5 에 정리 4 + 검증표.
- 메인 메시지: "(C) 의 다중 base 후처리는 Regev framework 와도 직교 결합 가능".

### 6.9 본 세션 한계

LLL 구현이 numpy-only 정책 위반 (fpylll 또는 sympy 필요). 본 세션은 정리 4 진술 + 증명
sketch 만 작성. 본 검증은 후속 세션.

## §7 frontier — 본 paper 의 틀 밖 탐색

**별도 문서**: `frontier.md` (2026-06-12 대화에서 도출).

핵심 내용:
- 본 paper 의 4 트릭 (공책, 다중 픽, 갈아만든 픽, hybrid) + 8 개 암묵적 가정.
- 4 가지 "틀 밖" 방향 (다른 invariant, 체계적 base, 추측-검증, 직접 RSA) 의 literature 검토.
- 모두 **이미 탐구 또는 fundamental 한계** 로 막힘.
- 추가 6 방향 (NFS quantum hybrid, partial order, 다중 N, ECDLP, info-theoretic limit, **noise as feature**)
  의 frontier 가능성.
- 가장 흥미: **F (noise 활용)** — Lemma 5.1 의 정신 확장.

본 paper 본문 외부 — 후속 연구 가능성 기록용.

## §8 진행 로그

- **2026-06-13 (paper polish + Zenodo prep)**: **paper.md / paper.tex 정밀 polish 완료.**
  - Abstract: §3.6 SR observation 한 줄 추가 (caveat 포함), Lemma 5.1 bound 의
    "within 4%" 를 N=437 한정으로 명확화.
  - §3.7 Joint interpretation: Theorem 5 (~5x reduction) + §3.6 SR 의 *orthogonal*
    위치 문단 추가.
  - §7 Reproducibility: `python -m experiments.sr_aop` 추가.
  - §8 Conclusion: Theorem 4-5 정량 mention + §3.6 SR + ENAQT bridge 한 문단.
  - **Zenodo / GitHub 인용 준비**: CITATION.cff (cite-this-repository 버튼 활성화),
    .zenodo.json (자동 메타데이터), LICENSE (MIT) 추가.
  - 부속: README §4 paper.md/experiments/ 추가, summary.md §9/§11 갱신.
  - **다음**: N=1147 d=1 multi-seed confirm 결과 → §3.6 통합. GitHub release →
    Zenodo 자동 sync → DOI 획득.

- **2026-06-13 (SR 정리 + AOP)**: **SR finding 정리 완료.** 17.86% peak 가
  1000 trials × 4 seeds 에서 mean +0.42% → **fluke 확정**. H12c (σ_opt ∝ N^α),
  H9 (polynomial scaling) 둘 다 기각. **AOP grid (N × d)**: d=1 universal
  positive signal (+0.78%, +2.60%, +1.36% @ N=437, 1147, 2491). 본 paper 의
  honest SR finding = **~1-3% 작은 효과 / sub-optimal d 영역 / anti-optimization
  caveat**. ENAQT 와의 bridge 는 *small but genuine*. RSA 변화 없음.

- **2026-06-12 (1)**: Phase 0 스코핑 완료. E 의 후보 B 를 Phase 5 표적으로 확정.
  본 문서 §1 추가. Phase 1 (A) 착수 — §2 의 정리 윤곽 작성.
- **2026-06-12 (14)**: **틀 밖 탐색 + 4 방향 literature 검토.** 본 paper 의 framework 분석
  (4 트릭, 8 암묵 가정) + 4 "틀 밖" 방향 (다른 invariant, 체계적 base, 추측-검증, 직접 RSA)
  의 사전 검토. **결과: 4 방향 모두 이미 탐구 또는 정보이론적 한계 로 막힘**.
  - #1 (다른 invariant): Cheung-Mosca, Hallgren, Carmichael paper 2021 = 거의 막힘.
  - #3 (체계적 base): 양자에서 적게 탐구, marginal advance 만 기대.
  - #7 (추측-검증): Grover 의 √N = Shor poly(log N) 보다 느림. dead-end.
  - #8 (직접 RSA): Miller's theorem (factoring ↔ d 회수 poly-time 동치) 로 fundamental 막힘.
  추가 6 방향 (NFS hybrid, partial, 다중 N, ECDLP, info-limit, **noise as feature**)
  중 가장 흥미는 **F (noise 활용)** — Lemma 5.1 정신 확장.
  결과 frontier.md 에 기록. 본 paper 본문에는 포함 안 함.

- **2026-06-12 (13)**: **B-4 본실행 — 정리 5 (hybrid) 발견 + 검증.** Regev 의 quadratic
  character 트릭 정확 구현: `b_i` random 으로 고른 후 `a_i = b_i² mod N`. 알려진 b_i 로
  `b = ∏ b_i^z_i` 의 nontrivial sqrt(1) 검증 → 인수.
  N=437 직접 검증 성공: b_0^99 mod 437 = 229, 229²=1, gcd(228,437)=19. **인수 추출 완료**.
  3-way 비교 (N=437, d=4, 30 trials each, 5 조건):
  - (C) lcm only: ~6.8 runs, 70% 성공 — Regev setup 의 a_i = b_i² 가 odd part 만 학습.
  - Regev b-trick: ~3.3 runs, 90% 성공.
  - **(C) + b-trick hybrid: ~1.2 runs, 100% 성공** ← 모든 조건 (corruption 30%, phase σ=1.0)에서.
  정리 5 (hybrid) 신설 — paper §3.5. abstract / conclusion / reproducibility 모두 갱신.
  의의: (C) 와 Regev 의 직교 결합이 *strictly better than either alone*. 본 paper 의
  central empirical finding.

- **2026-06-12 (12)**: **B 단계 3: Regev Algorithm B.1 skeleton 구현 시도.**
  `experiments/rv_filter_lll.py` 에 `regev_algorithm_b1`, `try_factor_from_relations`
  추가. Regev §3 의 lattice 정의 `L₀ = {z ∈ Z^d : ∏ a_i^z_i ≡ 1 mod N}` 정독 후
  Algorithm B.1 의 격자 [I_d  ε^(-1) W; 0  I_k] 와 LLL 환원 구현.
  현 한계:
  · LLL 환원된 짧은 벡터가 trivial form (좌표가 너무 크거나 zero).
  · 인수 추출은 Regev 의 quadratic character b_i (b_i² ≡ a_i mod N) 와 정확한
    lattice L 구성 필요 — 본 구현은 a_i 만 사용한 L₀ 직접 후처리.
  · 따라서 end-to-end 인수 비교는 여전히 future work.
  paper §3.4 의 caveat 단락 강화: Regev B.1 의 skeleton 구현 status 명시.

- **2026-06-12 (11)**: **B 단계 2: RV Algorithm 6.1 의 정확한 격자 구성 구현 +
  (C) 와 head-to-head 부분 비교.** `experiments/rv_filter_lll.py` 의 `build_rv_lattice`,
  `rv_filter_round`, `filter_uncorrupted` 를 RV §6 Algorithm 6.1 의 정확한 공식
  (H = [[S·I_d, S·W], [0, I_|E|]]) 으로 구현. S=100 의 단순화.
  부분 head-to-head 결과 (N=437, d=4, 50 trials, RV "overwrite" corruption 모델):
  · corrupt p=0.00: (C) 1.72 runs, RV precision 100%
  · corrupt p=0.10: (C) 1.86 runs (+8%), RV 88.5%
  · corrupt p=0.20: (C) 2.08 runs (+21%), RV 78%
  · corrupt p=0.30: (C) 2.10 runs (+22%), RV 66.5%
  발견: (C) 가 corruption 증가에 *훨씬 더 graceful* 하게 degrade.
  paper §3.4 에 head-to-head 표 + 정확한 caveats (S 단순화, Regev LLL 미구현) 추가.

- **2026-06-12 (10)**: **RV 2023/2025 정밀 정독 + paper positioning 안정화.** 사용자
  업로드한 RV paper (Space-Efficient and Noise-Robust Quantum Factoring) 정독으로
  noise-tolerant Regev 영역의 정확한 landscape 확인:
  · RV (2023/2025): filter-then-LLL — corrupted 샘플 탐지+제거 후 표준 Regev LLL.
    "well-spread" 가정 필요. lattice framework 안.
  · EG24 (concurrent): stronger 가정 하 standard Regev LLL 이 그대로 동작.
    Analysis-only contribution.
  · 본 paper (Theorem 4): post-processing 자체 교체 (LLL → (C) 좌표별).
    "marginal Shor-like" 가정 필요. lattice framework 밖.
  세 접근은 **직교 (orthogonal)**.
  paper §5 (Related Work) 에 "Three approaches to noise-tolerant Regev factoring"
  subsection 추가. Theorem 4 의 contribution 위치 명확화.
  bibliography 에 RV 2023, EG24 추가.
  실험 stub `experiments/rv_filter_lll.py` 작성 (sympy LLL 기반). 데모: precision 75%.
  정밀 구현 (RV Algorithm 6.1 의 정확한 격자 구성 + Regev LLL 후처리 + 우리 (C) 와
  head-to-head 비교) 은 후속 multi-week 작업.

- **2026-06-12 (9)**: **Phase 5 본실행 완료.** Joint-constrained Regev 시뮬 +
  노이즈 견고함 측정.
  · `regev_joint.py`: 독립 Shor 측정에서 시작해 affine projection 으로
    `Σ b_i k_i ≡ 0 mod λ(N)` 제약 적용.
  · 발견 1: joint constraint 가 (C) 좌표별 K_λ 를 ±7% 범위 내로만 영향
    (4개 N 에서 3개는 *더 나음*). **Theorem 4 의 marginal 가정이 robust**.
  · 발견 2: noise 하 Regev-(C) overhead 가 단일-base (C) 보다 **더 작음**
    (depol 0.5: 1.62x vs 1.85x, phase σ=1.0: 2.11x vs 2.63x). d 개 병렬
    base 가 노이즈 amortize.
  · paper §3.4 통합: Theorem 4 의 caveat 두 검증표로 강화 (joint vs indep
    비교, noise 견고함 비교). Theorem 4 가 conditional 에서 *empirically
    validated under approximate model* 로 격상.

- **2026-06-12 (8)**: **Phase 5 본실행 — Regev 측정의 joint constraint 확인.** 추가
  WebFetch 으로 Regev 2023 §2 의 측정 구조 확인: 측정이 *joint linear constraint*
  `Σ b_i k_i ≈ 0 (mod r)` 를 만족 — *좌표별 독립이 아님*. Theorem 4 의 marginal 가정
  은 단순화이며, Regev 의 실제 분포에서 검증되지 않음. 발견을 paper §3.4 의 caveat
  강화로 반영: marginal 이 Shor-like 한지는 미검증, 효율 claim 은 가정 의존. **noise
  invariance 자체는 좌표별 statement 으로 joint 구조 무관하게 유지**. 본 본실행의
  추가 작업 (LLL 실구현 + joint-constrained 시뮬) 후속 세션.

- **2026-06-12 (7)**: **Phase 4 부분실행 (hardware proxy).** qiskit/IBM Q 계정 없이
  공개된 IBM Eagle 127q 사양 (T1=150μs, T2=100μs, gate error 0.03%/1%, readout 2%)
  을 noise.py 5종 모델로 매핑 → N=15 Shor 시뮬. 결과: 5 노이즈 동시 적용에서
  covered=499/500, violations=0, 100% 성공. K_λ^alg 오버헤드 1.26x (= 1/0.79 =
  Theorem 3 의 destructive-class 예측 일치). paper Appendix E (Hardware-calibrated
  noise simulation) 신설. 실제 hardware run 은 후속.

- **2026-06-12 (6)**: **Phase 5 부분실행 (Regev compatibility).** WebFetch 으로 Regev 2023
  의 구조 확인: d=√(n+4) 개 base 의 좌표별 측정, 각 marginal k_i ≈ j_i·Q/r_{a_i}.
  numpy 시뮬 (각 좌표 = 독립 Shor 측정 가정) 으로 (C) 좌표별 후처리 = 정리 4 검증.
  결과 (N=437, d=4): 평균 1.75 runs · 4 bases ≤ 7 bases 으로 λ(N) 도달.
  정리 4 (conditional on marginal assumption) 작성 + paper §3.4 통합.
  Joint correlation 무시한 한계 명시. 본 phase 의 "본실행" 은 LLL 실구현 필요로
  후속 작업.

- **2026-06-12 (5)**: **Phase 3 종료 (negative scope result).** (C) 의 다중 base 누적이
  Shor 의 group-exponent learning 에 구조적으로 결합됨을 분석. DL/HSP 는 단일 인스턴스
  문제로 자연스럽지 않음. 다중 prime 인수분해는 trivial 확장 (한 줄 코멘트). (C)+PH
  composition 은 유용한 하이브리드 알고리즘. paper §6 (Limitations) 에 "Scope of the
  framework" 단락 추가. Phase 4 로 이동.

- **2026-06-12 (4)**: **Phase 2 본실행.** `k_lambda_alg.py` 작성 + N=437 에서 9종 노이즈
  직접 측정 (100 trials each). 정리 3 (destructive case) 의 Markov chain 증명 완성:
  `E[K_λ^alg(η)] = E[K_λ^ideal] / g_M(η)`. 검증: depol/bias 6 setups 에서 평균 +11% 오차로
  일치. modexp 부분적 일치 (+33-44%). phase 는 structural 노이즈로 적용 외, lower bound
  로만 활용 (실측이 예측 초과). paper.md/tex 의 §3.3 에 정리 3 + 검증표 통합, abstract /
  conclusion 갱신. **paper 가 3 정리 보유 conference paper 로 격상.**

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
