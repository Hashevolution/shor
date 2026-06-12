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

## §2 Phase 1 — A: 정량적 K_λ 정리 (착수)

### 2.1 목표 진술

반소수 `N = pq` 에 대해, 균등 무작위 base 를 `K` 개 뽑아 `L = lcm(r_{a_1}, …, r_{a_K})`
를 누적할 때:

> **정리 2 (안)**: `P[L < λ(N) | K]  ≤  Σ_{ℓ | λ(N)} ℓ^{-s_ℓ(K)}`
>
> 여기서 `s_ℓ` 는 `(Z/N)*` 에서 ℓ-Sylow 의 cover 빈도. 명시적 형태 derivation 필요.

대략 `K = O(ω(λ(N)) · log(1/ε))` 으로 `P ≥ 1-ε`.

### 2.2 증명 윤곽

1. `(Z/N)* ≅ C_{p-1} × C_{q-1}` (CRT, `gcd(p-1, q-1)` 신경써야 함).
2. 각 prime ℓ | λ(N) 에 대해 `ℓ^{v_ℓ(λ(N))} | L` 이어야 `L = λ(N)`.
3. 균등 무작위 a 가 `ℓ^v` 를 cover 할 확률을 컴포넌트 별 cyclic-group 분석으로 계산.
   - cyclic `C_n` 에서 `v_ℓ(n) = v`. 균등 무작위 x 에 대해
     `P[v_ℓ(ord(x)) ≥ k] = 1 − ℓ^{−(v−k+1)}`.
4. `(Z/N)*` 컴포넌트에서: `v_p = v_ℓ(p-1)`, `v_q = v_ℓ(q-1)`. `v = max(v_p, v_q) = v_ℓ(λ(N))`.
   - `s_ℓ` := `|{x ∈ {p, q} : v_x = v}|` ∈ {1, 2}.
   - 한 base 당 cover 확률 = `1 − (1/ℓ)^{s_ℓ}`.
5. coupon-collector 부등식 + `ω(λ(N)) ≤ ω(p-1) + ω(q-1) ≤ O(log log N)` (Hardy-Ramanujan).

### 2.3 실험 검증 계획

- `N` 30-40 개 (반소수, p≈q size 다양) 에서 100-trial 씩 `L` 누적 진행률 측정
- 경험 `K_λ` 분포와 정리의 상한 비교
- 분포 꼬리가 정리 안으로 들어오는지 확인

### 2.4 다음 세션 시작점

- `roadmap.md` §2.2 의 증명 단계 (3) 의 정확한 진술 검증 (n 이 ℓ-power 가 아닐 때).
- `experiments/k_lambda_dist.py` 작성: N 다양한 반소수에서 K_λ 경험 분포.

---

## §3 진행 로그

- **2026-06-12**: Phase 0 스코핑 완료. E 의 후보 B 를 Phase 5 표적으로 확정.
  본 문서 §1 추가. Phase 1 (A) 착수 — §2 의 정리 윤곽 작성.
