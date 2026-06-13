# 쇼어 알고리즘 연구 (Shor's Algorithm)

양자 인수분해 알고리즘을 numpy 기반 상태벡터 시뮬레이션으로 직접 구현하고 분석한다.

---

## 1. 문제와 의의

### 1.1 정수 인수분해 문제
- 입력: 합성수 N
- 출력: 자명하지 않은 약수 p (1 < p < N)
- 고전 최선 알고리즘인 GNFS의 시간 복잡도: `exp(O((log N)^(1/3) (log log N)^(2/3)))` — 준지수(sub-exponential).
- **RSA 안전성의 근거**가 바로 이 어려움이다.

### 1.2 쇼어의 기여 (1994)
- 양자 컴퓨터에서 N을 시간 복잡도 `O((log N)^3)` 안에 인수분해 가능.
- 다항 시간(polynomial) → 지수적 속도 향상(exponential speedup).
- RSA, DH, ECC 등 현재 공개키 암호계의 안전성을 정면으로 위협 → 포스트양자 암호(PQC) 연구 동기.

---

## 2. 알고리즘 구조

쇼어 알고리즘은 **인수분해 ↔ 주기 찾기(order finding) 환원** 위에 양자 주기 찾기 서브루틴을 얹은 구조다.

```
            ┌───────────────────────────────┐
N (합성수) →│  고전 환원 (Miller, gcd, 등)  │
            └──────────────┬────────────────┘
                           ↓ (a, N)
            ┌───────────────────────────────┐
            │  양자 주기 찾기 서브루틴      │ ← 핵심
            │  (Hadamard + ModExp + QFT⁻¹)  │
            └──────────────┬────────────────┘
                           ↓ r (a^r ≡ 1 mod N)
            ┌───────────────────────────────┐
            │  gcd(a^(r/2) ± 1, N)          │
            └──────────────┬────────────────┘
                           ↓
                       p, q (약수)
```

### 2.1 고전 부분 (환원)

1. N이 짝수면 2를 반환.
2. N = m^k 인 prime power면 m을 반환 (반복 제곱근으로 확인).
3. 1 < a < N에서 a를 무작위로 선택.
4. `g = gcd(a, N)`이 1보다 크면, 운 좋게 약수를 찾음 — 반환.
5. **양자 서브루틴 호출** → a의 N에 대한 곱셈 위수(order) r 획득.
6. r이 홀수면 다시 시도.
7. `a^(r/2) ≡ -1 (mod N)`이면 다시 시도.
8. 그 외에는 `gcd(a^(r/2) - 1, N)`와 `gcd(a^(r/2) + 1, N)`가 자명하지 않은 약수.

> **왜 동작하나?** `a^r ≡ 1 mod N`이므로 `(a^(r/2) - 1)(a^(r/2) + 1) ≡ 0 mod N`. 두 인자가 모두 N의 배수가 아니라면 N과 비자명한 gcd를 갖는다.

### 2.2 양자 부분 (주기 찾기)

목표: `f(x) = a^x mod N`의 최소 주기 r 찾기.

#### 회로 구성
- **계산 레지스터 (counting)**: t큐비트, `t ≈ 2·⌈log₂ N⌉`
- **작업 레지스터 (work)**: n큐비트, `n = ⌈log₂ N⌉`
- Q = 2^t로 둔다.

#### 단계
1. 두 레지스터 모두 |0⟩으로 초기화.
2. 계산 레지스터에 Hadamard:
   `|ψ₁⟩ = (1/√Q) Σₓ |x⟩|0⟩`
3. **모듈러 거듭제곱 (Modular Exponentiation)** `U_a: |x⟩|y⟩ → |x⟩|a^x · y mod N⟩`:
   `|ψ₂⟩ = (1/√Q) Σₓ |x⟩|a^x mod N⟩`
4. 작업 레지스터 측정 — 어떤 값 y₀ = a^(x₀) mod N으로 붕괴.
   계산 레지스터는 `a^x ≡ y₀ (mod N)`인 모든 x들의 균등 중첩이 됨. 그 x들은 `x₀, x₀+r, x₀+2r, …`로 등차수열.
5. **역 QFT** 적용:
   `QFT⁻¹: |x⟩ → (1/√Q) Σₖ e^(-2πixk/Q) |k⟩`
6. 계산 레지스터 측정 → k 얻음.

#### 측정 결과의 의미
측정된 k는 다음을 만족한다 (높은 확률로):
`k/Q ≈ j/r` (j는 정수)

**연분수 전개(continued fractions)** 로 k/Q를 분수로 근사 → 분모로 r 후보를 얻는다.

---

## 3. 핵심 수학

### 3.1 QFT와 주기 검출
함수 f가 주기 r을 가질 때, f의 푸리에 변환은 `1/r`의 정수배에서 강하게 집중된다.
계산 레지스터 상태가 `Σₘ |x₀ + mr⟩`일 때, 역 QFT 후 진폭:

`A(k) = (1/√M) Σₘ e^(-2πi(x₀+mr)k/Q)`

`= (e^(-2πix₀k/Q)/√M) Σₘ e^(-2πimrk/Q)`

`Σₘ e^(-2πimrk/Q)` 항이 `rk/Q ≈ 정수`일 때 보강 간섭 → `k ≈ jQ/r` 부근에서 피크.

### 3.2 연분수
QFT가 정확히 jQ/r을 주지 않더라도, `|k/Q - j/r| ≤ 1/(2Q) ≤ 1/(2r²)` 이면 연분수 전개의 수렴값 중 하나가 정확히 j/r이다 (Q ≥ N²이면 보장).

### 3.3 성공 확률
- 한 번 측정에서 올바른 j와 r을 얻을 확률: `Ω(1/log log r)`.
- 따라서 O(log log N)번 반복하면 높은 확률로 r 회수.
- r이 짝수이고 `a^(r/2) ≢ -1 mod N` 일 확률: `≥ 1/2` (a가 균등 무작위일 때).

---

## 4. 이 저장소의 구성

| 파일 | 역할 |
|------|------|
| `README.md` | 이론 정리 (이 문서) |
| `paper.md` / `paper.tex` | **정리 1-5 + Lemma 5.1 의 정식 논문 (영문)** |
| `compile-notes.md` | LaTeX 컴파일 / arXiv 업로드 가이드 |
| `classical.py` | 고전 환원 + 고전적 위수 계산 (검증/비교용 baseline) |
| `shor.py` | 양자 주기 찾기의 numpy 상태벡터 시뮬레이션 |
| `multi_base.py` | 다중 base λ(N) 누적으로 측정당 회수율 향상 (§7) |
| `noise.py` | 6종 노이즈 모델 + 5종 hybrid b-trick 부패 모델 |
| `demo.py` | 인수분해 데모 + 회수율 비교 실험 |
| `experiments/` | 정리 2-5 + §3.6 multi-boundary mechanism 재현 스크립트 |
| `experiments/sigma_scan_437.py` | §3.6 baseline σ-scan (3 seeds × 12 σ) |
| `experiments/sigma_scan_437_extend.py` | §3.6 extended seeds 4-13 + histogram backfill |
| `experiments/sigma_scan_general.py` | 임의 (N, d) cell 에서 σ-scan + cross-cell verification |
| `experiments/analyze_histograms.py` | §3.6 per-seed K-bin flip 식별 + boundary distribution |
| `section_3_6_draft.md` | §3.6 reframe 초안 ([TBD] 채워진 후 paper 에 반영됨) |
| `summary.md` / `hypotheses.md` / `frontier.md` | 연구 진행 메타 문서 |
| `roadmap.md` | Phase 1-7 격상 로드맵 |

### 4.1 시뮬레이션 한계
- 큐비트 수 ≈ 3·⌈log₂ N⌉. 상태벡터 차원 2^(t+n) → N이 커지면 메모리 폭발.
- 본 코드는 효율을 위해 **두 단계 트릭**을 사용:
  1. 작업 레지스터 측정을 먼저 수행하여 계산 레지스터의 부분상태로 축소.
  2. 계산 레지스터에는 numpy FFT를 그대로 사용 (역 QFT가 곧 정규화된 DFT).
- 이 트릭은 측정 통계상 완전한 양자 회로와 동일한 분포를 산출한다.

---

## 5. 실행

```powershell
python demo.py                         # 기본 데모 (15, 21, 35)
python demo.py 21                      # N=21 인수분해
python demo.py --multi 91              # 다중 base 모드
python demo.py --compare 33 35 77      # 인수분해 측정 횟수 비교
python demo.py --period 33 77 143      # 위수 회수 확률 비교 (§7)
python classical.py 21                 # 고전 baseline
```

---

## 6. 참고 문헌

- Shor, P. W. (1994). *Algorithms for quantum computation: discrete logarithms and factoring*. FOCS.
- Nielsen & Chuang. *Quantum Computation and Quantum Information*, Ch. 5.
- Kitaev, A. (1995). *Quantum measurements and the Abelian Stabilizer Problem*. arXiv:quant-ph/9511026.
- Knill, E. & Mosca, M. — 측정 후처리에서 모든 연분수 수렴값을 시도하는 기법.
- Bach, E. & Shallit, J. *Algorithmic Number Theory*, Vol 1. — Carmichael 함수와 위수 분포.

---

## 7. 위수 회수 확률 향상: 다중 base λ(N) 누적

### 7.1 동기

오일러 정리 `a^φ(N) ≡ 1 mod N` 에서, 더 정확히는 Carmichael 함수에 대해:
> **모든** a ∈ (Z/N)* 의 위수 r_a 는 λ(N) 을 나눈다.

φ(N) 이나 λ(N) 자체를 알면 인수분해와 동치 (φ(N), N → p+q, p·q → p,q) 라서
"먼저 λ(N) 부터 알아내자" 는 순환적이지만, **여러 base 의 위수 lcm**

```
L_k = lcm(r_{a_1}, r_{a_2}, ..., r_{a_k})
```

은 λ(N) 의 약수이며 점근적으로 λ(N) 에 수렴한다 (대개 2~4 개 base 면 같아짐).
일단 L 이 (Z/N)* 의 exponent 가 되면 (`b^L ≡ 1` for random b), 다음 두 가지가 가능:

1. 새 base a 의 위수 r_a = min{ d | L : a^d ≡ 1 mod N } — *측정 없이* 고전적으로 회수
2. L 자체로 Miller-Rabin 식 인수 추출 — `L = 2^t · m`, `a^m → a^(2m) → ... → 1` 시퀀스에서 1 직전이 ±1 아니면 그 값으로 gcd(x±1, N)

### 7.2 알고리즘

각 측정 k 에 대한 후처리:

```
candidates ← convergents(k/Q) ∪ divisors(L_누적)
유효 d 들 ← { d ∈ candidates : a^d ≡ 1 mod N }
r_a ← minimize_order(a, N, min(유효 d 들))
L ← lcm(L, r_a)
```

세부 구현은 `multi_base.py` (`quantum_order_multi`, `MultiBaseState`, `factor_from_exponent`).

### 7.3 실험 결과

`demo.py --period N` 으로 측정. 한 번의 양자 측정에서 r 을 회수하는 확률:

| N | (A) `limit_denominator` | (B) 모든 연분수 수렴값 | (C) (B) + 누적 L |
|---|---:|---:|---:|
| 21 | 43.0% | 43.0% | **100.0%** |
| 33 | 54.5% | 54.5% | **93.5%** |
| 77 | 36.0% | 36.5% | **98.0%** |
| 91 | 35.5% | 36.5% | **98.0%** |
| 143 | 36.5% | 36.5% | **94.5%** |
| 209 | 34.0% | 34.5% | **99.5%** |

(200 trials, base a 균등 무작위)

**해석.**
- (A) vs (B): 거의 차이 없음. `limit_denominator(N-1)` 가 수렴값 중 최선의 분모와 사실상 같음.
- (B) vs (C): 측정당 성공률이 **약 35% → 95% 이상**으로 점프. 한 번의 양자 측정에서 r 을 회수하는 확률이 **2~3배** 증가.
- 누적 L 이 λ(N) 의 약수만 거치므로 false positive 도 없음 — 새 base 의 r 후보가 L 을 나누면 검증이 필연.

`demo.py --compare` 의 인수분해 측정 횟수 비교에서는 차이가 작은데, 이는 작은 N 에서는
이미 단일 측정만으로도 인수분해에 성공할 만큼 운이 좋기 때문 (gcd shortcut 이나 한 번에 r 회수).
N 이 커지고 회로 비용이 측정에 지배될수록 위수 회수 확률 차이가 누적 효율로 직결된다.

### 7.4 한계

- 현재 `divisors(L)` 은 trial division. L ≲ N 범위에서는 충분하나 더 큰 N 에서는 소인수 분해 필요.
- (C) 의 초기 트라이얼들 (L 이 아직 작을 때) 은 (B) 와 동등. λ(N) 회수까지의 transient 가 있음.
- 양자 회로 자체의 노이즈/디코히어런스 가정 없음 (이상 시뮬레이션). 실제 하드웨어에서는 측정 분포가 흐려져 (A)/(B)/(C) 모두 성능이 떨어지는데, 그 환경에서 (C) 의 견고함은 향후 실험 대상.

### 7.5 선행연구와의 관계 (정직한 attribution)

위 알고리즘의 구성요소를 학계 선행연구에 맞춰 솔직히 정리한다.

| 우리 구현 | 선행 |
|---|---|
| 같은 base 다중 측정 → 분모 lcm (`quantum_order` 의 일부) | **Knill (1995)** Los Alamos tech report. 표준 기법. |
| 1회 측정에서 모든 연분수 수렴값 분모 시도 | **McAnally (2001)** [arXiv:quant-ph/0112055](https://arxiv.org/abs/quant-ph/0112055): "all denominators of convergents up to the first ≥ n" + Fourier 모듈러스 Q ≈ 2wn³ 확대로 단일 run 확신도 →1. |
| `r | λ(N)` 와 lcm of orders → λ(N) | **Bach-Shallit** *Algorithmic Number Theory*. 교과서. |
| 단일 measurement 위수의 sharp bound (~94%) | **Bourdon-Williams (2007)** [arXiv:quant-ph/0607148](https://arxiv.org/abs/quant-ph/0607148). 이론적 lower bound (알고리즘 변경 X). |
| **다중 base** 측정·order들의 lcm으로 λ(N) 직접 계산 | **2021 Carmichael paper** [arXiv:2111.02488](https://arxiv.org/abs/2111.02488). Algorithm 1: `for k in 1..K: a←random; λ←lcm(λ, ord(a))`. K = O((log N)²). |
| 단일 base 위수만으로 N 완전 인수분해 | **Ekerå (2021)** [arXiv:2007.10044](https://arxiv.org/abs/2007.10044). 다른 방향의 단일 측정 노선. |
| 위수 회수 확률 종합 분석 | **Ekerå (2024)** [arXiv:2201.07791](https://arxiv.org/abs/2201.07791). Survey: 연분수, 격자, 오프셋 탐색, 약수 탐색, Seifert joint solving. **다중 base 누적은 survey에 부재**. |

**우리 구현의 차별점.** 우리가 한 일은 위 요소들의 *특정 조합*:

1. **다중 base 누적 L을 새 측정의 CF 후처리 후보풀에 합치기** — `candidates ← convergents(k/Q) ∪ divisors(L)`. McAnally 의 "all convergents" 후처리와 Carmichael paper 의 "multi-base lcm" 을 묶음. 명시적 문헌은 찾지 못함 (자명한 조합이라 folklore 일 가능성).

2. **L 이 새 base 의 exponent가 되면 측정 생략** — `if pow(a, L, N) == 1: return min divisor d with a^d ≡ 1`. 표준 reduction (위수 = exponent 의 최소 양호한 약수) 의 알고리즘적 활용. 새로운 수학은 아님.

3. **N ≲ 209 작은 N에서의 정량 검증** — 측정당 회수율 ~35% → ~95% 점프 (200 trials × 6 N). 같은 metric으로 측정한 선행 표를 찾지 못했으나 시뮬레이션 규모가 작아 학술 기여로 보기는 어려움.

**솔직한 평가.** 알고리즘적으로 새로운 것은 거의 없다. 가치는:
- 여러 분산된 기법을 하나의 일관된 시뮬레이션 환경에서 결합·비교
- 측정당 회수율 점프를 시각적으로 보여주는 정량 데이터
- 교육·연구 entry point 로서 numpy 200줄 분량의 self-contained 구현

학술 기여를 노린다면 §7.4 의 한계 — 특히 **노이즈 모델 하 견고함** 또는 **베이지안 후처리와의 결합** — 으로 가야 한다.

### 7.6 Ekerå 단일-측정 + smoothness extension 과의 비교 실험

`multi_base.shor_quantum_ekera` 는 Ekerå 2021 의 핵심 아이디어를 구현:
> 한 base 의 위수 r 에 작은 소수들의 거듭제곱을 곱해 r′ = r · ∏(q^⌊log_q(m')⌋) 를 만든 뒤,
> r′ 을 exponent 로 가정해 `factor_from_exponent` (Miller-Rabin 식) 실행.

다중 base 누적과의 비교 (`demo.py --compare3`, shots=1/base):

| N | 원본 single-base | multi-base | Ekerå smoothness |
|---|---|---|---|
| 33 | 100% / 1.98 meas | 100% / 1.04 meas | 100% / 1.04 meas |
| 77 | 100% / 2.68 meas | 100% / 1.70 meas | 98% / 1.34 meas |
| 143 | 100% / 2.24 meas | 100% / 1.94 meas | 94% / 1.36 meas |
| 209 | 100% / 3.68 meas | 98% / 2.34 meas | 90% / 1.58 meas |
| 247 | 100% / 2.74 meas | 100% / 2.74 meas | 88% / 1.48 meas |
| 299 | 100% / 2.74 meas | 100% / 1.84 meas | 98% / 1.60 meas |
| 323 | 100% / 2.56 meas | 98% / 2.14 meas | 96% / 1.86 meas |
| 377 | 100% / 2.98 meas | 98% / 2.28 meas | 94% / 1.76 meas |

**해석.**
- **Ekerå** 가 가장 적은 측정으로 끝남 (1.3~1.8 meas/trial) 이지만 88~98% 성공률. 단일 측정 + smoothness 가 N 이 커질수록 부족.
- **Multi-base** 는 측정을 조금 더 쓰고 (1.7~2.7 meas) 98~100% 성공률. λ(N) 회수까지 가는 transient 덕분에 안정적.
- **원본 single-base** 는 내부 shots=8 까지 보강해 가장 견고하지만 항상 비싼 budget.

이는 두 접근이 **상호 배타적이지 않다**는 것을 시사: r 회수 후 `ekera_extend(L)` 으로 한 번 더 확장하면 multi-base 의 견고함 + Ekerå 의 압축성을 결합 가능. 향후 실험 후보.

### 7.7 노이즈 견고함 실험

`noise.simulate_period_finding_noisy` 는 **depolarizing 채널** 을 모델:
> 확률 p 로 측정 결과 k 가 균등 무작위 (0..Q−1), 그렇지 않으면 표준 측정.

p 는 거시적 노이즈 강도. 실제 NISQ 디바이스의 게이트 오류·디코히어런스를 한 점으로 압축한 거친 근사.

`demo.py --noise N` 으로 측정한 (A)/(B)/(C) 회수율 (300 trials 마다):

**N = 33** (λ=10)

| p | (A) limit_denom | (B) 수렴값 | (C) +누적 L |
|---:|---:|---:|---:|
| 0.0 | 53.0% | 53.0% | 95.7% |
| 0.2 | 45.7% | 49.7% | 98.3% |
| 0.4 | 40.0% | 44.0% | 98.3% |
| 0.6 | 29.3% | 40.7% | 98.3% |
| 0.8 | 22.3% | 37.7% | 96.0% |

**N = 77** (λ=30)

| p | (A) | (B) | (C) |
|---:|---:|---:|---:|
| 0.0 | 37.0% | 37.3% | 98.7% |
| 0.4 | 37.0% | 40.3% | 100.0% |
| 0.6 | 22.0% | 28.0% | 100.0% |
| 0.8 | 16.3% | 27.0% | 99.0% |

**N = 143** (λ=60)

| p | (A) | (B) | (C) |
|---:|---:|---:|---:|
| 0.0 | 36.3% | 36.3% | 96.3% |
| 0.4 | 31.7% | 35.3% | 98.3% |
| 0.6 | 21.0% | 25.3% | 96.7% |
| 0.8 | 16.0% | 22.0% | 97.3% |

**핵심 관찰.**
1. (A), (B) 는 p 가 커지면서 단조 감소 (예상대로 노이즈에 직접 노출).
2. **(C) 는 노이즈 강도에 거의 무관하게 96~100% 유지**. p=0.8 인 거의 무작위 측정에서도 견고.
3. 이유: 누적 L 의 약수 집합 자체가 측정과 무관하게 r_a 후보를 포함. 후처리에서 `a^d ≡ 1` 검증이 r_a 를 골라내기에 충분.

**해석과 한계.**
- (C) 는 사실상 **측정이 거의 의미 없어진 상태에서 고전적 약수 검색** 으로 환원. L 이 일단 의미 있게 누적되면 (위 표에서 L_final ∈ {10, 30, 60} = λ(N)) 추가 측정은 보강용.
- 단, **초기 transient** — 첫 1~2 trial 에서 L=1 이라 (C) ≡ (B). 노이즈가 크면 L 누적 자체가 늦어져 transient 가 길어짐.
- depolarizing 외 노이즈 (위상 디코히어런스로 인한 peak 분산, gate error 로 인한 f(x) 잘못 계산) 는 별도 실험 필요. (C) 의 견고함이 모든 노이즈 모드에서 유지된다는 보장은 없음.

**향후 연구 가치.** 실제 NISQ 하드웨어의 depolarizing 모델이 유의미한 부분 (p ≲ 0.3) 에서 (C) 가 ~98% 회수율을 유지한다는 것은, 다중 base 누적 후처리가 단순한 수학적 최적화 이상의 **실용적 노이즈 완화 기법** 임을 시사. 이 방향이 §7.5 에서 언급한 "학술 기여" 의 잠재 노선.

### 7.8 가설 H4 검증 — (C) 견고함의 보편성

**원본 가설 H4** (§7.4 에서 던진 한계):
> "위상 디코히어런스 (peak 분산) 하에서 (C) 의 견고함은 depolarizing 만큼 견고하지 않을 수 있다."

검증 위해 세 노이즈 모델 구현 (`noise.simulate_period_finding_noisy`):

| 모델 | 영향 지점 | 파라미터 |
|---|---|---|
| **Depolarizing** | 측정 k 자체 | p ∈ [0, 1], 균등 무작위 교체 확률 |
| **Phase decoherence** | iQFT 직전 amp | σ (라디안), 가우시안 위상 노이즈 |
| **ModExp error** | f(x) = a^x mod N 자체 | q ∈ [0, 1], 잘못 계산 확률 |

`demo.py --noise3 N` 으로 측정. 초기 5 trials (transient) 와 나머지 (steady-state) 분리.

**N = 77, steady-state 회수율:**

| 모델 | level | (A) | (B) | **(C)** |
|---|---:|---:|---:|---:|
| depol | 0.5 | 27.1% | 30.2% | **100.0%** |
| depol | 0.8 | 16.3% | 26.8% | **99.7%** |
| phase | 0.5 | 34.2% | 37.3% | **100.0%** |
| phase | 0.8 | 20.7% | 26.1% | **100.0%** |
| modexp | 0.5 | 13.2% | 23.7% | **99.7%** |
| modexp | 0.8 | 8.5% | 16.3% | **94.2%** |

**N = 143:**

| 모델 | level | (A) | (B) | **(C)** |
|---|---:|---:|---:|---:|
| depol | 0.8 | 15.9% | 21.4% | **98.0%** |
| phase | 0.8 | 17.6% | 22.0% | **98.6%** |
| modexp | 0.8 | 5.1% | 13.9% | **98.3%** |

**결론: H4 는 false.** (C) 는 *측정-층 노이즈 (depol, phase) 뿐 아니라 구조-층 노이즈 (modexp) 에서도 95%+ 견고*.

#### 왜? — 단순 정리

세 모델 모두에서 (C) 가 견고한 이유를 한 줄로:

> **(C) 의 회수 성공은 측정 분포 P(k) 와 무관하며, 누적 L 의 다음 조건에만 의존한다: `r_a | L`.**

**증명 (간단):**
- 후보 풀: `candidates = convergents(k/Q) ∪ divisors(L)`.
- `r_a | L` 가정. 그러면 `r_a ∈ divisors(L) ⊆ candidates`.
- `a^{r_a} ≡ 1 mod N` 이므로 `r_a` 는 검증 통과 (`pow(a, r_a, N) == 1`).
- 어떤 `d ∈ candidates` 가 `pow(a, d, N) == 1` 만족하면 `d` 는 `r_a` 의 배수.
- 따라서 `min{d ∈ candidates : a^d ≡ 1}` 가 `r_a` 의 양의 배수 중 최소 → `r_a` 자체.
- `minimize_order` 는 이미 최소면 그대로 반환. ∎

즉 측정 결과 k 가 어떻든 — 균등 무작위든, blurred peak 든, modexp 가 90% 부패하든 — `r_a | L` 만 성립하면 (C) 는 *결정적으로* `r_a` 를 회수한다.

#### 자명한 한계: 'r_a | L' 조건

위 정리는 *조건부* 결정성. 조건 `r_a | L` 의 확률을 분리해보면:
- 첫 측정 시 L = 1, `r_a | 1` 성립할 수가 없음 (r_a ≥ 1 이지만 r_a = 1 은 a = 1 뿐) → (C) ≡ (B).
- L 이 자라면서 `r_a | L` 의 확률 상승.
- L = λ(N) 이면 모든 a 에 대해 `r_a | L` 보장 → (C) = 100%.

표의 transient 구간 (초기 5 trials) 에서 (C) 가 (B) 와 비슷한 이유가 이것.

#### 학술적 시사

이 정리는 다음을 의미:
1. **(C) 의 노이즈 견고함은 측정-층 노이즈 모델에 *불변*.** 어떤 새로운 노이즈 모델을 들고 와도 *같은 보편성* 이 성립. — H4 가 false 임의 깊은 이유.
2. **연구 가치 있는 질문은 노이즈 견고함이 아닌 "L 누적 속도".** 즉 §1 에서 던진 **가설 H2** ("무작위 semiprime N=pq 에서 다중 base 누적 L 이 λ(N) 에 도달하는 base 수의 평균 분포") 가 이제 *진짜* 핵심 미해결 질문으로 남음.
3. modexp 노이즈가 일반 노이즈 모델보다 (C) 를 더 빨리 떨어뜨리지도 *못함* — 이는 (C) 의 견고함 메커니즘이 본질적으로 *quantum-classical 경계의 비대칭성* (양자 측정은 노이즈, 고전 검증은 결정적) 에 기댄 결과임을 시사.

### 7.9 H2 학계 검증 결과 — 무게중심 재이동

H2 ("E[K_λ] = O(log log N) for semiprime") 의 학계 선행연구 종합 검토 결과:

- **Carmichael paper (2021, arXiv 2111.02488)** 가 K = O((log N)² · log d) 를 증명, semiprime (d=2) 경우 O(log log N) 으로 환원. → **H2 사실상 이미 증명됨**.
- **Pomerance et al. (2017)** e(G) ≤ d + 2.752 ⇒ semiprime (Z/N)\* 에서 평균 상수 회 만에 생성 가능.
- **Erdős-Pomerance-Schmutz (1991)** 및 후속 (Kurlberg-Pomerance, Pollack) 가 λ(N) 분포 / 평균 위수의 정밀 분석 완료.

→ H2 는 결과를 numpy 로 재확인하는 수준이 되어 학술 가치 낮음. **무게중심을 §7.8 의 (C)-determinism 정리** 로 재이동.

---

## 8. (C)-determinism 정리: 정식화와 보편성 검증

§7.8 에서 자연스럽게 유도된 정리를 **explicit 형태로 정식화**하고, *어떤 측정-층 노이즈 하에서도 위배되지 않음* 을 직접 검증한다.

### 8.1 정리 statement

**가정.**
- N: 양의 합성수 (소수 거듭제곱 아님).
- a ∈ (Z/N)\*, r_a := ord_N(a).
- L ∈ ℕ: 알고리즘이 누적한 exponent 후보 (이전 base 들의 위수 lcm).
- k ∈ {0, 1, ..., Q−1}: 측정값 (Q = 2^t, t = 계산 레지스터 큐비트 수).
- **후보 풀**: candidates(k, L) := convergents(k/Q) ∪ divisors(L).
- **후처리 출력 (C)**: `(C)(k, L) := minimize_order(a, N, min{d ∈ candidates : pow(a,d,N) == 1})`, 유효 후보 없으면 0.

**정리 (C-determinism).** 만약 `r_a | L`, 그러면 `(C)(k, L) = r_a` (k 와 무관하게 결정적).

### 8.2 증명

가정: `r_a | L`. 따라서 `r_a ∈ divisors(L) ⊆ candidates(k, L)`.
`pow(a, r_a, N) = 1` 이므로 `r_a` 는 유효 후보 (검증 통과).

유효 후보 집합을 `V := {d ∈ candidates : pow(a, d, N) = 1}` 라 두자.
`r_a ∈ V` 이므로 `V ≠ ∅`.

위수의 정의: `r_a` 는 `pow(a, m, N) = 1` 을 만족하는 *최소* 양의 정수.
따라서 임의 `d ∈ V` 에 대해 `r_a ≤ d`. 그리고 `r_a ∈ V` 이므로 `min(V) = r_a`.

`minimize_order` 는 입력이 이미 최소 (`pow(a, r_a/p, N) ≠ 1` for all prime p | r_a, by definition of order) 이므로 `r_a` 그대로 반환. ∎

### 8.3 따름정리

**Corollary 1 (노이즈 불변성).** 위 정리는 측정 분포 `P(k)` 와 무관하다. 따라서 **어떤 측정-층 노이즈 모델 하에서도**, `r_a | L` 만 성립하면 (C) 는 결정적으로 `r_a` 를 회수한다.

**Corollary 2 (L 의 단조 무결성).** 알고리즘이 (C) 의 출력으로 L 을 업데이트할 때, 매 업데이트는 `r_a` (정확한 위수) 로 이뤄진다. 따라서 알고리즘 진행 중 누적 L 은 λ(N) 의 약수만을 lcm 하여 *결코 corrupted 되지 않는다* — 이는 노이즈 하에서도 보장.

**Corollary 3 (실패 영역).** (C) 의 실패는 `r_a ∤ L` 인 경우에만 발생하며, 이 경우 성공은 convergents(k/Q) 에 r_a 의 양의 배수가 있는지에 의존 (일반 (B)-경로).

### 8.4 직접 검증 — 11 종 노이즈 × 3 N × 500 trials

`demo.py --verify N` 으로 정리를 직접 시험. 각 측정마다 (r_a, L_before, condition, (C) success) 로깅.

**측정 메트릭:**
- `covered`: `r_a | L_before` 인 trial 수 (정리의 가정 영역)
- `violations`: covered 중 (C) 가 r_a 회수 실패한 수 (**0 이어야 정리 성립**)
- `lucky`: ¬covered 중 (C) 성공 (convergent 경로)
- `missed`: ¬covered 중 (C) 실패

**결과 (N=77 / 143 / 209, trials=500):**

| 노이즈 | N=77 violations | N=143 violations | N=209 violations |
|---|:-:|:-:|:-:|
| noise-free | 0 | 0 | 0 |
| depolarizing p=0.3 | 0 | 0 | 0 |
| depolarizing p=0.8 | 0 | 0 | 0 |
| readout flip p=0.3 | 0 | 0 | 0 |
| bias_zero p=0.5 | 0 | 0 | 0 |
| phase σ=1.0 | 0 | 0 | 0 |
| phase σ=2.5 | 0 | 0 | 0 |
| amp_damp γ=0.01 | 0 | 0 | 0 |
| amp_damp γ=0.05 | 0 | 0 | 0 |
| modexp q=0.3 | 0 | 0 | 0 |
| modexp q=0.8 | 0 | 0 | 0 |

**총 16,500 측정 중 violations = 0**. 정리는 *모든* 시험 노이즈 모델·강도에서 성립.

### 8.5 부속 통계 (success_C 와 missed)

전체 (C) 성공률은 covered 의 100% 와 lucky 의 일부 합으로 결정. 가장 가혹한 케이스:

| N | 노이즈 | covered/500 | lucky | missed | success_C |
|---|---|---:|---:|---:|---:|
| 77 | modexp q=0.8 | 476 | 3 | 21 | 95.8% |
| 143 | bias_zero p=0.5 | 480 | 4 | 16 | 96.8% |
| 209 | phase σ=2.5 | 437 | 3 | 60 | 88.0% |

이 미스율은 모두 **초기 trial 에서 L 이 미축적된 상태 (¬covered) + 노이즈 강함 → convergents 도 실패** 의 조합. 일단 L 이 λ(N) 에 가까워지면 covered 의 100% 가 보장됨.

### 8.6 학술적 함의

1. **(C) 의 견고함은 우연이 아니라 정리.** §7.8 의 경험적 관찰이 4줄 증명으로 확정.
2. **노이즈 모델 일반화 가능.** 어떤 새 측정-층 노이즈 (정신적 모델, 하드웨어 특이 채널) 에서도 자동으로 적용.
3. **알고리즘 설계 함의.** L 누적이 진행될수록 (C) 가 *결정적으로* 노이즈 면역이 됨 → adaptive shot allocation 의 이론적 토대.
4. **선행연구와의 차별점.** 위 정리 statement 의 *explicit 형태* (수렴값 ∪ 약수 + 검증 게이트의 결정적 분리) 는 검색한 문헌 (McAnally 2001, Carmichael 2021, Ekerå 2024 survey) 에서 찾지 못함. *folklore optimization 의 정식화* 로 워크숍/short paper 후보.

### 8.7 한계와 후속 질문

1. **양자 회로 자체가 깨졌을 때.** 모듈러 거듭제곱이 *완전히* 무작위면 (modexp_error → 1) L 이 누적 안 되므로 covered 영역도 못 만들어짐. 이 경우 (C) 는 무력 (정리 적용 안 됨).
2. **고전 검증의 노이즈.** 정리는 `pow(a, d, N)` 이 정확하다는 가정. 실제 하드웨어에서는 고전 검증 자체에 hardware 노이즈가 있을 수 있으나, 통상 양자 노이즈 대비 무시 가능.
3. **adaptive base 선택 (H5).** 정리가 보장하는 것은 covered → success. covered 의 비율을 어떻게 최대화하는가? base 선택 전략 설계 여지.
4. **격자 후처리 (Knill-Mosca) 와의 결합.** convergents 경로를 격자 enum 으로 강화하면 ¬covered 의 lucky 비율 증가 가능.

→ §9 에서 위 후속 후보 중 H5 (adaptive base selection) 와 격자 후처리의 결과 수록.

## 9. 후속 실험: Adaptive base + 격자 후처리 (negative results)

§8 정리의 함의를 알고리즘적으로 끝까지 밀어붙여, 측정 횟수를 더 줄일 수 있는지 점검.

### 9.1 H5 — Adaptive base 선택

**아이디어.** 정리에 의해 `r_a | L` 인 base 는 *고전적으로* (측정 없이) 회수 가능. 따라서 양자 측정 가치는 `r_a ∤ L` ⇔ `a^L ≢ 1 mod N` 인 base 에서만 발생. 사전 필터링으로 정보-greedy base 선택.

**구현 (`multi_base.shor_quantum_adaptive`).**
```python
def adaptive_base_select(N, L, rng):
    for _ in range(max_tries):
        a = rng.randrange(2, N)
        if pow(a, L, N) != 1:  # L 을 확장할 가능성
            return a, "extending"
    return None, "saturated"  # L ≈ λ(N) 추정
```

**비교 (`demo.py --adaptive`, shots=2/base, 100 trials):**

| N | random base (meas / iter) | adaptive (meas / iter) | 차이 |
|---|---|---|---|
| 77 | 2.01 / 1.60 | 2.04 / 1.60 | 0 |
| 143 | 2.11 / 1.62 | 2.11 / 1.62 | 0 |
| 209 | 2.78 / 1.99 | 2.78 / 1.98 | 0 |
| 323 | 2.43 / 1.68 | 2.46 / 1.68 | 0 |
| 377 | 2.16 / 1.55 | 2.20 / 1.56 | 0 |
| 437 | 2.53 / 1.91 | 2.61 / 1.89 | 0 |

**결과: 차이 없음.** 측정 횟수도 iteration 횟수도 동일.

**원인 (정리의 corollary).** `quantum_order_multi` 의 **fast-path** 가 이미 정리 적용 영역을 자동 검출하고 측정을 건너뛴다 (`is_exponent_for(a, N)` 체크). Adaptive 사전 필터링은 같은 영역을 다른 방식으로 식별할 뿐, 추가 절약할 측정이 없음.

**Corollary 4 (Adaptive 무용성).** §8 의 (C)-determinism 정리와 fast-path 메커니즘의 조합 하에서, 무작위 base + fast-path 는 측정 횟수 면에서 *어떤 사전 필터링 base 선택 전략과도 동등*. 측정 절약을 위해 다른 메커니즘 (격자 후처리, 측정 budget 적응 등) 이 필요.

### 9.2 격자 후처리 (Knill-Mosca식) — 분석적 평가

**아이디어.** 단일 측정에서 `Fraction(k, Q).limit_denominator(N-1)` 대신, multiple 측정의 k₁/Q, k₂/Q 를 동시 격자 분석하면 r 회수 안정성 증가 가능.

**현 구현의 위치.** `multi_base.convergent_denominators` 는 이미 단일 측정의 모든 연분수 수렴값 분모를 열거 (Knill-Mosca 의 single-measurement post-processing). 다중 측정에서는 `quantum_order_multi` 의 shots 루프가 분모들의 lcm 을 누적 (Mosca 의 multi-measurement extension).

**격자 추가 효과.** Joint diophantine approximation 으로 r 회수의 *failure mode* 를 분석:
- `gcd(j, r) > 1` 인 측정 → 수렴값 분모는 r 의 진약수
- 다중 측정 lcm 은 이를 보완 (서로 다른 j 들이 서로 다른 약수 정보 제공)
- 격자 분석은 lcm 보다 marginal 한 추가 정보만 제공 (개별 측정의 수렴값 분모가 이미 모든 가능한 r 후보를 망라하기 때문)

**결론: 격자 후처리는 우리 framework 에 marginal.** Knill-Mosca 의 가치는 single-base, single-measurement 시나리오에서 가장 큼. 우리 multi-base + (C)-determinism framework 에서는 이미 lcm 으로 대부분 capture.

### 9.3 §9 종합 — Negative results 의 의미

두 후속 실험 모두 측정 횟수 감소 면에서 효과 없음. 이는:
1. **§8 정리의 framework 가 이미 측정-효율 한계에 근접** 함을 시사 (lower bound 후보).
2. 추가 개선은 다른 차원 — 더 큰 N 으로의 scaling, 노이즈 적응, 격자-Regev 노선 등 — 에서 모색해야.
3. 본 codebase 의 학술적 기여는 §8 의 정리 + 그 verification + 자명한 한계 명시에 머무름. 정직하게 평가.
