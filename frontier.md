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

진행 로그:
- 2026-06-12: 초안 — §1-§8 작성.
