# SR 발견 → arXiv preprint 진행 계획

본 문서는 SR (phase-noise stochastic resonance for factoring) 발견을
arXiv preprint 까지 가져가는 3-stage 계획.

## Stage 1: 진위 확정 (V) — 1-3 일

**목표**: SR effect 가 fluke 가 아닌 robust 현상임을 확정.

### V1: 다른 seed 재현 (1일)

- N=437 d=4, N=1147 d=2 각각 seed ∈ {1, 2, 3, 4, 5}
- 5 seed 모두 SR ≥ 1% 면 진짜
- 평균 SR + 표준편차로 statistical significance 계산

### V2: 다른 RNG (보류 — V1 + V3 결과 보고 결정)

만약 V1 통과면 V2 생략. fail 시 numpy vs Python random 비교.

### V3: 2000 trials 통계 강화 (1일)

- N=437 d=4, σ ∈ {0.025, 0.05, 0.075, 0.10, 0.20} 각 2000 trials
- SE ≈ 0.034 → 1% 효과가 2 SE 이상이면 확정

### 통과 기준

- V1: 5 seed 모두 mean SR > 1%, 1-sample t-test p < 0.05
- V3: σ=0.025-0.10 의 reduction 가 2 SE 이상 (≥ 0.07)

→ 모두 통과 시 Stage 2 진행.

## Stage 2: 방향 A — N 확장 scaling law — 2-3 주

**목표**: SR % 의 N-의존성 정확화. Scaling law 도출.

### A1: N 시리즈 확장 (1주)

- N ∈ {437, 1147, 2491, 4087, 8009, 16001}
- 각 N 에서 d ∈ {2, 4, 8} × σ ∈ {0, 0.025, 0.05, 0.075, 0.10, 0.20}
- 500 trials per cell
- *N 시간 폭증 주의* — N=16001 의 Q=2^28, 시뮬 memory 한계 가능

### A2: 최적 (d*, σ*) 추적 (1주)

- 각 N 에서 SR 최대 위치 (d*, σ*) 와 SR* 정확 측정
- 회귀: SR*(N) = f(log N) ?
- d*(N), σ*(N) 의 패턴

### A3: 데이터 + figure 정리 (1주)

- SR % vs N 그래프 (primary figure)
- (d, σ) heatmap per N
- Scaling law 식 (가설 + fit)

### Stage 2 deliverable

- 정량적 statement: "SR % grows as O(?(N))"
- 명시적 phase diagram

## Stage 3: arXiv preprint 정리 — 1주

**목표**: SR-focused paper draft 완성, arXiv 제출.

### 구조

```
Title: Stochastic Resonance in Quantum Factoring — Phase Noise as a Resource
       for Hybrid (C) + Regev Algorithm

Authors: (author name TBD)

Abstract: ~200 words. Empirical finding + scale + mechanism + limitations.

§1 Introduction
   - Quantum factoring landscape (Shor, Regev, RV, EG24)
   - Noise as resource (SR, ENAQT, dissipative QC 인용)
   - Our finding: SR observed in hybrid factoring at phase noise

§2 Setup
   - Hybrid algorithm = (C) + Regev b-trick (앞 paper 의 정리 5 인용)
   - Phase noise model

§3 Empirical SR
   - V1, V3 결과
   - σ scan
   - Sign-test significance

§4 Scaling
   - Stage 2 의 N 시리즈 결과
   - Scaling law

§5 Mechanism (preliminary)
   - Phase distributional smearing
   - Sweet spot in d
   - Other noise types null

§6 Related Work
   - Quantum SR (channels)
   - ENAQT
   - Noise-tolerant Shor (Cai 2023, RV 2023, EG24, our paper)
   - *No prior on SR for factoring* (positioning)

§7 Limitations + Open Questions
   - 작은 N 만 검증 (cryptographic scale 미검증)
   - 실 hardware 미검증
   - Mechanism formal proof 없음

§8 Conclusion
```

### 길이 & venue

- ~12 페이지 (workshop 길이) 또는 ~6 페이지 (PRL)
- arXiv 제출시 cs.CR 또는 quant-ph 카테고리

### 부록

- Companion code: github.com/Hashevolution/shor (private → public 시점은 publish 후)

## 시점

| 시점 | 진행 |
|---|---|
| Now ~ +1 day | V 진행 중 |
| +1 day | V 결과 보고 Stage 2 진행 결정 |
| +1-3 weeks | A 진행 |
| +4-5 weeks | arXiv draft 완성 |
| +5-6 weeks | arXiv 제출 |

## 만일 V 실패 시

- effect 가 artifact 면: 본 paper 의 §3.6 만 단순 관찰로 남기기
- arXiv 제출은 정리 1-5 중심으로 진행 (SR 은 minor footnote)
- 새로운 frontier 방향 모색

## 만일 V 통과 + A 강한 결과 시

- SR-focused 별도 paper 제출 (본 paper 에 footnote 형식 인용)
- 또는 본 paper 의 §3.6 을 main result 로 격상
