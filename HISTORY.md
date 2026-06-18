# PROJECT JAMES 연구 연대기 (v0.1.0 → 현재)

*최초 쇼어 후처리 연구부터 현재 magic(비안정자성) 자원이론까지 시간순 정리. 기준일 2026-06-18.*
*릴리스별 상세는 `release_notes_v*.md`, 명제는 `magic-results.md`, 통합 로드맵은 `magic-program-overview.md`.*

---

## 큰 줄기 (한눈에)

**쇼어 후처리의 "잡음=자원"(SR) 가설** → *정직한 자기수정*으로 닫힌형·불가능정리화(SR은 막다른 길로
확정) → **"진짜 자원은 magic(비안정자성)"으로 방향 전환** → **부호이론 융합**으로 심화.
일관된 메타-특징: **claim → 검증 → 자기수정(retraction)** 의 반복(누적 9+건 철회).

| 단계 | 핵심 질문 | 결론 |
|---|---|---|
| v0.1–0.2.1 | 잡음이 쇼어 후처리의 *자원*인가? | 메커니즘은 보편, **방향은 확률적**(우위 아님) |
| v0.3–0.4 | 그 SR 효과의 정체는? | **표준 dephasing 닫힌형**, Theorem 6로 **가속 불가 확정** |
| v0.5+ | 그럼 우위의 *진짜* 자원은? | **magic의 양/밀도** — 다항(유한) vs 지수(최대)를 가름 |

---

## ① v0.1.0 — 최초 릴리스: 쇼어 차수찾기 다중베이스 후처리 이론
*"A Noise-Invariant Determinism Theorem for Multi-Base Post-Processing in Shor's Order Finding"*
- **정리 1–5 + 보조정리 5.1**: T1 잡음 불변 결정성, T2 로그 커버리지 시간, T3 정확 잡음 스케일링,
  T4 조건부 Regev 호환성, T5 하이브리드 (C)+Regev b-trick 인수분해, 보조정리 5.1(베이스별 비자명
  제곱근 확률 ≥ 1/2).
- 6개 합성수 크기 **17,700 측정**으로 T1 검증, 하드웨어 보정 잡음 시뮬(Appendix E).
- 순수 numpy, 양자 라이브러리 불요.

## ② v0.2.0 (2026-06-13) — §3.6 다중경계 메커니즘 (SR 관찰) · DOI zenodo.20679807
- **13 seeds × 200 trials × 12 σ = 31,200 측정** @ (437,4): **13/13 boundary-flip** 보편 관찰.
- (1147,2) **high-K rescue**(K=15→5 등), (C) 증강=잡음 버퍼; 고전 SR(Benzi–Buchleitner) 형태 대응.
- 알고리즘-구조 regime map(예측). ★ **과대주장 6건 철회**(17.86% peak, SR∝N^α, σ_opt∝N^α,
  AOP, V3 p=0.03, Goldilocks).

## ③ v0.2.1 (2026-06-14) — regime map 5/5 측정 + 보편적 방향 확률성 · DOI zenodo.20681847
- 다섯 항목 전부 직접 측정(Pure Shor, Pure Regev/충실한 LLL, cross-cells).
- **핵심**: 방향은 알고리즘 무관하게 *base-set 확률적*; 차이는 방향이 아니라 *seed별 크기 분포*.

## ④ v0.3.0 (2026-06-14) — 닫힌형 σ-곡선 (errata) · DOI zenodo.20685015
- §3.6 boundary-flip을 **단일 닫힌형**으로 대체: `p(σ)=ρ+(p₀−ρ)·exp(−σ²)`,
  `E[K(σ)]=(1−(1−p)^M)/p` (FFT 위상잡음 dephasing에서 직접 유도).
- **5개 알고리즘 교차검증**(Grover/Shor/QPE/Simon/Hybrid, R²∈[0.88,0.99]).
- Yang–Markidis(arXiv:2605.16074) 보완 위치; §3.6 내부 **3건 철회**.

## ⑤ v0.4.0 (2026-06-14~17) — 정리 6: SR 불가능정리 · DOI zenodo.20688069
- **Theorem 6 (부호 불문 no-go)**: coherence-loss류 잡음에서 ① 내부 공명 없음(최적=끝점),
  ② 닫힌형 swing `|ΔK|=E[K]·|1/g_∞−1/g₀|`, ③ 점근 우위 없음(O(log log N) 상수배).
- "잡음은 항상 해롭다"(g_∞≤g₀) 가정을 N-scaling 실측이 반증 → 부호 불문 재서술.
- Yang–Markidis 포지셔닝을 "수치검증+경계지도"로 정직하게 하향; `paper.md` 단일 정본 선언.
- → **SR 방향 "정직하게 닫힘"**(가속 불가 확정).

## ⑥ v0.5.0 (2026-06-17) — 방향 전환: magic(비안정자성) · DOI zenodo.20725965
*기존 order-finding 정리(T1–T6)는 유지, 새 방향 추가.*
- **`magic.py`**: SRE $M_2$ via XOR-FWHT, 3중 교차검증 10⁻¹⁵, **회귀 42 assert**.
- **명제**: 보조정리 1(평탄상태 $M_2{=}0\iff$ 아핀), 따름정리 1(오라클-가림: Simon 0 / Shor >0),
  따름정리 2(comb), 명제 2·3(Grover **정점 정확히 3 bit**), 명제 2′(일반 M).
- **속도우위 사다리**: Simon(0) → Grover(3bit, 밀도→0) → Shor(∝t→최대).
- 주 선점 2605.05347(EPFL), 양자걷기 2편 DOI 확인.

## ⑦ v0.5.0+ (2026-06-18, 현재 세션) — 부호이론 융합 + 선행조사 + 통합본
- **인수인계서 검토**: "최소 해밍거리" 지표 부정확 규명 → 정확 객체 = **자기상관 $A_W$/Walsh 4차모멘트**.
- **`marker_code_magic.py`**: 정확형 항등식(4e-15), 아핀⟺$A_W\in\{0,M\}$ 영점판정, $d_{\min}$ 반례,
  $\tau$ 지표($\tau{=}0\iff M_2{=}0$, $r{=}0.76$).
- **선행조사 `magic-prior-art.md` §5c**: 부호↔magic 3갈래(증류/하이퍼그래프/Dicke) 점유 정리,
  **(B) 하이퍼그래프 2편 전문 확인**(위상 인코딩·RM(2)/nonquadraticity, support 미취급) →
  "위상 vs support / RM(2) vs RM(1)" 차별화 확정.
- **통합본 `magic-program-overview.md`**: 그동안+앞으로 + M1–M5 로드맵.

---

## 다음 (로드맵 M3~)
과제 A 통계적 근사식 $E[M_2]\approx F(M,n,\tau)$ ★, B 일반-M 닫힌형($A_W$), C Grover T-count,
D 오라클-가림 정식화, E (A)(C) 선행 전문 대조. 상세는 `magic-program-overview.md` Part II–III.
